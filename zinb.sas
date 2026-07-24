/*Copyright (c) 2026 Community Health and Aging Outcomes (CHAO) Lab, Rutgers University.
Released under the MIT License. Full text in LICENSE at
https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model;*/

*PROJECT NAME: Traj2 Zero-Inflated Negative Binomial Latent-Class Trajectories
/*LAST UPDATED DATE: 24 JUL 2026
DATA SOURCES: NONE. This file defines macros only. INPUT IS the optional simulator OR your own wide TABLE
STATUS: PROTOTYPE, undergoing quality assurance. NOT part of the current Traj2 release
PURPOSE: Zero-Inflated Negative Binomial growth-mixture trajectory models fitted BY PROC NLMIXED, with
user-controlled starting values. ZINB extends ZIP BY replacing the Poisson component with an NB2
distribution, so the variance IS mu + mu squared over k rather than mu, which accommodates the
overdispersion routinely seen IN Medicare claims. Provides
(a)an optional ZINB simulator
(b)a fitting MACRO over a common latent-class structure
(c)class-specific MEAN curves AND a mixture-MEAN plot
AUTHOR: Anum Zafar
##########################################################################################################################
*Execution Environment: SAS 9.4 OR later, PROC NLMIXED FROM SAS/STAT. No compiled components         *
*This file defines macros only, apart FROM the demo IN STEP 5. RUN STEP 5 only WHEN you want the demo*
*Comments inside the code-generator macros must stay IN slash-star form. A star-semicolon comment    *
*inside a text-returning MACRO IS emitted INTO the statement it builds AND IS a syntax error         *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 1:Optional simulator. SIM_DATA builds SIM_LONG, SIM_WIDE AND BASE_FILE_SRS
#STEP 2:Starting values. Edit only this block
#STEP 3:Code-generator library AND the fitting MACRO CT_ZINB_NLMIXED
#STEP 4:Plotting. CT_ZINB_PLOTS draws class curves AND the mixture MEAN
#STEP 5:Example RUN, a three-class demo
##########################################################################################################################;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 1: OPTIONAL SIMULATOR
##########################################################################################################################
Builds SIM_LONG, SIM_WIDE and BASE_FILE_SRS for the worked example. Structural zeros are drawn with
probability p(t), otherwise a Gamma-Poisson mixture supplies the count, so the simulated data carry the
same overdispersion the model is meant to absorb.
##########################################################################################################################*/
%MACRO sim_data(
  class=5, N=500, T=12, seed=2026,
  miss_pattern=balanced, p_obs_min=0.6, ORDER=2, p_order=0
);

  DATA sim_long;
    CALL STREAMINIT(&seed);
    DO ID=1 TO &n;
      class = CEIL(RAND('uniform') * &class);

      IF "&miss_pattern"="balanced" THEN DO; first_t=1; last_t=&T; END;
      ELSE DO;
        frac=RAND('uniform')*(1-&p_obs_min) + &p_obs_min;
        n_obs=CEIL(&T*frac); first_t=1; last_t=n_obs;
      END;

      SELECT (class);
        WHEN (1) DO; b0=-0.4; b1=0.05; b2= 0.00; g0=-1.5; g1=0; g2=0; k=1.2; END;
        WHEN (2) DO; b0=-0.2; b1=0.10; b2= 0.00; g0=-1.0; g1=0; g2=0; k=0.8; END;
        WHEN (3) DO; b0= 0.0; b1=0.06; b2= 0.01; g0=-0.8; g1=0; g2=0; k=0.6; END;
        WHEN (4) DO; b0=-0.6; b1=0.16; b2=-0.01; g0=-1.2; g1=0; g2=0; k=1.0; END;
        OTHERWISE DO; b0=-0.1; b1=0.03; b2=0.02; g0=-0.6; g1=0; g2=0; k=0.9; END;
      END;

      DO t=1 TO &T;
        qtr=t; obs=(t>=first_t AND t<=last_t); y=.;
        IF obs THEN DO;
          eta = b0 + b1*t %IF &order>=2 %THEN + b2*(t*t); ;
          mu  = EXP(eta);

          logitp = g0
                   %IF &p_order>=1 %THEN + g1*t;
                   %IF &p_order>=2 %THEN + g2*(t*t);
                   ;
          p = 1/(1+EXP(-logitp));

          u = RAND('uniform');
          IF u < p THEN y=0;
          ELSE DO;
            lambda = RAND('gamma', k, mu/k);  /* shape=k, scale=mu/k => mean=mu */
y = RAND('poisson', lambda);
END;
END;
OUTPUT;
END;
END;
KEEP ID class qtr y obs;
RUN;

PROC SORT DATA=sim_long; BY ID qtr; RUN;

PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=sim_wide PREFIX=Y_;
BY ID; ID qtr; VAR y;
RUN;

/* [B] Build BASE_FILE_SRS with SUM_Q1..SUM_Q12, BENE_ID, and quar1..quarT=1..T */
DATA BASE_FILE_SRS;
SET sim_wide;
RENAME
ID = BENE_ID
Y_1 = SUM_Q1  Y_2 = SUM_Q2  Y_3 = SUM_Q3  Y_4 = SUM_Q4  Y_5 = SUM_Q5  Y_6 = SUM_Q6
Y_7 = SUM_Q7  Y_8 = SUM_Q8  Y_9 = SUM_Q9  Y_10= SUM_Q10 Y_11= SUM_Q11 Y_12= SUM_Q12
;
RUN;

DATA BASE_FILE_SRS;
SET BASE_FILE_SRS;
ARRAY quar[&T] quar1-quar&T;
DO _i=1 TO &T; quar[_i]=_i; END;
DROP _i;
RUN;

%MEND sim_data;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 2: STARTING VALUES, EDIT ONLY THIS BLOCK
##########################################################################################################################
Any parameter left unset defaults to 0. logk starts at 0, meaning k = 1 and moderate dispersion.
##########################################################################################################################*/
%LET alpha0_B = -0.50;
%LET alpha0_C = -0.80;

%LET beta0_A = -0.40; %LET beta1_A = 0.06; %LET beta2_A = 0;
%LET beta0_B = -0.20; %LET beta1_B = 0.10; %LET beta2_B = 0;

%LET gamma0_A = -1.20;
%LET gamma0_B = -0.80;

%LET logk_A  = 0.00;
%LET logk_B  = 0.18;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 3: CODE-GENERATOR LIBRARY AND FITTING MACRO
##########################################################################################################################
These macros return TEXT that is spliced into the PROC NLMIXED program. Comments inside them must be
slash-star only. A star-semicolon comment here lands inside the PARMS or model statement being built
and is a syntax error.
##########################################################################################################################*/
%MACRO _ct_nwords(list);
  %SYSFUNC(COUNTW(%SUPERQ(list), %STR( )))
%MEND _ct_nwords;

%MACRO _ct_abort(msg);
  %PUT ERROR: &msg;
  %ABORT cancel;
%MEND _ct_abort;

%MACRO _emit_parm(name, default);
  %IF %SYMEXIST(&name) %THEN %DO;
    %IF %LENGTH(%SUPERQ(&name)) %THEN %DO;
      &name = %SUPERQ(&name)
    %END;
    %ELSE %DO;
      &name = &default
    %END;
  %END;
  %ELSE %DO;
    &name = &default
  %END;
%MEND _emit_parm;

%MACRO _ct_array_from_list(name, list, T);
  ARRAY &name.[&T] &list.;
%MEND _ct_array_from_list;

%MACRO _ct_declare_mu_arrays(nclass, class_labels, T);
  %LOCAL k lab;
  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    ARRAY mu_&lab.[&T] _temporary_;
  %END;
%MEND _ct_declare_mu_arrays;

%MACRO _ct_declare_pi_arrays(nclass, class_labels, T);
  %LOCAL k lab;
  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    ARRAY pi_&lab.[&T] _temporary_;
  %END;
%MEND _ct_declare_pi_arrays;

%MACRO _ct_declare_zip_arrays(nclass, class_labels, T);
  %LOCAL k lab;
  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    ARRAY logitp_&lab.[&T] _temporary_;
    ARRAY p_&lab.[&T]      _temporary_;
  %END;
%MEND _ct_declare_zip_arrays;

%MACRO _ct_parms_all(nclass, class_labels, ORDER, p_order);
  PARMS
  %LOCAL k lab;

  %DO k=2 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    %_emit_parm(%SYSFUNC(catx(_,alpha0,&lab)), 0)
  %END;

  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    %_emit_parm(%SYSFUNC(catx(_,beta0,&lab)), 0)
    %_emit_parm(%SYSFUNC(catx(_,beta1,&lab)), 0)
    %IF &order>=2 %THEN %_emit_parm(%SYSFUNC(catx(_,beta2,&lab)), 0);
    %IF &order>=3 %THEN %_emit_parm(%SYSFUNC(catx(_,beta3,&lab)), 0);
  %END;

  %IF &p_order>=0 %THEN %DO;
    %DO k=1 %TO &nclass;
      %LET lab=%SCAN(&class_labels, &k, %STR( ));
      %_emit_parm(%SYSFUNC(catx(_,gamma0,&lab)), 0)
      %IF &p_order>=1 %THEN %_emit_parm(%SYSFUNC(catx(_,gamma1,&lab)), 0);
      %IF &p_order>=2 %THEN %_emit_parm(%SYSFUNC(catx(_,gamma2,&lab)), 0);
      %IF &p_order>=3 %THEN %_emit_parm(%SYSFUNC(catx(_,gamma3,&lab)), 0);
    %END;
  %END;

  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    %_emit_parm(%SYSFUNC(catx(_,logk,&lab)), 0)
  %END;
  ;
%MEND _ct_parms_all;

%MACRO _ct_fill_mu_poly(nclass, class_labels, ORDER, T);
  %LOCAL k lab;
  DO i=1 TO &T; t=i;
    %DO k=1 %TO &nclass;
      %LET lab=%SCAN(&class_labels, &k, %STR( ));
      eta_&lab = beta0_&lab + beta1_&lab*t
                 %IF &order>=2 %THEN + beta2_&lab*(t*t);
                 %IF &order>=3 %THEN + beta3_&lab*(t*t*t);
                 ;
      mu_&lab.[i] = EXP(eta_&lab);
    %END;
  END;
%MEND _ct_fill_mu_poly;

%MACRO _ct_fill_zip_poly(nclass, class_labels, p_order, T);
  %LOCAL k lab;
  DO i=1 TO &T; t=i;
    %DO k=1 %TO &nclass;
      %LET lab=%SCAN(&class_labels, &k, %STR( ));
      logitp_&lab.[i] = gamma0_&lab
                        %IF &p_order>=1 %THEN + gamma1_&lab*t;
                        %IF &p_order>=2 %THEN + gamma2_&lab*(t*t);
                        %IF &p_order>=3 %THEN + gamma3_&lab*(t*t*t);
                        ;
      p_&lab.[i] = 1/(1+EXP(-logitp_&lab.[i]));
    %END;
  END;
%MEND _ct_fill_zip_poly;

%MACRO _ct_accumulate_zinb_ll(nclass, class_labels, T);
  %LOCAL k lab;
  DO i=1 TO &T;
    IF NOT MISSING(Y[i]) THEN DO;

      %DO k=1 %TO &nclass;
        %LET lab=%SCAN(&class_labels, &k, %STR( ));

        k_&lab = EXP(logk_&lab);

        lp = -LOG(1 + EXP(-logitp_&lab.[i]));
        lq = -LOG(1 + EXP( logitp_&lab.[i]));

        logNB = LGAMMA(Y[i] + k_&lab) - LGAMMA(k_&lab) - LGAMMA(Y[i]+1)
                + k_&lab*(LOG(k_&lab) - LOG(k_&lab + mu_&lab.[i]))
                + Y[i]*(LOG(mu_&lab.[i]) - LOG(k_&lab + mu_&lab.[i]));

        logNB0 = k_&lab*(LOG(k_&lab) - LOG(k_&lab + mu_&lab.[i]));

        IF Y[i]=0 THEN DO;
          a = lp;
          b = lq + logNB0;
          m0 = MAX(a,b);
          pi_&lab.[i] = m0 + LOG(EXP(a-m0) + EXP(b-m0));
        END;
        ELSE DO;
          pi_&lab.[i] = lq + logNB;
        END;

      %END;

    END;
  END;
%MEND _ct_accumulate_zinb_ll;

%MACRO _ct_mixture_ll(nclass, class_labels, T);
  %LOCAL k lab;

  den = 1;
  %DO k=2 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    den = den + EXP(alpha0_&lab);
  %END;

  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    %IF &k=1 %THEN %DO; w_&lab = 1/den; %END;
    %ELSE %DO;          w_&lab = EXP(alpha0_&lab)/den; %END;
  %END;

  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    prod_&lab = 0;
    DO i=1 TO &T;
      IF NOT MISSING(pi_&lab.[i]) THEN prod_&lab + pi_&lab.[i];
    END;
  %END;

  m = prod_%SCAN(&class_labels, 1, %STR( ));
  %DO k=2 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    m = MAX(m, prod_&lab);
  %END;

  sum_exp = 0;
  %DO k=1 %TO &nclass;
    %LET lab=%SCAN(&class_labels, &k, %STR( ));
    sum_exp + w_&lab * EXP(prod_&lab - m);
  %END;

  ll = m + LOG(sum_exp);
%MEND _ct_mixture_ll;

%MACRO ct_zinb_nlmixed(
  DATA=BASE_FILE_SRS,
  ID=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=5,
  class_labels=A B C D E,
  ORDER=2,
  p_order=0,
  T=12,
  tech=newrap,
  maxiter=500,
  pe_out=pe_zinb,
  fit_out=fit_zinb,
  BOUNDS=
);

  %IF %EVAL(%_ct_nwords(&class_labels) NE &nclass) %THEN %DO;
    %_ct_abort(nclass=&nclass but class_labels=&class_labels has %_ct_nwords(&class_labels) labels. Fix mismatch.);
  %END;

  ODS LISTING;
  ODS OUTPUT ParameterEstimates=&pe_out FitStatistics=&fit_out;

  PROC NLMIXED DATA=&data qpoints=1 tech=&tech maxiter=&maxiter;
    %_ct_parms_all(&nclass, &class_labels, &order, &p_order);

    %IF %LENGTH(&bounds) %THEN %DO; BOUNDS &bounds; %END;

    %_ct_array_from_list(Y, &yvars, &T);
    %_ct_declare_mu_arrays(&nclass, &class_labels, &T);
    %_ct_declare_pi_arrays(&nclass, &class_labels, &T);
    %_ct_declare_zip_arrays(&nclass, &class_labels, &T);

    %_ct_fill_mu_poly(&nclass, &class_labels, &order, &T);
    %_ct_fill_zip_poly(&nclass, &class_labels, &p_order, &T);

    %_ct_accumulate_zinb_ll(&nclass, &class_labels, &T);
    %_ct_mixture_ll(&nclass, &class_labels, &T);

    one = 1;
    MODEL one ~ GENERAL(ll);
    ID &id;
  RUN;

  ODS OUTPUT CLOSE;

%MEND ct_zinb_nlmixed;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 4: PLOTTING
##########################################################################################################################
Reconstructs E[Y | class, t] = (1 - p(t)) * mu(t) per class, then overlays the mixture mean.
##########################################################################################################################*/
%MACRO ct_zinb_plots(pe=pe_zinb, T=12, out_traj=traj_zinb, out_mix=mix_zinb);

  ODS LISTING;

  PROC SQL;
    CREATE TABLE _betas AS
    SELECT SCAN(Parameter,2,'_') AS class LENGTH=32,
           INPUT(COMPRESS(SUBSTR(Parameter,5),,'kd'), best.) AS deg,
           ESTIMATE
    FROM &pe
    WHERE UPCASE(SUBSTR(Parameter,1,4))='BETA';
  QUIT;
  PROC SORT DATA=_betas; BY class deg; RUN;
  PROC TRANSPOSE DATA=_betas OUT=_betas_w PREFIX=b;
    BY class; ID deg; VAR ESTIMATE;
  RUN;

  PROC SQL;
    CREATE TABLE _alphas AS
    SELECT SCAN(Parameter,2,'_') AS class LENGTH=32,
           ESTIMATE AS alpha
    FROM &pe
    WHERE UPCASE(SUBSTR(Parameter,1,6))='ALPHA0';
  QUIT;
  PROC SORT DATA=_alphas; BY class; RUN;

  /* FIX: gamma degree digit sits at position 6 ("gamma" is 5 chars), not 7 */
PROC SQL;
CREATE TABLE _gammas AS
SELECT SCAN(Parameter,2,'_') AS class LENGTH=32,
INPUT(COMPRESS(SUBSTR(Parameter,6),,'kd'), best.) AS deg,
ESTIMATE
FROM &pe
WHERE UPCASE(SUBSTR(Parameter,1,5))='GAMMA';
QUIT;
PROC SORT DATA=_gammas; BY class deg; RUN;
PROC TRANSPOSE DATA=_gammas OUT=_gammas_w PREFIX=g;
BY class; ID deg; VAR ESTIMATE;
RUN;

PROC SQL;
CREATE TABLE _logk AS
SELECT SCAN(Parameter,2,'_') AS class LENGTH=32,
ESTIMATE AS logk
FROM &pe
WHERE UPCASE(SUBSTR(Parameter,1,4))='LOGK';
QUIT;
PROC SORT DATA=_logk; BY class; RUN;

DATA _classes;
MERGE _betas_w(IN=b) _alphas(IN=a) _gammas_w(IN=g) _logk(IN=k);
BY class;
IF MISSING(alpha) THEN alpha=0;
exp_alpha = EXP(alpha);
k_nb = EXP(COALESCE(logk,0));
RUN;

PROC SQL NOPRINT;
SELECT SUM(exp_alpha) INTO :_den FROM _classes;
QUIT;

DATA &out_traj;
SET _classes;
LENGTH class $32;
DO t=1 TO &T;
eta    = COALESCE(b0,0) + COALESCE(b1,0)*t + COALESCE(b2,0)*(t*t) + COALESCE(b3,0)*(t*t*t);
mu     = EXP(eta);

logitp = COALESCE(g0,0) + COALESCE(g1,0)*t + COALESCE(g2,0)*(t*t) + COALESCE(g3,0)*(t*t*t);
p      = 1/(1+EXP(-logitp));

w      = exp_alpha / &_den;

mu_zinb = (1-p)*mu;
OUTPUT;
END;

KEEP class t mu p mu_zinb w k_nb;
RUN;

PROC SQL;
CREATE TABLE &out_mix AS
SELECT t, SUM(w*mu_zinb) AS mu_mix_zinb
FROM &out_traj
GROUP BY t;
QUIT;

PROC SGPLOT DATA=&out_traj;
SERIES x=t y=mu_zinb / GROUP=class lineattrs=(thickness=2);
XAXIS integer LABEL="Quarter" MIN=1 MAX=&T;
YAXIS LABEL="Expected count (ZINB mean)";
TITLE "Latent-Class ZINB Trajectories";
RUN;

PROC SORT DATA=&out_traj; BY t; RUN;
DATA _traj_all;
MERGE &out_traj &out_mix;
BY t;
RUN;

PROC SGPLOT DATA=_traj_all;
SERIES x=t y=mu_mix_zinb / lineattrs=(pattern=shortdash thickness=3) name="mix" legendlabel="Mixture mean";
SERIES x=t y=mu_zinb     / GROUP=class lineattrs=(thickness=2);
KEYLEGEND / position=topright;
XAXIS integer LABEL="Quarter" MIN=1 MAX=&T;
YAXIS LABEL="Expected count (ZINB mean)";
TITLE "Latent-Class ZINB Trajectories with Mixture Mean";
RUN;

%MEND ct_zinb_plots;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 5: EXAMPLE RUN, THREE-CLASS DEMO
##########################################################################################################################
Executes on submit. Comment this block out to use the file as a macro library only.
##########################################################################################################################*/
%LET T=12;

%sim_data(class=3, N=500, T=&T, seed=1, miss_pattern=balanced);

%ct_zinb_nlmixed(
  DATA=BASE_FILE_SRS,
  ID=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=3,
  class_labels=A B C,
  ORDER=2,
  p_order=0,
  T=&T
);

%ct_zinb_plots(pe=pe_zinb, T=&T);

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT DATASETS
sim_long, sim_wide      simulated data from SIM_DATA
BASE_FILE_SRS           the fitting contract, wide format with SUM_Q1 to SUM_QT and BENE_ID
pe_zinb                 parameter estimates from CT_ZINB_NLMIXED
fit_zinb                fit statistics including AIC and BIC for choosing the class count
traj_zinb, mix_zinb     class curves and mixture mean from CT_ZINB_PLOTS

NOTE ON STATUS
  This is a prototype undergoing quality assurance and is not part of the current Traj2 release.
  Interfaces and behaviour may change. The released outcome families are ordinal-probit and
  censored-normal continuous.
##########################################################################################################################*/
