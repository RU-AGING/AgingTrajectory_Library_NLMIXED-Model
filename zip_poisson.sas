/*Copyright (c) 2026 Community Health and Aging Outcomes (CHAO) Lab, Rutgers University.
Released under the MIT License. Full text in LICENSE at
https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model;*/

/*PROJECT NAME: Traj2 Zero-Inflated Poisson Latent-Class Trajectories 
LAST UPDATED DATE: 24 JUL 2026
DATA SOURCES: NONE. This file defines macros only. INPUT IS the optional simulator OR your own wide TABLE
STATUS: PROTOTYPE, undergoing quality assurance. NOT part of the current Traj2 release
PURPOSE: Zero-Inflated Poisson growth-mixture trajectory models fitted BY PROC NLMIXED with
user-controlled starting values. ZIP separates structural zeros, subjects who would NOT use the service
at all, FROM sampling zeros, subjects who happen NOT TO use it IN a given period. Provides
(a)an optional bounded-count simulator
(b)a fitting MACRO over a common latent-class structure
(c)class-specific MEAN curves AND a mixture-MEAN plot
KEY FEATURES: no carryover of starting values between runs, exactly one PARMS statement with no
duplicates, starting values SET through percent-LET lines with unset values defaulting TO 0, AND
post-estimation parsing of ParameterEstimates TO build the trajectories AND plots.
AUTHOR: Anum Zafar
##########################################################################################################################
*Execution Environment: SAS 9.4 OR later, PROC NLMIXED FROM SAS/STAT. No compiled components         *
*This file defines macros only, apart FROM the demo IN STEP 6. RUN STEP 6 only WHEN you want the demo*
*Comments inside the code-generator macros must stay IN slash-star form. A star-semicolon comment    *
*inside a text-returning MACRO IS emitted INTO the statement it builds AND IS a syntax error         *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 1:Optional simulator. SIM_DATA builds SIM_LONG AND SIM_WIDE
#STEP 2:INPUT prep. RENAME the wide outcomes TO SUM_Q1 TO SUM_QT AND the ID TO BENE_ID
#STEP 3:Starting values. Edit only this block
#STEP 4:Code-generator library AND the fitting MACRO CT_ZIP_POIS_NLMIXED
#STEP 5:Plotting. CT_ZIP_PLOTS draws class curves AND the mixture MEAN
#STEP 6:Example RUN
##########################################################################################################################;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 1: OPTIONAL SIMULATOR
##########################################################################################################################
Builds SIM_LONG and SIM_WIDE, bounded counts on 0 to cap, for the worked example.
##########################################################################################################################*/

%MACRO sim_data(
  dist=poisson,             /* poisson | BIN4 | TPOIS4 | CAT5                 */
class=5,                  /* number of latent classes                         */
N=10000,                  /* subjects                                         */
T=12,                     /* time points (quarters)                           */
seed=2025,                /* RNG seed                                         */
miss_pattern=balanced,    /* balanced | unbalanced (monotone)                 */
p_obs_min=0.6,            /* min fraction observed if unbalanced              */
sigma_re=0.4,             /* subject random effect SD                         */
cap=7                      /* upper bound of support (0..cap)                  */
);

%LOCAL _DIST;
%LET _DIST=%SYSFUNC(UPCASE(%SYSFUNC(STRIP(&dist.))));

DATA sim_long;
CALL STREAMINIT(&seed.);
DO ID = 1to &n.;

/* 1) latent class (equal probs) */
class = CEIL(RAND('uniform') * &class.);

/* 2) follow-up window / monotone missing */
IF "&miss_pattern." = "balanced" THEN DO;
first_t = 1; last_t = &T.;
END;
ELSE DO;
frac   = RAND('uniform')*(1-&p_obs_min.) + &p_obs_min.;
n_obs  = CEIL(&T.*frac);
first_t=1 ; last_t  = n_obs;
END;

/* 3) subject random effect (induces correlation) */
bi = RAND('normal', 0, &sigma_re.);

/* 4) loop over quarters */
DO t = 1to &T.;
qtr = t;
obs = (t >= first_t AND t <= last_t);
y1 = .; y2 = .;

IF obs THEN DO;
/* class-specific mean trajectories (just to diversify) */
SELECT (class);
WHEN (1) DO; mu1 = 0.6+ 0.06*t;            mu2 = 0.8+ 0.03*t; END;
WHEN (2) DO; mu1 = 0.5+ 0.30*t;            mu2 = 0.7+ 0.10*t; END;
WHEN (3) DO; mu1 = 1.0+ 0.10*t + 0.01*t*t; mu2 = 0.9+ 0.08*t; END;
OTHERWISE DO; mu1 = 0.5+ 0.20*t;           mu2 = 0.5+ 0.15*t; END;
END;

/* Map to bounded support \0..cap\ by chosen family */
%IF "&_DIST." = "BIN4" %THEN %DO;
m1 = MAX(0,  mu1*EXP(bi));
m2 = MAX(0,  mu2*EXP(bi));
p1 = MIN(MAX(m1/&cap., 0), 1);
p2 = MIN(MAX(m2/&cap., 0), 1);
y1 = RAND('binomial', p1, &cap.);
y2 = RAND('binomial', p2, &cap.);
%END;
%ELSE %IF "&_DIST." = "TPOIS4" %THEN %DO;
/* Truncated Poisson (0..cap) via renormalized probabilities */
m1 = MAX(1e-6, mu1*EXP(bi));
m2 = MAX(1e-6, mu2*EXP(bi));
ARRAY p1[%EVAL(&cap.+1)];
ARRAY p2[%EVAL(&cap.+1)];
s1 = 0; s2 = 0;
DO k = 0to &cap.; p1[k+1]=PDF('poisson', k, m1); s1 + p1[k+1]; END;
DO k = 0to &cap.; p2[k+1]=PDF('poisson', k, m2); s2 + p2[k+1]; END;
DO k = 0to &cap.; p1[k+1]=p1[k+1]/s1; p2[k+1]=p2[k+1]/s2; END;
y1 = RAND('table', of p1[*]) - 1;
y2 = RAND('table', of p2[*]) - 1;
%END;
%ELSE %IF "&_DIST." = "CAT5" %THEN %DO;
/* Categorical uniform on \0..cap\ */
y1 = RAND('integer', %EVAL(&cap.+1)) - 1;
y2 = RAND('integer', %EVAL(&cap.+1)) - 1;
%END;
%ELSE %DO;
/* Fallback: default to BIN4 */
IF ID=1  AND t=1  THEN PUT "WARNING: dist=&dist not recognized. Defaulting to BIN4.";
m1 = MAX(0,  mu1*EXP(bi));
m2 = MAX(0,  mu2*EXP(bi));
p1 = MIN(MAX(m1/&cap., 0), 1);
p2 = MIN(MAX(m2/&cap., 0), 1);
y1 = RAND('binomial', p1, &cap.);
y2 = RAND('binomial', p2, &cap.);
%END;
END;

OUTPUT;
END;
END;
KEEP ID class qtr y1 y2 obs;
RUN;

/* reshape to wide (quarters as columns) */
PROC SORT DATA=sim_long; BY ID qtr; RUN;

PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=y1_w PREFIX=Y1_;
BY ID; ID qtr; VAR y1;
RUN;

PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=y2_w PREFIX=Y2_;
BY ID; ID qtr; VAR y2;
RUN;

DATA sim_wide;
MERGE y1_w y2_w; BY ID;
RUN;

%MEND sim_data;

/* Demo: build a toy dataset (comment out if you already have real data) */
/* %sim_data(dist=poisson, n=5000, T=12, seed=1); */

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 2: INPUT PREP
##########################################################################################################################
Renames the wide outcomes to SUM_Q1 to SUM_QT and the identifier to BENE_ID, which is the shape
the fitting macro expects. Point it at your own table to use real data instead.
##########################################################################################################################*/

DATA
 BASE_FILE_SRS;
  SET WORK.Y1_W;  /* replace with your source if needed */
RENAME
Y1_1 = SUM_Q1  Y1_2 = SUM_Q2  Y1_3  = SUM_Q3  Y1_4  = SUM_Q4  Y1_5  = SUM_Q5
Y1_6 = SUM_Q6  Y1_7 = SUM_Q7  Y1_8  = SUM_Q8  Y1_9  = SUM_Q9  Y1_10 = SUM_Q10
Y1_11= SUM_Q11 Y1_12= SUM_Q12 ID    = BENE_ID
;
RUN
;

/* Add time index variables (0..11) if you want them for other tooling */
DATA
BASE_FILE_SRS(KEEP=BENE_ID SUM_Q1-SUM_Q12 quar1-quar12);
SET BASE_FILE_SRS;
quar1=0 ; quar2=1 ; quar3=2 ; quar4=3 ; quar5=4 ; quar6=5 ; quar7=6 ; quar8=7 ;
quar9=8 ; quar10=9 ; quar11=10 ; quar12=11 ;
RUN
;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 3: STARTING VALUES, EDIT ONLY THIS BLOCK
##########################################################################################################################
Any parameter left unset defaults to 0. Class A is the mixing reference, so never supply alpha0_A.
##########################################################################################################################*/

/* mixing weights */
%LET alpha0_B = -0.50;
%LET alpha0_C = -0.80;
/* %let alpha0_D = ; %let alpha0_E = ; %let alpha0_F = ; */

/* mean curve betas up to ORDER=2 for each class */
%LET beta0_A =  0.10; %LET beta1_A = 0.02; %LET beta2_A = 0;
%LET beta0_B = -0.20; %LET beta1_B = 0.01; %LET beta2_B = 0;
/* %let beta0_C = ; %let beta1_C = ; %let beta2_C = ; */
/* %let beta0_D = ; %let beta1_D = ; %let beta2_D = ; */
/* %let beta0_E = ; %let beta1_E = ; %let beta2_E = ; */
/* %let beta0_F = ; %let beta1_F = ; %let beta2_F = ; */

/* ZIP logits p(t) up to P_ORDER=0 (gamma0 only here) */
%LET gamma0_A = -1.00;
%LET gamma0_B = -0.50;
/* %let gamma0_C = ; %let gamma0_D = ; %let gamma0_E = ; %let gamma0_F = ; */

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 4: CODE-GENERATOR LIBRARY AND FITTING MACRO
##########################################################################################################################
These macros return TEXT that is spliced into the PROC NLMIXED program. Comments inside them must be
slash-star only. A star-semicolon comment here lands inside the PARMS or model statement being built
and is a syntax error.
##########################################################################################################################*/

/* helper to emit one value per parameter into a single PARMS statement */
%MACRO _emit_parm(name, default);
/* name is a token like alpha0_B */
%IF %SYMEXIST(&name) %THEN %DO;                /* variable exists? */
%IF %LENGTH(%SUPERQ(&name)) %THEN %DO;       /* nonblank value? */
&name = %SUPERQ(&name)
%END;
%ELSE %DO;                                   /* blank -> default */
&name = &default
%END;
%END;
%ELSE %DO;                                     /* does not exist -> default */
&name = &default
%END;
%MEND _emit_parm;

/* array and bookkeeping helpers */
%MACRO _ct_array_from_list(name, list, T); ARRAY &name.[&T] &list.;
%MEND _ct_array_from_list;

%MACRO _ct_declare_mu_arrays(nclass, class_labels, T);
%LOCAL k lab;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
ARRAY mu_&lab.[&T];
%END;
%MEND _ct_declare_mu_arrays;

%MACRO _ct_declare_pi_arrays(nclass, class_labels, T);
%LOCAL k lab;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
ARRAY pi_&lab.[&T];
%END;
%MEND _ct_declare_pi_arrays;

%MACRO _ct_declare_zip_arrays(nclass, class_labels, T);
%LOCAL k lab;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
ARRAY logitp_&lab.[&T];
ARRAY p_&lab.[&T];
%END;
%MEND _ct_declare_zip_arrays;

/* build one PARMS statement using either %let starts or 0 defaults */
%MACRO _ct_parms_all(nclass, class_labels, ORDER, p_order);
PARMS
%LOCAL k lab;
/* Mixing alphas: class 1 is reference -> no alpha0_A */
%DO k=2  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
%_emit_parm(%SYSFUNC(catx(_,alpha0,&lab)), 0)
%END;
/* Mean curve betas */
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
%_emit_parm(%SYSFUNC(catx(_,beta0,&lab)), 0)
%_emit_parm(%SYSFUNC(catx(_,beta1,&lab)), 0)
%IF &order>=2%THEN %_emit_parm(%SYSFUNC(catx(_,beta2,&lab)), 0);
%IF &order>=3%THEN %_emit_parm(%SYSFUNC(catx(_,beta3,&lab)), 0);
%END;
/* ZIP logits p_k(t) */
%IF &p_order>=0%THEN %DO;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
%_emit_parm(%SYSFUNC(catx(_,gamma0,&lab)), 0)
%IF &p_order>=1%THEN %_emit_parm(%SYSFUNC(catx(_,gamma1,&lab)), 0);
%IF &p_order>=2%THEN %_emit_parm(%SYSFUNC(catx(_,gamma2,&lab)), 0);
%IF &p_order>=3%THEN %_emit_parm(%SYSFUNC(catx(_,gamma3,&lab)), 0);
%END;
%END;
;
%MEND _ct_parms_all;

/* fill class-specific mu(t) over T using polynomial of degree order */
%MACRO _ct_fill_mu_poly(nclass, class_labels, ORDER, T);
%LOCAL k lab;
DO i = 1to &T; t = i;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
eta_&lab = beta0_&lab. + beta1_&lab.*t
%IF &order>=2%THEN + beta2_&lab.*(t*t);
%IF &order>=3%THEN + beta3_&lab.*(t*t*t);
;
mu_&lab.[i] = EXP(eta_&lab);
%END;
END;
%MEND _ct_fill_mu_poly;

/* fill class-specific logit p_k(t) over T using polynomial of degree p_order */
%MACRO _ct_fill_zip_poly(nclass, class_labels, p_order, T);
%LOCAL k lab;
DO i = 1to &T; t = i;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
logitp_&lab.[i] = gamma0_&lab.
%IF &p_order>=1%THEN + gamma1_&lab.*t;
%IF &p_order>=2%THEN + gamma2_&lab.*(t*t);
%IF &p_order>=3%THEN + gamma3_&lab.*(t*t*t);
;
p_&lab.[i] = 1/(1+EXP(-logitp_&lab.[i]));
%END;
END;
%MEND _ct_fill_zip_poly;

/* accumulate ZIP log-likelihood per class across observed time points */
%MACRO _ct_accumulate_zip_ll(nclass, class_labels, T);
%LOCAL k lab;
DO i = 1to &T;
IF NOT MISSING(Y[i]) THEN DO;
%DO k=1  %TO &nclass;
%LET lab=%SCAN(&class_labels., &k, %STR( ));
lp = -LOG(1+ EXP(-logitp_&lab.[i]));   /* log p   */
lq = -LOG(1+ EXP( logitp_&lab.[i]));   /* log(1-p) */
IF Y[i]=0then DO;                      /* log-sum-exp for 0 mass */
a = lp; b = lq - mu_&lab.[i]; m = MAX(a,b);
pi_&lab.[i] = m + LOG(EXP(a-m) + EXP(b-m));
END;
ELSE DO;                                /* Poisson mass for Y>0 */
pi_&lab.[i] = lq + (Y[i]*LOG(mu_&lab.[i]) - mu_&lab.[i] - LGAMMA(Y[i]+1));
END;
%END;
END;
END;
%MEND _ct_accumulate_zip_ll;

/* mix across classes via softmax weights from alphas */
%MACRO _ct_mixture_ll_zip(nclass, class_labels, T);
%LOCAL k lab;
den = 1;
%DO k=2  %TO &nclass; %LET lab=%SCAN(&class_labels., &k, %STR( )); den = den + EXP(alpha0_&lab.); %END;

%DO k=1  %TO &nclass; %LET lab=%SCAN(&class_labels., &k, %STR( ));
%IF &k=1  %THEN %DO; w_&lab.=1/den; %END;
%ELSE %DO;          w_&lab.=EXP(alpha0_&lab.)/den; %END;
%END;

%DO k=1  %TO &nclass; %LET lab=%SCAN(&class_labels., &k, %STR( ));
prod_&lab.=0;
DO i=1  TO &T; IF NOT MISSING(pi_&lab.[i]) THEN prod_&lab. + pi_&lab.[i]; END;
%END;

m = prod_%SCAN(&class_labels.,1,%STR( ));
%DO k=2  %TO &nclass; %LET lab=%SCAN(&class_labels., &k, %STR( )); m = MAX(m, prod_&lab.); %END;

sum_exp = 0;
%DO k=1  %TO &nclass; %LET lab=%SCAN(&class_labels., &k, %STR( ));
sum_exp + w_&lab.*EXP(prod_&lab. - m);
%END;

ll = m + LOG(sum_exp);
%MEND _ct_mixture_ll_zip;

/* main modeling macro: runs NLMIXED and captures ParameterEstimates and FitStatistics */
%MACRO ct_zip_pois_nlmixed(
DATA=BASE_FILE_SRS,
ID=BENE_ID,
yvars=SUM_Q1-SUM_Q12,
nclass=6,
class_labels=A B C D E F,
ORDER=2,
p_order=0,
T=12,
tech=newrap,
maxiter=500,
pe_out=pe_zip,
fit_out=fit_zip,
BOUNDS=
);
ODS LISTING;
ODS OUTPUT ParameterEstimates=&pe_out FitStatistics=&fit_out;

PROC NLMIXED DATA=&data qpoints=1  tech=&tech maxiter=&maxiter;
/* Build single PARMS with user starts or 0 defaults */
%_ct_parms_all(&nclass, &class_labels, &order, &p_order);

/* Optional stability bounds (example: bounds %str(-6 < beta0_A < 6)) */
%IF %LENGTH(&bounds) %THEN %DO; BOUNDS &bounds; %END;

/* Build arrays and deterministic pieces */
%_ct_array_from_list(Y, &yvars, &T);
%_ct_declare_mu_arrays(&nclass, &class_labels, &T);
%_ct_declare_pi_arrays(&nclass, &class_labels, &T);
%_ct_declare_zip_arrays(&nclass, &class_labels, &T);

%_ct_fill_mu_poly(&nclass, &class_labels, &order, &T);
%_ct_fill_zip_poly(&nclass, &class_labels, &p_order, &T);

/* Accumulate per-class log-likelihood and mix */
%_ct_accumulate_zip_ll(&nclass, &class_labels, &T);
%_ct_mixture_ll_zip(&nclass, &class_labels, &T);

one = 1;
MODEL one ~ GENERAL(ll);
ID &id;
RUN;

ODS OUTPUT CLOSE;
%MEND ct_zip_pois_nlmixed;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 5: PLOTTING
##########################################################################################################################
Reconstructs the class mean curves from ParameterEstimates and overlays the mixture mean.
##########################################################################################################################*/

%MACRO ct_zip_plots(pe=pe_zip, T=12, out_traj=traj, out_mix=mix);
  ODS LISTING;

  /* extract class-specific betas */
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

/* mixing weights (alpha logits) */
PROC SQL;
CREATE TABLE _alphas AS
SELECT SCAN(Parameter,2,'_') AS class LENGTH=32,
ESTIMATE AS alpha
FROM &pe
WHERE UPCASE(SUBSTR(Parameter,1,6))='ALPHA0';
QUIT;
PROC SORT DATA=_alphas; BY class; RUN;

/* ZIP logits p(t) coefficients (gamma) */
PROC SQL;
CREATE TABLE _gammas AS
SELECT SCAN(Parameter,2,'_') AS class LENGTH=32,
INPUT(COMPRESS(SUBSTR(Parameter,7),,'kd'), best.) AS deg,
ESTIMATE
FROM &pe
WHERE UPCASE(SUBSTR(Parameter,1,5))='GAMMA';
QUIT;
PROC SORT DATA=_gammas; BY class deg; RUN;
PROC TRANSPOSE DATA=_gammas OUT=_gammas_w PREFIX=g;
BY class; ID deg; VAR ESTIMATE;
RUN;

/* combine class parameters; alpha missing -> 0 for reference class */
DATA _classes;
MERGE _betas_w(IN=b) _alphas(IN=a) _gammas_w(IN=g);
BY class;
IF MISSING(alpha) THEN alpha=0 ;
exp_alpha = EXP(alpha);
RUN;

/* denominator for class weights */
PROC SQL NOPRINT;
SELECT SUM(exp_alpha) INTO :_den FROM _classes;
QUIT;

/* build trajectories over time */
DATA &out_traj;
SET _classes;
LENGTH class $32;
DO t = 1to &T;
eta    = COALESCE(b0,0) + COALESCE(b1,0)*t + COALESCE(b2,0)*(t*t) + COALESCE(b3,0)*(t*t*t);
mu     = EXP(eta);
logitp = COALESCE(g0,0) + COALESCE(g1,0)*t + COALESCE(g2,0)*(t*t) + COALESCE(g3,0)*(t*t*t);
p      = 1/(1+EXP(-logitp));
w      = exp_alpha / &_den;
mu_zip = (1-p)*mu;
OUTPUT;
END;
KEEP class t mu p mu_zip w;
RUN;

/* mixture mean across classes */
PROC SQL;
CREATE TABLE &out_mix AS
SELECT t, SUM(w*mu_zip) AS mu_mix_zip
FROM &out_traj
GROUP BY t;
QUIT;

/* Plot 1: class ZIP means */
PROC SGPLOT DATA=&out_traj;
SERIES x=t y=mu_zip / GROUP=class lineattrs=(thickness=2);
XAXIS integer LABEL="Quarter" MIN=1  MAX=&T;
YAXIS LABEL="Expected count (ZIP mean)";
TITLE "Latent-Class ZIP Trajectories";
RUN;

/* Plot 2: class curves + mixture mean */
PROC SORT DATA=&out_traj; BY t; RUN;
DATA traj_all; MERGE &out_traj &out_mix; BY t; RUN;

PROC SGPLOT DATA=traj_all;
SERIES x=t y=mu_mix_zip / lineattrs=(pattern=shortdash thickness=3) name="mix" legendlabel="Mixture mean";
SERIES x=t y=mu_zip     / GROUP=class lineattrs=(thickness=2);
KEYLEGEND / position=topright;
XAXIS integer LABEL="Quarter" MIN=1  MAX=&T;
YAXIS LABEL="Expected count (ZIP mean)";
TITLE "Latent-Class ZIP Trajectories with Mixture Mean";
RUN;
%MEND ct_zip_plots;

/*##########################################################################################################################
/*##########################################################################################################################
*STEP 6: EXAMPLE RUN
##########################################################################################################################
Executes on submit. Comment this block out to use the file as a macro library only.
##########################################################################################################################*/

%LET T=12;

%ct_zip_pois_nlmixed(
  DATA=BASE_FILE_SRS,
  ID=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=6,
  class_labels=A B C D E F,
  ORDER=2,
  p_order=0,
  T=&T
  /* bounds=%str(-6 < alpha0_B alpha0_C < 6) */
);

/* Produce the plots from ParameterEstimates (default pe_zip) */
%ct_zip_plots(pe=pe_zip, T=&T);

/*##########################################################################################################################
 Optional: Save plots as PNG to WORK and embed in Results (uncomment to use)

%let outdir=%sysfunc(getoption(work));
ods _all_ close;
ods html5 path="&outdir" (url=none) gpath="&outdir";
ods graphics / reset imagename="ZIP_Trajectories" imagefmt=png;
%ct_zip_plots(pe=pe_zip, T=&T);
ods html5 close;

##########################################################################################################################*/

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT DATASETS
sim_long, sim_wide      simulated data from SIM_DATA
BASE_FILE_SRS           the fitting contract, wide format with SUM_Q1 to SUM_QT and BENE_ID
pe_zip                  parameter estimates from CT_ZIP_POIS_NLMIXED
fit_zip                 fit statistics including AIC and BIC for choosing the class count
traj, mix               class curves and mixture mean from CT_ZIP_PLOTS

NOTE ON STATUS
  This is a prototype undergoing quality assurance and is not part of the current Traj2 release.
  Interfaces and behaviour may change. The released outcome families are ordinal-probit and
  censored-normal continuous.
##########################################################################################################################*/
