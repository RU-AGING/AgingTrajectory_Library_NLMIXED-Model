/*===============================================================================
 Program:    ZINB_LatentClass_Trajectories.sas
 Purpose:    Fit Zero-Inflated Negative Binomial (ZINB) latent-class trajectory
             models (growth-mixture) in PROC NLMIXED with user-controlled starting
             values, and produce class curves plus a mixture-mean plot.

 Key features:
   - No carryover of starting values from previous runs.
   - Exactly one PARMS statement (no duplicates).
   - User sets starting values via %let lines (unset values default to 0).
   - Post-estimation parsing of ParameterEstimates to build trajectories and plots.
   - Optional simulator (%sim_data) included to generate toy data.

 Requirements:
   - SAS 9.4+ recommended.
   - Wide data: 1 row per person, repeated outcome columns.
   - For the example call at the bottom: a wide dataset BASE_FILE_SRS with:
       BENE_ID, SUM_Q1 ... SUM_Q12

 Sections:
   [A] Optional simulator (generates sim_long/sim_wide)
   [B] Input prep example (prepare BASE_FILE_SRS)
   [C] Starting values (user edits here)
   [D] Modeling macros (helpers, PARMS builder, likelihood, main fit macro)
   [E] Plotting macro (parses estimates and plots)
   [F] Example run: fit + plots

 Author:     Anum Zafar 
 Last edit:  2026-02-09
===============================================================================*/


/*===============================================================================
 [A] OPTIONAL SIMULATOR: ZINB toy data (sim_long/sim_wide)
    - Generates a mixture of K latent classes with class-specific polynomial mean
      trajectories and class-specific zero-inflation trajectories.
    - Outcome is drawn from a ZINB:
        With prob p(t): structural zero
        Else: NB2(y | mean=mu(t), size=k)   where Var = mu + mu^2/k

 Outputs:
   sim_long: id, class, qtr, y, obs
   sim_wide: Y1..YT (named Y_1..Y_T)  (you can rename in [B])
===============================================================================*/

%macro sim_data(
  class=5,
  n=500,
  T=12,
  seed=2026,
  miss_pattern=balanced,     /* balanced | unbalanced (monotone) */
  p_obs_min=0.6,             /* if unbalanced: minimum fraction observed */
  order=2,                   /* polynomial order for mu(t) used in sim */
  p_order=0                  /* polynomial order for logit(p(t)) used in sim */
);

  data sim_long;
    call streaminit(&seed);

    do id=1 to &n;

      /* 1) latent class membership (equal probs for simulator) */
      class = ceil(rand('uniform') * &class);

      /* 2) follow-up window / monotone missingness */
      if "&miss_pattern"="balanced" then do;
        first_t=1; last_t=&T;
      end;
      else do;
        frac   = rand('uniform')*(1-&p_obs_min) + &p_obs_min;
        n_obs  = ceil(&T*frac);
        first_t=1; last_t=n_obs;
      end;

      /* 3) class-specific parameters for simulation */
      /*    (simple defaults just to create distinct curves) */
      select (class);
        when (1) do; b0= -0.4; b1= 0.05; b2= 0.00; g0= -1.5; g1=0; g2=0; k=1.2; end;
        when (2) do; b0= -0.2; b1= 0.10; b2= 0.00; g0= -1.0; g1=0; g2=0; k=0.8; end;
        when (3) do; b0=  0.0; b1= 0.06; b2= 0.01; g0= -0.8; g1=0; g2=0; k=0.6; end;
        when (4) do; b0= -0.6; b1= 0.16; b2=-0.01; g0= -1.2; g1=0; g2=0; k=1.0; end;
        otherwise do; b0= -0.1; b1= 0.03; b2= 0.02; g0= -0.6; g1=0; g2=0; k=0.9; end;
      end;

      do t=1 to &T;
        qtr=t;
        obs = (t>=first_t and t<=last_t);
        y = .;

        if obs then do;

          /* mean trajectory mu(t) */
          eta = b0 + b1*t
                %if &order>=2 %then + b2*(t*t);
                ;
          mu  = exp(eta);

          /* zero inflation p(t) */
          logitp = g0
                   %if &p_order>=1 %then + g1*t;
                   %if &p_order>=2 %then + g2*(t*t);
                   ;
          p = 1/(1+exp(-logitp));

          /* draw: structural zero vs NB2 */
          u = rand('uniform');
          if u < p then y=0;
          else do;
            /* NB2 draw:
               SAS has RAND('NEGBINOMIAL', p, r) in some versions, but parameterizations vary.
               To avoid ambiguity, we simulate NB via Gamma-Poisson mixture:
                 lambda ~ Gamma(shape=k, scale=mu/k)
                 y ~ Poisson(lambda)
            */
            lambda = rand('gamma', k, mu/k);  /* shape=k, scale=mu/k => mean=mu */
            y = rand('poisson', lambda);
          end;
        end;

        output;
      end;
    end;

    keep id class qtr y obs;
  run;

  /* reshape to wide: Y_1..Y_T */
  proc sort data=sim_long; by id qtr; run;

  proc transpose data=sim_long(where=(obs=1)) out=sim_wide prefix=Y_;
    by id;
    id qtr;
    var y;
  run;

  /*===============================================================================
 [B] INPUT PREP EXAMPLE
    - Expect a wide dataset with ID + repeated outcomes.
    - This block renames sim_wide into BASE_FILE_SRS with SUM_Q1..SUM_Q12 and BENE_ID.
    - Comment out if your data is already prepared.
===============================================================================*/

  data BASE_FILE_SRS;
  set sim_wide;
  rename
    id = BENE_ID
    Y_1 = SUM_Q1  Y_2 = SUM_Q2  Y_3 = SUM_Q3  Y_4 = SUM_Q4  Y_5 = SUM_Q5  Y_6 = SUM_Q6
    Y_7 = SUM_Q7  Y_8 = SUM_Q8  Y_9 = SUM_Q9  Y_10= SUM_Q10 Y_11= SUM_Q11 Y_12= SUM_Q12
  ;
run;

%mend sim_data;

/* Example simulator call (optional) */
/* %sim_data(class=5, n=5000, T=12, seed=1, miss_pattern=balanced); */




/*===============================================================================
 [C] STARTING VALUES (USER EDITS HERE)
    - Define only what you want; all others default to 0.
    - Class 1 is reference for mixing weights (no alpha0_A).
    - New for ZINB: logk_CLASS controls NB size k=exp(logk)>0.
===============================================================================*/

/* Mixing logits (optional) */
%let alpha0_B = -0.50;
%let alpha0_C = -0.80;
/* %let alpha0_D = ; %let alpha0_E = ; */

/* Mean curve coefficients (mu) */
%let beta0_A = -0.40; %let beta1_A = 0.06; %let beta2_A = 0;
%let beta0_B = -0.20; %let beta1_B = 0.10; %let beta2_B = 0;
/* %let beta0_C = ; %let beta1_C = ; %let beta2_C = ; */

/* Zero-inflation logits p(t) */
%let gamma0_A = -1.20;
%let gamma0_B = -0.80;
/* %let gamma0_C = ; */

/* NB size (dispersion): k = exp(logk). Larger k => closer to Poisson. */
%let logk_A  = 0.00;   /* k=1 */
%let logk_B  = 0.18;   /* k~1.20 */
/* %let logk_C = ; */


/*===============================================================================
 [D] MODELING MACROS
===============================================================================*/

/* ---- helper: count words in a space-separated list ---- */
%macro _ct_nwords(list);
  %sysfunc(countw(%superq(list), %str( )))
%mend;

/* ---- helper: stop with message ---- */
%macro _ct_abort(msg);
  %put ERROR: &msg;
  %abort cancel;
%mend;

/* ---- emit one parameter value into PARMS: uses %let or default ---- */
%macro _emit_parm(name, default);
  %if %symexist(&name) %then %do;
    %if %length(%superq(&name)) %then %do;
      &name = %superq(&name)
    %end;
    %else %do;
      &name = &default
    %end;
  %end;
  %else %do;
    &name = &default
  %end;
%mend;

/* arrays */
%macro _ct_array_from_list(name, list, T);
  array &name.[&T] &list.;
%mend;

%macro _ct_declare_mu_arrays(nclass, class_labels, T);
  %local k lab;
  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    array mu_&lab.[&T] _temporary_;
  %end;
%mend;

%macro _ct_declare_pi_arrays(nclass, class_labels, T);
  %local k lab;
  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    array pi_&lab.[&T] _temporary_;
  %end;
%mend;

%macro _ct_declare_zip_arrays(nclass, class_labels, T);
  %local k lab;
  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    array logitp_&lab.[&T] _temporary_;
    array p_&lab.[&T]      _temporary_;
  %end;
%mend;

/* ---- build ONE PARMS statement (alphas + betas + gammas + logk) ---- */
%macro _ct_parms_all(nclass, class_labels, order, p_order);
  parms
  %local k lab;

  /* mixing alphas: class 1 reference => no alpha0_<class1> */
  %do k=2 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    %_emit_parm(%sysfunc(catx(_,alpha0,&lab)), 0)
  %end;

  /* mu(t) polynomial coefficients */
  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    %_emit_parm(%sysfunc(catx(_,beta0,&lab)), 0)
    %_emit_parm(%sysfunc(catx(_,beta1,&lab)), 0)
    %if &order>=2 %then %_emit_parm(%sysfunc(catx(_,beta2,&lab)), 0);
    %if &order>=3 %then %_emit_parm(%sysfunc(catx(_,beta3,&lab)), 0);
  %end;

  /* p(t) = logistic(gamma poly) */
  %if &p_order>=0 %then %do;
    %do k=1 %to &nclass;
      %let lab=%scan(&class_labels, &k, %str( ));
      %_emit_parm(%sysfunc(catx(_,gamma0,&lab)), 0)
      %if &p_order>=1 %then %_emit_parm(%sysfunc(catx(_,gamma1,&lab)), 0);
      %if &p_order>=2 %then %_emit_parm(%sysfunc(catx(_,gamma2,&lab)), 0);
      %if &p_order>=3 %then %_emit_parm(%sysfunc(catx(_,gamma3,&lab)), 0);
    %end;
  %end;

  /* NB size per class: k=exp(logk_*) */
  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    %_emit_parm(%sysfunc(catx(_,logk,&lab)), 0)
  %end;
  ;
%mend;

/* fill mu_k(t) */
%macro _ct_fill_mu_poly(nclass, class_labels, order, T);
  %local k lab;
  do i=1 to &T; t=i;
    %do k=1 %to &nclass;
      %let lab=%scan(&class_labels, &k, %str( ));
      eta_&lab = beta0_&lab + beta1_&lab*t
                 %if &order>=2 %then + beta2_&lab*(t*t);
                 %if &order>=3 %then + beta3_&lab*(t*t*t);
                 ;
      mu_&lab.[i] = exp(eta_&lab);
    %end;
  end;
%mend;

/* fill logitp_k(t) and p_k(t) */
%macro _ct_fill_zip_poly(nclass, class_labels, p_order, T);
  %local k lab;
  do i=1 to &T; t=i;
    %do k=1 %to &nclass;
      %let lab=%scan(&class_labels, &k, %str( ));
      logitp_&lab.[i] = gamma0_&lab
                        %if &p_order>=1 %then + gamma1_&lab*t;
                        %if &p_order>=2 %then + gamma2_&lab*(t*t);
                        %if &p_order>=3 %then + gamma3_&lab*(t*t*t);
                        ;
      p_&lab.[i] = 1/(1+exp(-logitp_&lab.[i]));
    %end;
  end;
%mend;

/* accumulate ZINB per-time loglik per class */
%macro _ct_accumulate_zinb_ll(nclass, class_labels, T);
  %local k lab;
  do i=1 to &T;
    if not missing(Y[i]) then do;

      %do k=1 %to &nclass;
        %let lab=%scan(&class_labels, &k, %str( ));

        k_&lab = exp(logk_&lab);                /* size > 0 */

        /* stable log(p) and log(1-p) */
        lp = -log(1 + exp(-logitp_&lab.[i]));   /* log p */
        lq = -log(1 + exp( logitp_&lab.[i]));   /* log(1-p) */

        /* NB2 log pmf */
        logNB = lgamma(Y[i] + k_&lab) - lgamma(k_&lab) - lgamma(Y[i]+1)
                + k_&lab*(log(k_&lab) - log(k_&lab + mu_&lab.[i]))
                + Y[i]*(log(mu_&lab.[i]) - log(k_&lab + mu_&lab.[i]));

        /* NB mass at zero */
        logNB0 = k_&lab*(log(k_&lab) - log(k_&lab + mu_&lab.[i]));

        if Y[i]=0 then do;
          /* log( p + (1-p)*NB0 ) via log-sum-exp */
          a = lp;
          b = lq + logNB0;
          m0 = max(a,b);
          pi_&lab.[i] = m0 + log(exp(a-m0) + exp(b-m0));
        end;
        else do;
          pi_&lab.[i] = lq + logNB;
        end;

      %end;

    end;
  end;
%mend;

/* mix across classes via softmax weights */
%macro _ct_mixture_ll(nclass, class_labels, T);
  %local k lab;

  den = 1;
  %do k=2 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    den = den + exp(alpha0_&lab);
  %end;

  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    %if &k=1 %then %do; w_&lab = 1/den; %end;
    %else %do;          w_&lab = exp(alpha0_&lab)/den; %end;
  %end;

  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    prod_&lab = 0;
    do i=1 to &T;
      if not missing(pi_&lab.[i]) then prod_&lab + pi_&lab.[i];
    end;
  %end;

  m = prod_%scan(&class_labels, 1, %str( ));
  %do k=2 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    m = max(m, prod_&lab);
  %end;

  sum_exp = 0;
  %do k=1 %to &nclass;
    %let lab=%scan(&class_labels, &k, %str( ));
    sum_exp + w_&lab * exp(prod_&lab - m);
  %end;

  ll = m + log(sum_exp);
%mend;


/* ---- MAIN fit macro: mirrors ct_zip_pois_nlmixed ---- */
%macro ct_zinb_nlmixed(
  data=BASE_FILE_SRS,
  id=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=5,
  class_labels=A B C D E,
  order=2,
  p_order=0,
  T=12,
  tech=newrap,
  maxiter=500,
  pe_out=pe_zinb,
  fit_out=fit_zinb,
  bounds=
);

  /* Safety check: nclass must match number of labels */
  %if %eval(%_ct_nwords(&class_labels) ne &nclass) %then %do;
    %_ct_abort(nclass=&nclass but class_labels=&class_labels has %_ct_nwords(&class_labels) labels. Fix mismatch.);
  %end;

  ods listing;
  ods output ParameterEstimates=&pe_out FitStatistics=&fit_out;

  proc nlmixed data=&data qpoints=1 tech=&tech maxiter=&maxiter;
    %_ct_parms_all(&nclass, &class_labels, &order, &p_order);

    /* optional bounds for stability */
    %if %length(&bounds) %then %do; bounds &bounds; %end;

    %_ct_array_from_list(Y, &yvars, &T);
    %_ct_declare_mu_arrays(&nclass, &class_labels, &T);
    %_ct_declare_pi_arrays(&nclass, &class_labels, &T);
    %_ct_declare_zip_arrays(&nclass, &class_labels, &T);

    %_ct_fill_mu_poly(&nclass, &class_labels, &order, &T);
    %_ct_fill_zip_poly(&nclass, &class_labels, &p_order, &T);

    %_ct_accumulate_zinb_ll(&nclass, &class_labels, &T);
    %_ct_mixture_ll(&nclass, &class_labels, &T);

    one = 1;
    model one ~ general(ll);
    id &id;
  run;

  ods output close;

%mend ct_zinb_nlmixed;


/*===============================================================================
 [E] PLOTTING MACRO (mirrors ct_zip_plots)
    - Reconstructs expected ZINB mean:
        E[Y | class k, t] = (1 - p_k(t)) * mu_k(t)
    - Computes mixture mean using softmax weights from alpha0_*
===============================================================================*/

%macro ct_zinb_plots(pe=pe_zinb, T=12, out_traj=traj_zinb, out_mix=mix_zinb);

  ods listing;

  /* betas */
  proc sql;
    create table _betas as
    select scan(Parameter,2,'_') as class length=32,
           input(compress(substr(Parameter,5),,'kd'), best.) as deg,
           Estimate
    from &pe
    where upcase(substr(Parameter,1,4))='BETA';
  quit;
  proc sort data=_betas; by class deg; run;
  proc transpose data=_betas out=_betas_w prefix=b;
    by class; id deg; var Estimate;
  run;

  /* alphas */
  proc sql;
    create table _alphas as
    select scan(Parameter,2,'_') as class length=32,
           Estimate as alpha
    from &pe
    where upcase(substr(Parameter,1,6))='ALPHA0';
  quit;
  proc sort data=_alphas; by class; run;

  /* gammas */
  proc sql;
    create table _gammas as
    select scan(Parameter,2,'_') as class length=32,
           input(compress(substr(Parameter,7),,'kd'), best.) as deg,
           Estimate
    from &pe
    where upcase(substr(Parameter,1,5))='GAMMA';
  quit;
  proc sort data=_gammas; by class deg; run;
  proc transpose data=_gammas out=_gammas_w prefix=g;
    by class; id deg; var Estimate;
  run;

  /* logk (optional to report; not needed for mean curve) */
  proc sql;
    create table _logk as
    select scan(Parameter,2,'_') as class length=32,
           Estimate as logk
    from &pe
    where upcase(substr(Parameter,1,4))='LOGK';
  quit;
  proc sort data=_logk; by class; run;

  /* combine */
  data _classes;
    merge _betas_w(in=b) _alphas(in=a) _gammas_w(in=g) _logk(in=k);
    by class;
    if missing(alpha) then alpha=0;
    exp_alpha = exp(alpha);
    k_nb = exp(coalesce(logk,0));
  run;

  proc sql noprint;
    select sum(exp_alpha) into :_den from _classes;
  quit;

  data &out_traj;
    set _classes;
    length class $32;
    do t=1 to &T;
      eta    = coalesce(b0,0) + coalesce(b1,0)*t + coalesce(b2,0)*(t*t) + coalesce(b3,0)*(t*t*t);
      mu     = exp(eta);

      logitp = coalesce(g0,0) + coalesce(g1,0)*t + coalesce(g2,0)*(t*t) + coalesce(g3,0)*(t*t*t);
      p      = 1/(1+exp(-logitp));

      w      = exp_alpha / &_den;

      mu_zinb = (1-p)*mu;
      output;
    end;

    keep class t mu p mu_zinb w k_nb;
  run;

  proc sql;
    create table &out_mix as
    select t, sum(w*mu_zinb) as mu_mix_zinb
    from &out_traj
    group by t;
  quit;

  /* plot 1: class curves */
  proc sgplot data=&out_traj;
    series x=t y=mu_zinb / group=class lineattrs=(thickness=2);
    xaxis integer label="Quarter" min=1 max=&T;
    yaxis label="Expected count (ZINB mean)";
    title "Latent-Class ZINB Trajectories";
  run;

  /* plot 2: add mixture mean */
  proc sort data=&out_traj; by t; run;
  data _traj_all;
    merge &out_traj &out_mix;
    by t;
  run;

  proc sgplot data=_traj_all;
    series x=t y=mu_mix_zinb / lineattrs=(pattern=shortdash thickness=3) name="mix" legendlabel="Mixture mean";
    series x=t y=mu_zinb     / group=class lineattrs=(thickness=2);
    keylegend / position=topright;
    xaxis integer label="Quarter" min=1 max=&T;
    yaxis label="Expected count (ZINB mean)";
    title "Latent-Class ZINB Trajectories with Mixture Mean";
  run;

%mend ct_zinb_plots;


/*===============================================================================
 [F] EXAMPLE RUN
===============================================================================*/

%let T=12;

/* (Optional) Generate toy data */
%sim_data(class=5, n=500, T=&T, seed=1, miss_pattern=balanced);

/* Fit model + capture pe/fit in one call */
%ct_zinb_nlmixed(
  data=BASE_FILE_SRS,
  id=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=5,
  class_labels=A B C D E,
  order=2,
  p_order=0,
  T=&T
  /* bounds=%str(-8 < beta0_A beta0_B beta0_C < 8) */
);

/* Plots */
%ct_zinb_plots(pe=pe_zinb, T=&T);


/*===============================================================================
 Optional: Save plots as PNG to WORK (uncomment)

%let outdir=%sysfunc(getoption(work));
ods _all_ close;
ods html5 path="&outdir" (url=none) gpath="&outdir";
ods graphics / reset imagename="ZINB_Trajectories" imagefmt=png;
%ct_zinb_plots(pe=pe_zinb, T=&T);
ods html5 close;

===============================================================================*/
