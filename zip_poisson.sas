/*===============================================================================
 Program:    ZIP_LatentClass_Trajectories.sas
 Purpose:    Fit Zero-Inflated Poisson (ZIP) latent-class trajectory models
             (growth-mixture) in PROC NLMIXED with user-controlled starting values,
             and produce class curves plus a mixture-mean plot.

 Key features:
   - No carryover of starting values from previous runs.
   - Exactly one PARMS statement (no duplicates).
   - User sets starting values via %let lines (unset values default to 0).
   - Post-estimation parsing of ParameterEstimates to build trajectories and plots.

 Requirements:
   - SAS 9.4+ recommended.
   - Input long or wide data of integer counts per time point (quarters in this example).
   - For the example call at the bottom: a wide dataset named BASE_FILE_SRS with:
       BENE_ID, SUM_Q1 ... SUM_Q12

 Optional:
   - A self-contained simulator (%sim_data) is included to generate toy data.
   - ODS HTML5 lines (commented) show how to save plots to PNGs.

 Sections:
   [A] Optional simulator (bounded outcomes 0..cap)
   [B] Input prep example (reshape/rename to BASE_FILE_SRS)
   [C] Starting values (user edits here)
   [D] Modeling macros (helpers, PARMS builder, likelihood, main fit macro)
   [E] Plotting macro (parses estimates and plots)
   [F] Example run: fit + plots

 Author:    Anum Zafar
 Last edit:  02/09/2026
===============================================================================*/

/*===============================================================================
 [A] OPTIONAL: Bounded (0..cap) Simulator to create toy data (sim_long/sim_wide)
    - dist=poisson (default, truncated to 0..cap), BIN4 (binomial), TPOIS4 (trunc. Poisson),
      CAT5 (categorical on 0..cap)
    - Outputs:
        sim_long(id, class, qtr, y1, y2, obs)
        sim_wide(Y1_1..Y1_&T, Y2_1..Y2_&T)
===============================================================================*/

%macro sim_data(
  dist=poisson,             /* poisson | BIN4 | TPOIS4 | CAT5                 */
  class=5,                  /* number of latent classes                         */
  n=10000,                  /* subjects                                         */
  T=12,                     /* time points (quarters)                           */
  seed=2025,                /* RNG seed                                         */
  miss_pattern=balanced,    /* balanced | unbalanced (monotone)                 */
  p_obs_min=
0.6
,            /* min fraction observed if unbalanced              */
  sigma_re=
0.4
,             /* subject random effect SD                         */
  cap=7                      /* upper bound of support (0..cap)                  */
);

  %local _DIST;
  %let _DIST=%sysfunc(upcase(%sysfunc(strip(&dist.))));

  data sim_long;
    call streaminit(&seed.);
    do id = 
1
 to &n.;

      /* 1) latent class (equal probs) */
      class = ceil(rand('uniform') * &class.);

      /* 2) follow-up window / monotone missing */
      if "&miss_pattern." = "balanced" then do;
        first_t = 
1
; last_t = &T.;
      end;
      else do;
        frac   = rand('uniform')*(
1
-&p_obs_min.) + &p_obs_min.;
        n_obs  = ceil(&T.*frac);
        first_t=1 ; last_t  = n_obs;
      end;

      /* 3) subject random effect (induces correlation) */
      bi = rand('normal', 
0
, &sigma_re.);

      /* 4) loop over quarters */
      do t = 
1
 to &T.;
        qtr = t;
        obs = (t >= first_t and t <= last_t);
        y1 = 
.
; y2 = 
.
;

        if obs then do;
          /* class-specific mean trajectories (just to diversify) */
          select (class);
            when (
1
) do; mu1 = 
0.6
 + 
0.06
*t;            mu2 = 
0.8
 + 
0.03
*t; end;
            when (
2
) do; mu1 = 
0.5
 + 
0.30
*t;            mu2 = 
0.7
 + 
0.10
*t; end;
            when (
3
) do; mu1 = 
1.0
 + 
0.10
*t + 
0.01
*t*t; mu2 = 
0.9
 + 
0.08
*t; end;
            otherwise do; mu1 = 
0.5
 + 
0.20
*t;           mu2 = 
0.5
 + 
0.15
*t; end;
          end;

          /* Map to bounded support \0..cap\ by chosen family */
          %if "&_DIST." = "BIN4" %then %do;
            m1 = max(
0
,  mu1*exp(bi));
            m2 = max(
0
,  mu2*exp(bi));
            p1 = min(max(m1/&cap., 
0
), 
1
);
            p2 = min(max(m2/&cap., 
0
), 
1
);
            y1 = rand('binomial', p1, &cap.);
            y2 = rand('binomial', p2, &cap.);
          %end;
          %else %if "&_DIST." = "TPOIS4" %then %do;
            /* Truncated Poisson (0..cap) via renormalized probabilities */
            m1 = max(
1e-6
, mu1*exp(bi));
            m2 = max(
1e-6
, mu2*exp(bi));
            array p1[%eval(&cap.+
1
)];
            array p2[%eval(&cap.+
1
)];
            s1 = 
0
; s2 = 
0
;
            do k = 
0
 to &cap.; p1[k+
1
]=pdf('poisson', k, m1); s1 + p1[k+
1
]; end;
            do k = 
0
 to &cap.; p2[k+
1
]=pdf('poisson', k, m2); s2 + p2[k+
1
]; end;
            do k = 
0
 to &cap.; p1[k+
1
]=p1[k+
1
]/s1; p2[k+
1
]=p2[k+
1
]/s2; end;
            y1 = rand('table', of p1[*]) - 
1
;
            y2 = rand('table', of p2[*]) - 
1
;
          %end;
          %else %if "&_DIST." = "CAT5" %then %do;
            /* Categorical uniform on \0..cap\ */
            y1 = rand('integer', %eval(&cap.+
1
)) - 
1
;
            y2 = rand('integer', %eval(&cap.+
1
)) - 
1
;
          %end;
          %else %do;
            /* Fallback: default to BIN4 */
            if id=1  and t=1  then put "WARNING: dist=&dist not recognized. Defaulting to BIN4.";
            m1 = max(
0
,  mu1*exp(bi));
            m2 = max(
0
,  mu2*exp(bi));
            p1 = min(max(m1/&cap., 
0
), 
1
);
            p2 = min(max(m2/&cap., 
0
), 
1
);
            y1 = rand('binomial', p1, &cap.);
            y2 = rand('binomial', p2, &cap.);
          %end;
        end;

        output;
      end;
    end;
    keep id class qtr y1 y2 obs;
  run;

  /* reshape to wide (quarters as columns) */
  proc sort data=sim_long; by id qtr; run;

  proc transpose data=sim_long(where=(obs=1)) out=y1_w prefix=Y1_;
    by id; id qtr; var y1;
  run;

  proc transpose data=sim_long(where=(obs=1)) out=y2_w prefix=Y2_;
    by id; id qtr; var y2;
  run;

  data sim_wide;
    merge y1_w y2_w; by id;
  run;

%mend sim_data;

/* Demo: build a toy dataset (comment out if you already have real data) */
/* %sim_data(dist=poisson, n=5000, T=12, seed=1); */

/*===============================================================================
 [B] INPUT PREP EXAMPLE
    - Expect a wide dataset with Y1_1..Y1_12 and an ID named id.
    - This block renames to BASE_FILE_SRS with SUM_Q1..SUM_Q12 and BENE_ID.
    - Comment out if your data is already prepared.
===============================================================================*/

data
 BASE_FILE_SRS;
  set WORK.Y1_W;  /* replace with your source if needed */
  rename
    Y1_1 = SUM_Q1  Y1_2 = SUM_Q2  Y1_3  = SUM_Q3  Y1_4  = SUM_Q4  Y1_5  = SUM_Q5
    Y1_6 = SUM_Q6  Y1_7 = SUM_Q7  Y1_8  = SUM_Q8  Y1_9  = SUM_Q9  Y1_10 = SUM_Q10
    Y1_11= SUM_Q11 Y1_12= SUM_Q12 id    = BENE_ID
  ;
run
;

/* Add time index variables (0..11) if you want them for other tooling */
data
 BASE_FILE_SRS(keep=BENE_ID SUM_Q1-SUM_Q12 quar1-quar12);
  set BASE_FILE_SRS;
  quar1=0 ; quar2=1 ; quar3=2 ; quar4=3 ; quar5=4 ; quar6=5 ; quar7=6 ; quar8=7 ;
  quar9=8 ; quar10=9 ; quar11=10 ; quar12=11 ;
run
;

/*===============================================================================
 [C] STARTING VALUES (USER EDITS HERE)
    - Example below matches nclass=6, order=2, p_order=0.
    - Define only the parameters you want; all others default to 0.
    - Class 1 is the reference for mixing weights (no alpha0_A).
===============================================================================*/

/* mixing weights */
%let alpha0_B = -0.50;
%let alpha0_C = -0.80;
/* %let alpha0_D = ; %let alpha0_E = ; %let alpha0_F = ; */

/* mean curve betas up to ORDER=2 for each class */
%let beta0_A =  0.10; %let beta1_A = 0.02; %let beta2_A = 0;
%let beta0_B = -0.20; %let beta1_B = 0.01; %let beta2_B = 0;
/* %let beta0_C = ; %let beta1_C = ; %let beta2_C = ; */
/* %let beta0_D = ; %let beta1_D = ; %let beta2_D = ; */
/* %let beta0_E = ; %let beta1_E = ; %let beta2_E = ; */
/* %let beta0_F = ; %let beta1_F = ; %let beta2_F = ; */

/* ZIP logits p(t) up to P_ORDER=0 (gamma0 only here) */
%let gamma0_A = -1.00;
%let gamma0_B = -0.50;
/* %let gamma0_C = ; %let gamma0_D = ; %let gamma0_E = ; %let gamma0_F = ; */

/*===============================================================================
 [D] MODELING MACROS
    - Helpers
    - Safe PARMS builder (uses %let values or defaults to 0)
    - Likelihood builders
    - Main fitting macro: %ct_zip_pois_nlmixed
===============================================================================*/

/* helper to emit one value per parameter into a single PARMS statement */
%macro _emit_parm(name, default);
  /* name is a token like alpha0_B */
  %if %symexist(&name) %then %do;                /* variable exists? */
    %if %length(%superq(&name)) %then %do;       /* nonblank value? */
      &name = %superq(&name)
    %end;
    %else %do;                                   /* blank -> default */
      &name = &default
    %end;
  %end;
  %else %do;                                     /* does not exist -> default */
    &name = &default
  %end;
%mend;

/* array and bookkeeping helpers */
%macro _ct_array_from_list(name, list, T); array &name.[&T] &list.; 
%mend;

%macro _ct_declare_mu_arrays(nclass, class_labels, T);
  %local k lab;
  %do k=1  %to &nclass;
    %let lab=%scan(&class_labels., &k, %str( ));
    array mu_&lab.[&T];
  %end;
%mend;

%macro _ct_declare_pi_arrays(nclass, class_labels, T);
  %local k lab;
  %do k=1  %to &nclass;
    %let lab=%scan(&class_labels., &k, %str( ));
    array pi_&lab.[&T];
  %end;
%mend;

%macro _ct_declare_zip_arrays(nclass, class_labels, T);
  %local k lab;
  %do k=1  %to &nclass;
    %let lab=%scan(&class_labels., &k, %str( ));
    array logitp_&lab.[&T];
    array p_&lab.[&T];
  %end;
%mend;

/* build one PARMS statement using either %let starts or 0 defaults */
%macro _ct_parms_all(nclass, class_labels, order, p_order);
  parms
  %local k lab;
  /* Mixing alphas: class 1 is reference -> no alpha0_A */
  %do k=2  %to &nclass;
    %let lab=%scan(&class_labels., &k, %str( ));
    %_emit_parm(%sysfunc(catx(_,alpha0,&lab)), 
0
)
  %end;
  /* Mean curve betas */
  %do k=1  %to &nclass;
    %let lab=%scan(&class_labels., &k, %str( ));
    %_emit_parm(%sysfunc(catx(_,beta0,&lab)), 
0
)
    %_emit_parm(%sysfunc(catx(_,beta1,&lab)), 
0
)
    %if &order>=
2
 %then %_emit_parm(%sysfunc(catx(_,beta2,&lab)), 
0
);
    %if &order>=
3
 %then %_emit_parm(%sysfunc(catx(_,beta3,&lab)), 
0
);
  %end;
  /* ZIP logits p_k(t) */
  %if &p_order>=
0
 %then %do;
    %do k=1  %to &nclass;
      %let lab=%scan(&class_labels., &k, %str( ));
      %_emit_parm(%sysfunc(catx(_,gamma0,&lab)), 
0
)
      %if &p_order>=
1
 %then %_emit_parm(%sysfunc(catx(_,gamma1,&lab)), 
0
);
      %if &p_order>=
2
 %then %_emit_parm(%sysfunc(catx(_,gamma2,&lab)), 
0
);
      %if &p_order>=
3
 %then %_emit_parm(%sysfunc(catx(_,gamma3,&lab)), 
0
);
    %end;
  %end;
  ;
%mend;

/* fill class-specific mu(t) over T using polynomial of degree order */
%macro _ct_fill_mu_poly(nclass, class_labels, order, T);
  %local k lab;
  do i = 
1
 to &T; t = i;
    %do k=1  %to &nclass;
      %let lab=%scan(&class_labels., &k, %str( ));
      eta_&lab = beta0_&lab. + beta1_&lab.*t
                 %if &order>=
2
 %then + beta2_&lab.*(t*t);
                 %if &order>=
3
 %then + beta3_&lab.*(t*t*t);
                 ;
      mu_&lab.[i] = exp(eta_&lab);
    %end;
  end;
%mend;

/* fill class-specific logit p_k(t) over T using polynomial of degree p_order */
%macro _ct_fill_zip_poly(nclass, class_labels, p_order, T);
  %local k lab;
  do i = 
1
 to &T; t = i;
    %do k=1  %to &nclass;
      %let lab=%scan(&class_labels., &k, %str( ));
      logitp_&lab.[i] = gamma0_&lab.
                        %if &p_order>=
1
 %then + gamma1_&lab.*t;
                        %if &p_order>=
2
 %then + gamma2_&lab.*(t*t);
                        %if &p_order>=
3
 %then + gamma3_&lab.*(t*t*t);
                        ;
      p_&lab.[i] = 
1
/(
1
+exp(-logitp_&lab.[i]));
    %end;
  end;
%mend;

/* accumulate ZIP log-likelihood per class across observed time points */
%macro _ct_accumulate_zip_ll(nclass, class_labels, T);
  %local k lab;
  do i = 
1
 to &T;
    if not missing(Y[i]) then do;
      %do k=1  %to &nclass;
        %let lab=%scan(&class_labels., &k, %str( ));
        lp = -log(
1
 + exp(-logitp_&lab.[i]));   /* log p   */
        lq = -log(
1
 + exp( logitp_&lab.[i]));   /* log(1-p) */
        if Y[i]=
0
 then do;                      /* log-sum-exp for 0 mass */
          a = lp; b = lq - mu_&lab.[i]; m = max(a,b);
          pi_&lab.[i] = m + log(exp(a-m) + exp(b-m));
        end;
        else do;                                /* Poisson mass for Y>0 */
          pi_&lab.[i] = lq + (Y[i]*log(mu_&lab.[i]) - mu_&lab.[i] - lgamma(Y[i]+
1
));
        end;
      %end;
    end;
  end;
%mend;

/* mix across classes via softmax weights from alphas */
%macro _ct_mixture_ll_zip(nclass, class_labels, T);
  %local k lab;
  den = 
1
;
  %do k=2  %to &nclass; %let lab=%scan(&class_labels., &k, %str( )); den = den + exp(alpha0_&lab.); %end;

  %do k=1  %to &nclass; %let lab=%scan(&class_labels., &k, %str( ));
    %if &k=1  %then %do; w_&lab.=
1
/den; %end;
    %else %do;          w_&lab.=exp(alpha0_&lab.)/den; %end;
  %end;

  %do k=1  %to &nclass; %let lab=%scan(&class_labels., &k, %str( ));
    prod_&lab.=
0
;
    do i=1  to &T; if not missing(pi_&lab.[i]) then prod_&lab. + pi_&lab.[i]; end;
  %end;

  m = prod_%scan(&class_labels.,
1
,%str( ));
  %do k=2  %to &nclass; %let lab=%scan(&class_labels., &k, %str( )); m = max(m, prod_&lab.); %end;

  sum_exp = 
0
;
  %do k=1  %to &nclass; %let lab=%scan(&class_labels., &k, %str( ));
    sum_exp + w_&lab.*exp(prod_&lab. - m);
  %end;

  ll = m + log(sum_exp);
%mend;

/* main modeling macro: runs NLMIXED and captures ParameterEstimates and FitStatistics */
%macro ct_zip_pois_nlmixed(
  data=BASE_FILE_SRS,
  id=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=6,
  class_labels=A B C D E F,
  order=2,
  p_order=0,
  T=12,
  tech=newrap,
  maxiter=500,
  pe_out=pe_zip,
  fit_out=fit_zip,
  bounds=
);
  ods listing;
  ods output ParameterEstimates=&pe_out FitStatistics=&fit_out;

  proc nlmixed data=&data qpoints=1  tech=&tech maxiter=&maxiter;
    /* Build single PARMS with user starts or 0 defaults */
    %_ct_parms_all(&nclass, &class_labels, &order, &p_order);

    /* Optional stability bounds (example: bounds %str(-6 < beta0_A < 6)) */
    %if %length(&bounds) %then %do; bounds &bounds; %end;

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

    one = 
1
;
    model one ~ general(ll);
    id &id;
  run;

  ods output close;
%mend;

/*===============================================================================
 [E] PLOTTING MACRO
    - Reads ParameterEstimates, reconstructs per-class trajectories:
        mu_zip(t) = (1 - p_k(t)) * exp(eta_k(t))
    - Computes mixture mean across classes using softmax weights.
    - Produces two plots (class curves; class curves + mixture mean).
===============================================================================*/

%macro ct_zip_plots(pe=pe_zip, T=12, out_traj=traj, out_mix=mix);
  ods listing;

  /* extract class-specific betas */
  proc sql;
    create table _betas as
    select scan(Parameter,
2
,'_') as class length=32,
           input(compress(substr(Parameter,
5
),,'kd'), best.) as deg,
           Estimate
    from &pe
    where upcase(substr(Parameter,
1
,
4
))='BETA';
  quit;
  proc sort data=_betas; by class deg; run;
  proc transpose data=_betas out=_betas_w prefix=b;
    by class; id deg; var Estimate;
  run;

  /* mixing weights (alpha logits) */
  proc sql;
    create table _alphas as
    select scan(Parameter,
2
,'_') as class length=32,
           Estimate as alpha
    from &pe
    where upcase(substr(Parameter,
1
,
6
))='ALPHA0';
  quit;
  proc sort data=_alphas; by class; run;

  /* ZIP logits p(t) coefficients (gamma) */
  proc sql;
    create table _gammas as
    select scan(Parameter,
2
,'_') as class length=32,
           input(compress(substr(Parameter,
7
),,'kd'), best.) as deg,
           Estimate
    from &pe
    where upcase(substr(Parameter,
1
,
5
))='GAMMA';
  quit;
  proc sort data=_gammas; by class deg; run;
  proc transpose data=_gammas out=_gammas_w prefix=g;
    by class; id deg; var Estimate;
  run;

  /* combine class parameters; alpha missing -> 0 for reference class */
  data _classes;
    merge _betas_w(in=b) _alphas(in=a) _gammas_w(in=g);
    by class;
    if missing(alpha) then alpha=0 ;
    exp_alpha = exp(alpha);
  run;

  /* denominator for class weights */
  proc sql noprint;
    select sum(exp_alpha) into :_den from _classes;
  quit;

  /* build trajectories over time */
  data &out_traj;
    set _classes;
    length class $
32
;
    do t = 
1
 to &T;
      eta    = coalesce(b0,
0
) + coalesce(b1,
0
)*t + coalesce(b2,
0
)*(t*t) + coalesce(b3,
0
)*(t*t*t);
      mu     = exp(eta);
      logitp = coalesce(g0,
0
) + coalesce(g1,
0
)*t + coalesce(g2,
0
)*(t*t) + coalesce(g3,
0
)*(t*t*t);
      p      = 
1
/(
1
+exp(-logitp));
      w      = exp_alpha / &_den;
      mu_zip = (
1
-p)*mu;
      output;
    end;
    keep class t mu p mu_zip w;
  run;

  /* mixture mean across classes */
  proc sql;
    create table &out_mix as
    select t, sum(w*mu_zip) as mu_mix_zip
    from &out_traj
    group by t;
  quit;

  /* Plot 1: class ZIP means */
  proc sgplot data=&out_traj;
    series x=t y=mu_zip / group=class lineattrs=(thickness=2);
    xaxis integer label="Quarter" min=1  max=&T;
    yaxis label="Expected count (ZIP mean)";
    title "Latent-Class ZIP Trajectories";
  run;

  /* Plot 2: class curves + mixture mean */
  proc sort data=&out_traj; by t; run;
  data traj_all; merge &out_traj &out_mix; by t; run;

  proc sgplot data=traj_all;
    series x=t y=mu_mix_zip / lineattrs=(pattern=shortdash thickness=3) name="mix" legendlabel="Mixture mean";
    series x=t y=mu_zip     / group=class lineattrs=(thickness=2);
    keylegend / position=topright;
    xaxis integer label="Quarter" min=1  max=&T;
    yaxis label="Expected count (ZIP mean)";
    title "Latent-Class ZIP Trajectories with Mixture Mean";
  run;
%mend;

/*===============================================================================
 [F] EXAMPLE RUN
    - Adjust nclass, class_labels, order, p_order, and T as needed.
    - Set your starts in Section [C].
    - Optionally add bounds via the bounds= parameter.
===============================================================================*/

%let T=12;

%ct_zip_pois_nlmixed(
  data=BASE_FILE_SRS,
  id=BENE_ID,
  yvars=SUM_Q1-SUM_Q12,
  nclass=6,
  class_labels=A B C D E F,
  order=2,
  p_order=0,
  T=&T
  /* bounds=%str(-6 < alpha0_B alpha0_C < 6) */
);

/* Produce the plots from ParameterEstimates (default pe_zip) */
%ct_zip_plots(pe=pe_zip, T=&T);

/*===============================================================================
 Optional: Save plots as PNG to WORK and embed in Results (uncomment to use)

%let outdir=%sysfunc(getoption(work));
ods _all_ close;
ods html5 path="&outdir" (url=none) gpath="&outdir";
ods graphics / reset imagename="ZIP_Trajectories" imagefmt=png;
%ct_zip_plots(pe=pe_zip, T=&T);
ods html5 close;

===============================================================================*/