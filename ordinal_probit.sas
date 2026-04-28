/*==============================================================================
TITLE:  Ordinal-Probit Latent Class (Single-Run: simulate -> fit -> plot)
WHAT:
  1) %sim_data (optional): make simulated data SIM_WIDE.
  2) %build_base_from_simwide: build BASE_FILE_SRS with quar1..quarT and BENE_ID.
  3) %ordprob_mix_fit_one: ONE NLMIXED run (K fixed). Class A alpha fixed to 0.
     - Ordered thresholds via cumulative exp() increments
     - Exactly ONE PARMS statement: defaults + user overrides (no duplicates)
  4) %ordprob_mix_plot_one: class proportions + predicted mean trajectories.
Author:     Haiqun Lin, Anum Zafar

==============================================================================*/

/*============================ OPTIONAL SIMULATOR ============================*/
%macro sim_data(
  dist=BIN4, class=4, n=2000, T=12, seed=1,
  miss_pattern=balanced, p_obs_min=0.6, sigma_re=0.4, cap=4
);
  %local _DIST; %let _DIST=%sysfunc(upcase(%sysfunc(strip(&dist.))));
  data sim_long;
    call streaminit(&seed.);
    do id = 1 to &n.;
      class = ceil(rand('uniform') * &class.);
      if "&miss_pattern."="balanced" then do; first_t=1; last_t=&T.; end;
      else do; frac=rand('uniform')*(1-&p_obs_min.) + &p_obs_min.; n_obs=ceil(&T.*frac); first_t=1; last_t=n_obs; end;
      bi = rand('normal', 0, &sigma_re.);
      do t = 1 to &T.;
        qtr=t; obs=(t>=first_t and t<=last_t); y1=.; y2=.;
        if obs then do;
          select (class);
            when (1) do; mu1=0.6+0.06*t;          mu2=0.8+0.03*t; end;
            when (2) do; mu1=0.5+0.30*t;          mu2=0.7+0.10*t; end;
            when (3) do; mu1=1.0+0.10*t+0.01*t*t; mu2=0.9+0.08*t; end;
            otherwise do; mu1=0.5+0.20*t;         mu2=0.5+0.15*t; end;
          end;

          %if "&_DIST."="BIN4" %then %do;
            m1=max(0,mu1*exp(bi)); m2=max(0,mu2*exp(bi));
            p1=min(max(m1/&cap.,0),1); p2=min(max(m2/&cap.,0),1);
            y1=rand('binomial',p1,&cap.); y2=rand('binomial',p2,&cap.);
          %end;
          %else %if "&_DIST."="TPOIS4" %then %do;
            m1=max(1e-6,mu1*exp(bi)); m2=max(1e-6,mu2*exp(bi));
            array p1[%eval(&cap.+1)]; array p2[%eval(&cap.+1)];
            s1=0; s2=0;
            do k=0 to &cap.; p1[k+1]=pdf('poisson',k,m1); s1 + p1[k+1];
                              p2[k+1]=pdf('poisson',k,m2); s2 + p2[k+1]; end;
            do k=0 to &cap.; p1[k+1]=p1[k+1]/s1; p2[k+1]=p2[k+1]/s2; end;
            y1=rand('table', of p1[*]) - 1; y2=rand('table', of p2[*]) - 1;
          %end;
          %else %if "&_DIST."="CAT5" %then %do;
            y1=rand('integer',%eval(&cap.+1))-1; y2=rand('integer',%eval(&cap.+1))-1;
          %end;
          %else %do;
            if id=1 and t=1 then put "WARNING: dist=&dist not recognized. Using BIN4.";
            m1=max(0,mu1*exp(bi)); m2=max(0,mu2*exp(bi));
            p1=min(max(m1/&cap.,0),1); p2=min(max(m2/&cap.,0),1);
            y1=rand('binomial',p1,&cap.); y2=rand('binomial',p2,&cap.);
          %end;
        end;
        output;
      end;
    end;
    keep id class qtr y1 y2 obs;
  run;

  proc sort data=sim_long; by id qtr; run;
  proc transpose data=sim_long(where=(obs=1)) out=y1_w prefix=Y1_; by id; id qtr; var y1; run;
  proc transpose data=sim_long(where=(obs=1)) out=y2_w prefix=Y2_; by id; id qtr; var y2; run;
  data sim_wide; merge y1_w y2_w; by id; run;
%mend sim_data;

/*============================== DATA PREP EXAMPLE ============================*/
%macro build_base_from_simwide(T=12);
  %if %sysfunc(exist(work.sim_wide)) %then %do;
    data BASE_FILE_SRS;
      set sim_wide;
      array quar[&T] quar1-quar&T;
      do i=1 to &T; quar[i]=i-1; end;
      drop i;
      rename id = BENE_ID;
    run;
  %end;
  %else %put NOTE: SIM_WIDE not found. Point DATA= in ordprob_mix_fit_one() to your own table.;
%mend build_base_from_simwide;

/*============================= SINGLE-RUN FIT ================================*/
/* Class index -> letter (A..O) */
%macro CL(c); %scan(A B C D E F G H I J K L M N O, &c., %str( )) %mend;

/* Mark user-specified overrides so defaults skip them */
%macro _index_override_pairs(pairs);
  %local i token eqpos name;
  %do i=1 %to %sysfunc(countw(&pairs,%str( )));
    %let token=%scan(&pairs,&i,%str( ));
    %let eqpos=%index(&token,=);
    %if &eqpos>1 %then %do;
      %let name=%substr(&token,1,%eval(&eqpos-1));
      %global OV_&name; %let OV_&name=1;
    %end;
  %end;
%mend;

/* Build defaults for (K,M,deg), skipping anything overridden */
%macro _make_default_parms(k, m, deg);
  %local c L l j s m1; %let s=; %let m1=%eval(&m-1);
  %do c=1 %to &k;
    %let L=%CL(&c); %let l=%lowcase(&L);
    /* mixing: skip A (fixed to 0) and skip if user overrode */
    %if &c>1 %then %if not %symexist(OV_alpha0_&L) %then %let s=&s alpha0_&L.=0;
    /* betas */
    %if not %symexist(OV_beta_&l.0) %then %let s=&s beta_&l.0=0;
    %if &deg>=1 %then %do; %if not %symexist(OV_beta_&l.1) %then %let s=&s beta_&l.1=0; %end;
    %if &deg>=2 %then %do j=2 %to &deg; %if not %symexist(OV_beta_&l.&j) %then %let s=&s beta_&l.&j.=0; %end;
    /* threshold increments i l1..i l(m-1) */
    %if &m1>=1 %then %do j=1 %to &m1; %if not %symexist(OV_i&l.&j) %then %let s=&s i&l.&j.=0; %end;
  %end;
  &s
%mend;

/* ONE NLMIXED run (guarded, single PARMS) */
%macro ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=BENE_ID,
  yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
  m=4, ycodes=0 1 2 3, deg=2,
  k=4, qpoints=40, maxiter=1000, tech=dbldog,
  start_values=,                     /* e.g., %str(alpha0_B=-0.5 beta_a0=0.1 ia1=0 ia2=0 ia3=0) */
  bounds=,                           /* e.g., %str(-6 < beta_a0 beta_b0 beta_c0 beta_d0 < 6)    */
  outestlib=work, prefix=ordprob_single
);
  %local m1 d j; %let m1=%eval(&m-1); %let d=&deg;

  /* Hard stop if data missing */
  %if not %sysfunc(exist(&data)) %then %do;
    %put ERROR: DATA=&data not found. Create it (e.g., run %sim_data and %build_base_from_simwide) or point DATA= to your table.;
    %return;
  %end;

  /* ycodes -> macro vars */
  %do j=1 %to &m; %global ycode&j; %let ycode&j=%scan(&ycodes,&j,%str( )); %end;

  /* Build override index so defaults skip user-specified params */
  %if %length(&start_values) %then %_index_override_pairs(%superq(start_values));

  ods listing;
  ods output
    ParameterEstimates  =&outestlib..&prefix._EST_K&k
    FitStatistics       =&outestlib..&prefix._FIT_K&k
    AdditionalEstimates =&outestlib..&prefix._ESTS_K&k;

  proc nlmixed data=&data qpoints=&qpoints noad maxiter=&maxiter tech=&tech;
    array yv[&ttotal] &yvars;
    array tv[&ttotal] &tvars;

    /* ONE PARMS: defaults (sans duplicates) + user overrides */
    parms %_make_default_parms(&k, &m, &deg)
          %superq(start_values);

    /* optional bounds */
    %if %length(&bounds) %then %do; bounds &bounds; %end;

    /* fix class A mixing */
    alpha0_A = 0;

    /* softmax mixing */
    denom=0;
    %do c=1 %to &k; %let U=%CL(&c); pin&c = exp(alpha0_&U); denom + pin&c; %end;
    %do c=1 %to &k; pie&c = pin&c/denom; %end;

    /* class log-likelihoods */
    %do c=1 %to &k; llik&c=0; %end;

    do t=1 to &ttotal;
      tval = tv[t];
      cat = .; %do j=1 %to &m; if yv[t] = &&ycode&j then cat=&j; %end;

      if cat>0 then do;
        %do c=1 %to &k;
          %let U=%CL(&c); %let L=%lowcase(&U);

          /* eta(t): polynomial in t */
          eta_&U = beta_&L.0;
          %if &d>=1 %then %do; eta_&U = eta_&U + beta_&L.1*(tval); %end;
          %if &d>=2 %then %do j=2 %to &d; eta_&U = eta_&U + beta_&L.&j.*((tval)**&j); %end;

          /* ordered thresholds: th1 = iL1; thj = th(j-1) + exp(iLj) */
          th1_&U = i&L.1;
          %if &m1>=2 %then %do j=2 %to &m1;
            %if &j=2 %then %do; th&j._&U = th1_&U + exp(i&L.&j); %end;
            %else %do;          th&j._&U = th%eval(&j-1)_&U + exp(i&L.&j); %end;
          %end;

          /* category probabilities */
          %do j=1 %to &m;
            %if &j=1 %then %do;
              z&j._&U = probnorm(th1_&U - eta_&U);
            %end; %else %if &j<&m %then %do;
              z&j._&U = probnorm(th&j._&U - eta_&U) - probnorm(th%eval(&j-1)_&U - eta_&U);
            %end; %else %do;
              z&j._&U = 1 - probnorm(th%eval(&m1)_&U - eta_&U);
            %end;
          %end;

          array z&U[&m] %do j=1 %to &m; z&j._&U %end;;
          p_obs&c = max(1e-12, z&U[cat]);
          llik&c  = llik&c + log(p_obs&c);
        %end;
      end;
    end;

    mix=0; %do c=1 %to &k; mix + pie&c*exp(llik&c); %end;
    mix=max(mix,1e-300);
    dummy=0;
    model dummy ~ general(log(mix));

    /* thresholds on natural scale for convenience */
    %do c=1 %to &k; %let U=%CL(&c); %let L=%lowcase(&U);
      estimate "threshold1_&L"  i&L.1;
      %if &m1>=2 %then %do _jj=2 %to &m1;
        %local _expr; %let _expr = i&L.1;
        %do j=2 %to &_jj; %let _expr = &_expr + exp(i&L.&j); %end;
        estimate "threshold&_jj._&L" &_expr;
      %end;
    %end;
  run;

  ods output close;
%mend;

/*============================= PLOTTING (SINGLE RUN) ========================*/
%macro ordprob_mix_plot_one(
  outestlib=work, prefix=ordprob_single,
  k=4, ttotal=12, m=4, ycodes=0 1 2 3
);
  %local m1 j; %let m1=%eval(&m-1);
  %do j=1 %to &m; %global ycode&j; %let ycode&j=%scan(&ycodes,&j,%str( )); %end;

  /* mixing proportions from alpha0_*; add A if missing (fixed to 0) */
  data _mix;
    set &outestlib..&prefix._EST_K&k(keep=Parameter Estimate);
    length class $1 alpha 8;
    if upcase(substr(Parameter,1,7))='ALPHA0_' then do;
      class = substr(Parameter,8,1); alpha = Estimate; output;
    end;
  run;
  proc sql noprint; select count(*) into :_hasA from _mix where upcase(class)='A'; quit;
  %if %sysevalf(&_hasA=0) %then %do;
    data _addA; length class $1 alpha 8; class='A'; alpha=0; run;
    proc append base=_mix data=_addA force; run;
    proc datasets lib=work nolist; delete _addA; quit;
  %end;

  data _mix; set _mix; class=upcase(class); run;
  proc sort data=_mix; by class; run;

  data mix_props; set _mix end=last; retain sumexp 0; expa=exp(alpha); sumexp + expa; if last then call symputx('_sumexp',sumexp); run;
  data mix_props; set mix_props; pie = expa / &_sumexp; run;

  proc sgplot data=mix_props;
    vbar class / response=pie datalabel;
    yaxis grid label='Mixing proportion pi_c';
    xaxis label='Class';
    title "Estimated Class Proportions (K=&k)";
  run;

  /* Predicted trajectories */
  proc sql;
    create table _betas as
    select upcase(substr(Parameter,6,1)) as CL length=1,
           input(substr(Parameter,7), best.) as IDX,
           Estimate
    from &outestlib..&prefix._EST_K&k
    where upcase(substr(Parameter,1,5))='BETA_';
  quit;
  proc sort data=_betas; by CL IDX; run;
  proc transpose data=_betas out=betas_w prefix=beta_; by CL; id IDX; var Estimate; run;

  proc sql;
    create table _ths as
    select upcase(scan(Label,2,'_')) as CL length=1,
           input(compress(substr(lowcase(Label),10),,'kd'), best.) as IDX,
           Estimate
    from &outestlib..&prefix._ESTS_K&k
    where upcase(substr(Label,1,9))='THRESHOLD';
  quit;
  proc sort data=_ths; by CL IDX; run;
  proc transpose data=_ths out=ths_w prefix=th_; by CL; id IDX; var Estimate; run;

  proc sort data=betas_w; by CL; run; proc sort data=ths_w; by CL; run;
  data coeffs; merge betas_w ths_w; by CL; run;

  data pred_long;
    set coeffs;
    length class $1; class=CL;
    array b  beta_0-beta_99;  /* only existing indices are used */
    array th th_1-th_99;
    array p  p1-p&m;

    do tval=0 to %eval(&ttotal-1);
      eta=0; do j=1 to dim(b); if not missing(b[j]) then eta + b[j]*(tval**(j-1)); end;

      p[1] = probnorm(th[1] - eta);
      %if &m1>=2 %then %do jj=2 %to &m1;
        p[&jj] = probnorm(th[&jj] - eta) - probnorm(th[%eval(&jj-1)] - eta);
      %end;
      p[&m] = 1 - probnorm(th[&m1] - eta);

      EY=0; %do j=1 %to &m; EY + %scan(&ycodes,&j,%str( ))*p[&j]; %end;
      output;
    end;
    keep class tval EY p1-p&m;
  run;

  proc sgplot data=pred_long;
    series x=tval y=EY / group=class markers;
    xaxis integer label='Quarter (t)';
    yaxis label='E[Y_t | class]';
    title "Predicted Mean Trajectories by Class (K=&k)";
  run;
%mend;

/*============================== EXAMPLE: SINGLE RUN =========================*/
/* EITHER: build simulated data… */
%sim_data();                  /* comment this out if using your own data */
%build_base_from_simwide(T=12);

/* …OR: point DATA= to your wide table and ensure:
   - yvars = your wide outcomes (e.g., Y1_1-Y1_12)
   - tvars = quar1-quar12 with values 0..11 (create them if needed)
   - id    = your subject ID (rename to BENE_ID or pass id=)      */

%let K   = 4;
%let M   = 4;
%let DEG = 2;

/* Optional: user starts (list only what you want; do NOT include alpha0_A) */
%let STARTS = %str(
  alpha0_B=-0.5 alpha0_C=-0.2 alpha0_D=0
  beta_a0=0  beta_a1=0  beta_a2=0
  ia1=0 ia2=0 ia3=0
);

/* Fit ONE model */
%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=BENE_ID,
  yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
  m=&M, ycodes=0 1 2 3, deg=&DEG, k=&K,
  qpoints=40, maxiter=1000, tech=dbldog,
  start_values=&STARTS,              /* leave blank to start all at zero */
  outestlib=work, prefix=ordprob_single
);

/* Plot */
%ordprob_mix_plot_one(outestlib=work, prefix=ordprob_single, k=&K, ttotal=12, m=&M, ycodes=0 1 2 3);
