/*****************************************************************************************************************************************
* Community Health and Aging Outcome (CHAO) Lab - Rutgers, The State University of New Jersey                                           *
* Title:   Group-Based Trajectory Modeling using PROC NLMIXED (Two Continuous Outcomes)                                                  *
* Purpose: Implements group-based trajectory modeling (latent-class) for TWO continuous outcomes using a censored-normal (Tobit) model. *
*          Includes: (1) single-outcome models (2) joint model (3) plotting + posterior summaries                                         *
* Data:    Wide format: 1 row/person with repeated measures for each outcome across T time points                                        *
* Outputs: work.nlm_fix_T1&class, work.nlm_fix_T2&class, work.nlm_fix_T1_T2&class, plots, avg_membership                                  *
* Authors: Weiyi Xia, Haiqun Lin, Anum Zafar                                                                                            *
* Last edit: 2026-02-09                                                                                                                 *
*****************************************************************************************************************************************/


/*=============================================================
  Step 0: GLOBAL SETTINGS (edit these only)
=============================================================*/
%let USE_SIM     = 1;            /* 1 = run simulator, 0 = use your real BASE_FILE_SRS */
%let T           = 12;           /* number of time points */
%let class       = 5;            /* number of latent classes */
%let order_model = 3;            /* 1=linear, 2=quadratic, 3=cubic */
%let equal       = T;            /* T=equal sigma across classes (per outcome), F=class-specific */

/* censoring caps, numeric list length T */
%let max_values = 10 10 10 10 10 10 10 10 10 10 10 10;

/* outcome variable lists in YOUR wide dataset */
%let y1vars = QHH1-QHH12;        /* Outcome 1 repeated measures */
%let y2vars = QINP1-QINP12;      /* Outcome 2 repeated measures */

/* ID + time index variables */
%let idvar = BENE_ID;
%let tvars = quar1-quar12;


/*=============================================================
  Step 1: OPTIONAL SIMULATOR (Wide + Long)
=============================================================*/
%macro sim_data_cont(
  class=5,
  n=500,
  T=12,
  seed=2026,
  miss_pattern=balanced,   /* balanced | unbalanced */
  p_obs_min=0.6,
  max_const=10,
  sigma1=1.0,
  sigma2=1.2
);
  data sim_long;
    call streaminit(&seed);

    do id = 1 to &n;

      class = ceil(rand('uniform') * &class);

      if "&miss_pattern" = "balanced" then do;
        first_t = 1; last_t = &T;
      end;
      else do;
        frac    = rand('uniform')*(1-&p_obs_min) + &p_obs_min;
        n_obs   = ceil(&T*frac);
        first_t = 1; last_t = n_obs;
      end;

      select (class);
        when (1) do; b10=-2.0; b11= 0.25; b12= 0.00;  b20=-1.0; b21= 0.10; b22=0.01; end;
        when (2) do; b10=-1.0; b11= 0.15; b12= 0.01;  b20=-1.8; b21= 0.22; b22=0.00; end;
        when (3) do; b10=-0.5; b11= 0.05; b12= 0.02;  b20=-0.7; b21= 0.08; b22=0.02; end;
        when (4) do; b10=-2.5; b11= 0.35; b12=-0.01;  b20=-1.2; b21= 0.18; b22=0.00; end;
        otherwise do; b10=-1.5; b11= 0.12; b12= 0.00;  b20=-1.3; b21= 0.12; b22=0.01; end;
      end;

      do qtr = 1 to &T;
        obs   = (qtr >= first_t and qtr <= last_t);
        cap_t = &max_const;

        y1 = .; y2 = .;

        if obs then do;
          mu1 = b10 + b11*qtr + b12*(qtr*qtr);
          mu2 = b20 + b21*qtr + b22*(qtr*qtr);

          z1 = mu1 + rand('normal', 0, &sigma1);
          z2 = mu2 + rand('normal', 0, &sigma2);

          y1 = max(0, min(cap_t, z1));
          y2 = max(0, min(cap_t, z2));
        end;

        output;
      end;

    end;
    keep id class qtr y1 y2 obs;
  run;

  proc sort data=sim_long; by id qtr; run;

  proc transpose data=sim_long(where=(obs=1)) out=_y1_w prefix=Y1_;
    by id; id qtr; var y1;
  run;

  proc transpose data=sim_long(where=(obs=1)) out=_y2_w prefix=Y2_;
    by id; id qtr; var y2;
  run;

  data sim_wide;
    merge _y1_w _y2_w;
    by id;
  run;
%mend;


/*=============================================================
  Step 2: BUILD / LOAD BASE_FILE_SRS
  + adds CAP1-CAP&T as real variables (important fix)
=============================================================*/
%macro add_caps(ds=BASE_FILE_SRS);
  data &ds;
    set &ds;
    array CAP[&T] CAP1-CAP&T;
    /* load from macro list */
    %do _i=1 %to &T;
      CAP[&_i] = %scan(&max_values, &_i);
    %end;
    drop _i;
  run;
%mend;

%macro build_base_file_srs;
  %if &USE_SIM = 1 %then %do;

    %sim_data_cont(class=&class, n=500, T=&T, seed=1, miss_pattern=balanced, max_const=10);

    data BASE_FILE_SRS;
      set sim_wide;
      rename
        id   = &idvar
        Y1_1 = QHH1   Y1_2 = QHH2   Y1_3 = QHH3   Y1_4 = QHH4   Y1_5 = QHH5   Y1_6 = QHH6
        Y1_7 = QHH7   Y1_8 = QHH8   Y1_9 = QHH9   Y1_10= QHH10  Y1_11= QHH11  Y1_12= QHH12
        Y2_1 = QINP1  Y2_2 = QINP2  Y2_3 = QINP3  Y2_4 = QINP4  Y2_5 = QINP5  Y2_6 = QINP6
        Y2_7 = QINP7  Y2_8 = QINP8  Y2_9 = QINP9  Y2_10= QINP10 Y2_11= QINP11 Y2_12= QINP12
      ;
    run;

    data BASE_FILE_SRS;
      set BASE_FILE_SRS;
      array quar[&T] &tvars;
      do _i=1 to &T;
        quar[_i] = _i;
      end;
      drop _i;
    run;

    %add_caps(ds=BASE_FILE_SRS);

  %end;
  %else %do;
    %put NOTE: USE_SIM=0, expecting BASE_FILE_SRS already exists with ID=&idvar, time=&tvars, y1=&y1vars, y2=&y2vars.;
    %add_caps(ds=BASE_FILE_SRS);
  %end;
%mend;

%build_base_file_srs;


/*=============================================================
  Step 3: MACRO LIBRARY
=============================================================*/
%let class_names = A B C D E F G H I J K L M N O P Q R S T;

%macro starting_value_alpha(class);
  %local s class_;
  %do s=2 %to &class;
    %let class_ = %scan(&class_names, &s);
    alpha0_&class_.=0
  %end;
%mend;

%macro starting_value_beta_sigma(class,outcome,order,equal_sigma);
  %local s i class_;

  %if %upcase(&equal_sigma)=T %then %do;
    sigma&outcome._ = 30
  %end;
  %else %do;
    sigma&outcome._A = 30
    %do s=2 %to &class;
      %let class_=%scan(&class_names, &s);
      sigma&outcome._&class_. = 30
    %end;
  %end;

  %do s=1 %to &class;
    %let class_=%scan(&class_names, &s);
    beta&outcome._&class_.0 = 0
    %do i=1 %to &order;
      beta&outcome._&class_.&i = 0
    %end;
  %end;
%mend;

%macro bounds_alpha(bounds_alpha,class);
  %local s class_;
  %do s=2 %to &class;
    %let class_=%scan(&class_names, &s);
    -&bounds_alpha.<alpha0_&class_.<&bounds_alpha.
    %if &s < &class %then ,;
  %end;
%mend;

%macro bounds_sigma(bounds_sigma,class,outcome,equal_sigma);
  %local s class_;
  %if %upcase(&equal_sigma)=T %then %do;
    sigma&outcome._ > &bounds_sigma.
  %end;
  %else %do;
    sigma&outcome._A > &bounds_sigma.
    %do s=2 %to &class;
      %let class_=%scan(&class_names, &s);
      , sigma&outcome._&class_. > &bounds_sigma.
    %end;
  %end;
%mend;

%macro initiation_universal(T);
  array X[&T] &tvars;
  array Y1[&T] &y1vars;
  array Y2[&T] &y2vars;
%mend;

%macro initiation(T,class,outcome);
  %local s class_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    array PI&class_.&outcome.[&T] PI&class_.&outcome._1-PI&class_.&outcome._&T;
    array mu&class_.&outcome.[&T] mu&class_.&outcome._1-mu&class_.&outcome._&T;
    PROD&class_.&outcome.=0;
  %end;
%mend;

%macro model(order,class,outcome);
  %local s i class_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    mu&class_.&outcome.[I] = beta&outcome._&class_.0
    %do i=1 %to &order;
      + beta&outcome._&class_.&i * (X[I]**&i)
    %end;
    ;
  %end;
%mend;

%macro residual(class,outcome);
  %local s class_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    e&outcome._&class_. = Y&outcome.[I] - mu&class_.&outcome.[I];
  %end;
%mend;

%macro float_control(float,class,outcome,equal_sigma);
  %local s class_ class2_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    %if %upcase(&equal_sigma)=T %then %let class2_=;
    %else %let class2_=&class_;
    e&outcome._&class_. = min(max(e&outcome._&class_.,
                          -&float.*sigma&outcome._&class2_.),
                           &float.*sigma&outcome._&class2_.);
  %end;
%mend;

/* censored-normal contribution uses CAP[I] (CAP array points to CAP1-CAP&T) */
%macro Prob_cnorm(class,outcome,equal_sigma);
  %local s class_ class2_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    %if %upcase(&equal_sigma)=T %then %let class2_=;
    %else %let class2_=&class_;

    PI&class_.&outcome.[I] = logpdf('NORMAL', e&outcome._&class_. / sigma&outcome._&class2_.)
                             - log(sigma&outcome._&class2_.);

    if Y&outcome.[I] = 0 then
      PI&class_.&outcome.[I] = logcdf('NORMAL', e&outcome._&class_. / sigma&outcome._&class2_.);
    else if Y&outcome.[I] = CAP[I] then
      PI&class_.&outcome.[I] = logcdf('NORMAL', (-e&outcome._&class_.) / sigma&outcome._&class2_.);

    PROD&class_.&outcome. = PROD&class_.&outcome. + PI&class_.&outcome.[I];
  %end;
%mend;

%macro Class_Membership(class);
  %local s class_;
  alpha0_A = 0;

  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    pinumer_&class_. = exp(alpha0_&class_.);
  %end;

  pideno = 0
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    + pinumer_&class_.
  %end;
  ;

  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    pie_&class_. = pinumer_&class_. / pideno;
  %end;
%mend;

%macro LogLike(class,outcome);
  %local s class_;
  l_latclass =
    %do s=1 %to &class;
      %let class_=%scan(&class_names,&s);
      %if &s>1 %then + ;
      pie_&class_. * exp(PROD&class_.&outcome.)
    %end;
  ;
%mend;

%macro LogLike_multi(class);
  %local s class_;
  ( %do s=1 %to &class;
      %let class_=%scan(&class_names,&s);
      %if &s>1 %then + ;
      pie_&class_. * exp(PROD&class_.1) * exp(PROD&class_.2)
    %end;
  )
%mend;

%macro Posterior(class);
  %local s class_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    post_&class_. = pie_&class_. * exp(PROD&class_.1) * exp(PROD&class_.2) / %LogLike_multi(class=&class);
  %end;

  keep
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    post_&class_.
  %end;
  ;
%mend;

%macro initiation_pred(T,class,outcome);
  %local s class_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    array Pred&class_.&outcome.[&T] Pred&class_.&outcome._1-Pred&class_.&outcome._&T;
    array mu&class_.&outcome.[&T]   mu&class_.&outcome._1-mu&class_.&outcome._&T;
  %end;
%mend;

%macro Pred_cnorm(class,outcome,equal_sigma);
  %local s class_ class2_;
  %do s=1 %to &class;
    %let class_=%scan(&class_names,&s);
    %if %upcase(&equal_sigma)=T %then %let class2_=;
    %else %let class2_=&class_;

    temp = (exp(logcdf('normal',(CAP[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)) -
            exp(logcdf('normal',(0      -mu&class_.&outcome.[I])/sigma&outcome._&class2_.)));
    if temp = 0 then temp = 0.000001;

    Pred&class_.&outcome.[I] =
      0
      + temp * (
          mu&class_.&outcome.[I]
          + sigma&outcome._&class2_. *
            (exp(logpdf('normal', -mu&class_.&outcome.[I]/sigma&outcome._&class2_.)) -
             exp(logpdf('normal', (CAP[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.))) / temp
        )
      + CAP[I] * exp(logcdf('normal', (-CAP[I]+mu&class_.&outcome.[I])/sigma&outcome._&class2_.))
    ;
  %end;
%mend;


/*=============================================================
  PLOTTING MACRO
=============================================================*/
%macro plot_prep(T,LC,result,order,equal_sigma);

  data parameter; set &result; keep parameter estimate; run;
  proc transpose data=parameter out=parameter; id parameter; run;

  proc sql;
    create table data_pred as
    select * from BASE_FILE_SRS, parameter;
  quit;

  data pred_membership_y;
    set data_pred;

    %initiation_universal(T=&T);
    array CAP[&T] CAP1-CAP&T;

    %initiation(T=&T, class=&LC, outcome=1);
    %initiation(T=&T, class=&LC, outcome=2);

    do I=1 to &T;
      %model(order=&order, class=&LC, outcome=1);
      %residual(class=&LC, outcome=1);
      %float_control(float=8, class=&LC, outcome=1, equal_sigma=&equal_sigma);
      %Prob_cnorm(class=&LC, outcome=1, equal_sigma=&equal_sigma);

      %model(order=&order, class=&LC, outcome=2);
      %residual(class=&LC, outcome=2);
      %float_control(float=8, class=&LC, outcome=2, equal_sigma=&equal_sigma);
      %Prob_cnorm(class=&LC, outcome=2, equal_sigma=&equal_sigma);
    end;

    %Class_Membership(class=&LC);
    %Posterior(class=&LC);
  run;

  data y_1; set BASE_FILE_SRS(keep=&y1vars); run;
  proc iml;
    start colsum(m); return(m[+,]); finish;
    use y_1; read all var _num_ into y;
    use pred_membership_y; read all var _num_ into p;
    sum = colsum(p);
    invsum = 1/sum;
    avg_y  = y`*p;
    avg_y2 = avg_y#invsum;
    create avg_y1 from avg_y2; append from avg_y2; close avg_y1;
  quit;

  data y_2; set BASE_FILE_SRS(keep=&y2vars); run;
  proc iml;
    start colsum(m); return(m[+,]); finish;
    use y_2; read all var _num_ into y;
    use pred_membership_y; read all var _num_ into p;
    sum = colsum(p);
    invsum = 1/sum;
    avg_y  = y`*p;
    avg_y2 = avg_y#invsum;
    create avg_y2 from avg_y2; append from avg_y2; close avg_y2;
  quit;

  /* build single-row dataset for prediction curves */
  data wide;
    set parameter;
    %do j=1 %to &T;
      quar&j = &j;
      CAP&j  = %scan(&max_values,&j);
    %end;
  run;

  data wide;
    set wide;
    array X[&T] &tvars;
    array CAP[&T] CAP1-CAP&T;

    %initiation_pred(T=&T, class=&LC, outcome=1);
    %initiation_pred(T=&T, class=&LC, outcome=2);

    do I=1 to &T;
      %model(order=&order, class=&LC, outcome=1);
      %model(order=&order, class=&LC, outcome=2);
      %Pred_cnorm(class=&LC, outcome=1, equal_sigma=&equal_sigma);
      %Pred_cnorm(class=&LC, outcome=2, equal_sigma=&equal_sigma);
    end;
  run;

  proc transpose data=wide out=pred_temp name=_NAME_;
    var _numeric_;
  run;

  data pred_temp;
    set pred_temp;
    rename COL1 = estimate;
  run;

  data pred;
    set pred_temp;
    length outcome $3 type $4 class $1;
    where _NAME_ contains 'Pred';
    class  = upcase(substr(_NAME_,5,1));
    y      = substr(_NAME_,6,1);
    if y='1' then outcome='HH';
    if y='2' then outcome='INP';
    type   = 'pred';
    quar   = input(substr(_NAME_,8), best.);
    keep outcome type class quar estimate;
  run;

  title 'Predicted Trajectory by Latent Class';
  proc sgpanel data=pred;
    panelby class / columns=4;
    series x=quar y=estimate / group=outcome;
  run;
  title;

  data avg_y1_plot;
    set avg_y1;
    length outcome $3 type $3 class $1;
    array col col1-col&LC;
    quar=_n_;
    do i=1 to &LC;
      estimate=col[i];
      class=upcase(byte(i+64));
      outcome="HH"; type="Avg";
      output;
    end;
    keep estimate class quar outcome type;
  run;

  data avg_y2_plot;
    set avg_y2;
    length outcome $3 type $3 class $1;
    array col col1-col&LC;
    quar=_n_;
    do i=1 to &LC;
      estimate=col[i];
      class=upcase(byte(i+64));
      outcome="INP"; type="Avg";
      output;
    end;
    keep estimate class quar outcome type;
  run;

  data avg; set avg_y1_plot avg_y2_plot; run;

  proc sort data=avg;  by outcome quar class; run;
  proc sort data=pred; by outcome quar class; run;

  data pred2;
    set pred;
    pred = estimate;
    keep outcome quar class pred;
  run;

  data data_plot;
    merge avg(rename=(estimate=avg)) pred2;
    by outcome quar class;
  run;

  title 'Predicted vs Averaged Observed by Latent Class';
  proc sgpanel data=data_plot;
    panelby outcome;
    series  x=quar y=pred / group=class name="pred";
    scatter x=quar y=avg  / group=class name="obs";
    keylegend "pred" / title="Predicted";
    keylegend "obs"  / title="Averaged Observed";
  run;
  title;

  title 'Averaged Posterior Class Membership';
  proc means data=pred_membership_y mean missing;
    var
    %do s=1 %to &LC;
      %let class_=%scan(&class_names,&s);
      post_&class_.
    %end;
    ;
    output out=avg_membership mean=Avg std=SD;
  run;
  title;

%mend;


/*=============================================================
  NLMIXED MACROS
=============================================================*/
%macro nlmixed_1(T,LC,Y,starting,output,order,equal_sigma);

  proc nlmixed data=BASE_FILE_SRS itdetails qpoints=40 noad maxiter=1000 tech=dbldog;

    bounds %bounds_alpha(bounds_alpha=3, class=&LC),
           %bounds_sigma(bounds_sigma=0, class=&LC, outcome=&Y, equal_sigma=&equal_sigma);

    parms &starting;

    %initiation_universal(T=&T);
    array CAP[&T] CAP1-CAP&T;

    %initiation(T=&T, class=&LC, outcome=&Y);

    do I=1 to &T;
      %model(order=&order, class=&LC, outcome=&Y);
      %residual(class=&LC, outcome=&Y);
      %float_control(float=8, class=&LC, outcome=&Y, equal_sigma=&equal_sigma);
      %Prob_cnorm(class=&LC, outcome=&Y, equal_sigma=&equal_sigma);
    end;

    %Class_Membership(class=&LC);
    %LogLike(class=&LC, outcome=&Y);

    ll_latclass = log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output;
  run;

%mend;

%macro nlmixed_MultiTraj(T,LC,starting,output,order,equal_sigma);

  proc nlmixed data=BASE_FILE_SRS itdetails qpoints=40 noad maxiter=1000 tech=dbldog;

    bounds %bounds_alpha(bounds_alpha=3, class=&LC),
           %bounds_sigma(bounds_sigma=0, class=&LC, outcome=1, equal_sigma=&equal_sigma),
           %bounds_sigma(bounds_sigma=0, class=&LC, outcome=2, equal_sigma=&equal_sigma);

    parms &starting;

    %initiation_universal(T=&T);
    array CAP[&T] CAP1-CAP&T;

    %initiation(T=&T, class=&LC, outcome=1);
    %initiation(T=&T, class=&LC, outcome=2);

    do I=1 to &T;
      %model(order=&order, class=&LC, outcome=1);
      %residual(class=&LC, outcome=1);
      %float_control(float=8, class=&LC, outcome=1, equal_sigma=&equal_sigma);
      %Prob_cnorm(class=&LC, outcome=1, equal_sigma=&equal_sigma);

      %model(order=&order, class=&LC, outcome=2);
      %residual(class=&LC, outcome=2);
      %float_control(float=8, class=&LC, outcome=2, equal_sigma=&equal_sigma);
      %Prob_cnorm(class=&LC, outcome=2, equal_sigma=&equal_sigma);
    end;

    %Class_Membership(class=&LC);

    l_latclass  = %LogLike_multi(class=&LC);
    ll_latclass = log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output;
  run;

%mend;


/*=============================================================
  Step 4: RUN PIPELINE
=============================================================*/
%nlmixed_1(
  T=&T, LC=&class, Y=1,
  starting=%starting_value_alpha(class=&class)
           %starting_value_beta_sigma(class=&class,outcome=1,order=&order_model,equal_sigma=&equal),
  output=nlm_fix_T1&class,
  order=&order_model,
  equal_sigma=&equal
);

%nlmixed_1(
  T=&T, LC=&class, Y=2,
  starting=%starting_value_alpha(class=&class)
           %starting_value_beta_sigma(class=&class,outcome=2,order=&order_model,equal_sigma=&equal),
  output=nlm_fix_T2&class,
  order=&order_model,
  equal_sigma=&equal
);

data work.nlm_2y_starting;
  set nlm_fix_T1&class nlm_fix_T2&class;
  if parameter =: 'alpha' then delete;
run;

%nlmixed_MultiTraj(
  T=&T, LC=&class,
  starting=%starting_value_alpha(class=&class) / data=work.nlm_2y_starting,
  output=nlm_fix_T1_T2&class,
  order=&order_model,
  equal_sigma=&equal
);

%plot_prep(
  T=&T, LC=&class,
  result=nlm_fix_T1_T2&class,
  order=&order_model,
  equal_sigma=&equal
);
