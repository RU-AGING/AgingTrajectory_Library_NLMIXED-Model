/*****************************************************************************************************************************************
* Community Health and Aging Outcome (CHAO) Lab - Rutgers, The State University of New Jersey                                           *
* Title:   Group-Based Trajectory Modeling using PROC NLMIXED (Two Continuous Outcomes)                                                  *
* Purpose: Implements group-based trajectory modeling (latent-class) for TWO continuous outcomes using a censored-normal (Tobit) model. *
*          Includes: (1) single-outcome models (2) joint model (3) plotting + posterior summaries                                         *
* Data:    Wide format: 1 row/person with repeated measures for each outcome across T time points                                        *
* Outputs: work.nlm_fix_T1&class, work.nlm_fix_T2&class, work.nlm_fix_T1_T2&class, plots, avg_membership                                  *                                                                                           *
LAST UPDATED DATE: 23 JUL 2026
DATA SOURCES: NONE. This file defines macros only. INPUT IS the optional simulator OR your own wide TABLE
PURPOSE: GROUP-based trajectory modelling for ONE OR TWO continuous outcomes under a censored-normal
Tobit likelihood, estimated BY PROC NLMIXED. Provides
(a)an optional simulator AND a wide-FORMAT DATA-prep helper
(b)single-outcome fits, one per outcome
(c)a joint two-outcome fit over a common latent-class structure
(d)observed-versus-predicted plots AND averaged posterior class-membership summaries
DATA: one row per subject, repeated measures for each outcome across T time points, a numeric
time INDEX quar1 TO quarT, AND censoring caps CAP1 TO CAPT which ADD_CAPS attaches AS real variables.
AUTHOR: Weiyi Xia, Haiqun Lin, Anum Zafar
##########################################################################################################################
*Execution Environment: SAS 9.4 OR later, PROC NLMIXED FROM SAS/STAT. No compiled components         *
*This file defines macros only. SET TRAJ2_RUN_DEMO=1 BEFORE the include TO also RUN the built-IN demo*
*Comments inside the code-generator macros must stay IN slash-star form. A star-semicolon comment    *
*inside a text-returning MACRO IS emitted INTO the statement it builds AND IS a syntax error         *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 0:GLOBAL settings AND the demo switch. Edit only this block
#STEP 1:Optional simulator. SIM_DATA_CONT builds SIM_WIDE
#STEP 2:DATA prep. BUILD_BASE_FILE_SRS AND ADD_CAPS build BASE_FILE_SRS including CAP1 TO CAPT
#STEP 3:Code-generator library. Starting values, BOUNDS, MODEL, likelihood, posterior, prediction
#STEP 4:Plotting. PLOT_PREP builds observed versus predicted AND the averaged posterior summary
#STEP 5:Fitting. NLMIXED_1 for one outcome, NLMIXED_MULTITRAJ for the joint two-outcome MODEL
#STEP 6:Built-IN demo pipeline, RUN only WHEN TRAJ2_RUN_DEMO=1
##########################################################################################################################;


/*##########################################################################################################################
*STEP 0: GLOBAL SETTINGS (EDIT THESE ONLY)

  TRAJ2_RUN_DEMO controls whether sourcing this file also RUNS the
  built-in demo (the Step 2 build plus the Step 4 pipeline). The default
  is 0, i.e. define macros only, which is the behaviour Section 5.1 of
  the paper describes. Set it to 1 BEFORE %include to get the old
  auto-run behaviour back.
##########################################################################################################################*/
%MACRO _traj2_demo_default;
%GLOBAL TRAJ2_RUN_DEMO;
%IF %LENGTH(%SUPERQ(TRAJ2_RUN_DEMO)) = 0 %THEN %LET TRAJ2_RUN_DEMO = 0;
%MEND _traj2_demo_default;
%_traj2_demo_default;

%LET USE_SIM     = 1;            /* 1 = run simulator, 0 = use your real BASE_FILE_SRS */
%LET T           = 12;           /* number of time points */
%LET class       = 5;            /* number of latent classes */
%LET order_model = 3;            /* 1=linear, 2=quadratic, 3=cubic */
%LET equal       = T;            /* T=equal sigma across classes (per outcome), F=class-specific */

/* censoring caps, numeric list length T */
%LET max_values = 10 10 10 10 10 10 10 10 10 10 10 10;

/* outcome variable lists in YOUR wide dataset */
%LET y1vars = QHH1-QHH12;        /* Outcome 1 repeated measures */
%LET y2vars = QINP1-QINP12;      /* Outcome 2 repeated measures */

/* ID + time index variables */
%LET idvar = BENE_ID;
%LET tvars = quar1-quar12;


/*##########################################################################################################################
*STEP 1: OPTIONAL SIMULATOR (WIDE + LONG)
##########################################################################################################################*/
%MACRO sim_data_cont(
class=5,
N=500,
T=12,
seed=2026,
  miss_pattern=balanced,   /* balanced | unbalanced */
p_obs_min=0.6,
max_const=10,
sigma1=1.0,
sigma2=1.2
);
DATA sim_long;
CALL STREAMINIT(&seed);

DO ID = 1 TO &n;

class = CEIL(RAND('uniform') * &class);

IF "&miss_pattern" = "balanced" THEN DO;
first_t = 1; last_t = &T;
END;
ELSE DO;
frac    = RAND('uniform')*(1-&p_obs_min) + &p_obs_min;
n_obs   = CEIL(&T*frac);
first_t = 1; last_t = n_obs;
END;

SELECT (class);
WHEN (1) DO; b10=-2.0; b11= 0.25; b12= 0.00;  b20=-1.0; b21= 0.10; b22=0.01; END;
WHEN (2) DO; b10=-1.0; b11= 0.15; b12= 0.01;  b20=-1.8; b21= 0.22; b22=0.00; END;
WHEN (3) DO; b10=-0.5; b11= 0.05; b12= 0.02;  b20=-0.7; b21= 0.08; b22=0.02; END;
WHEN (4) DO; b10=-2.5; b11= 0.35; b12=-0.01;  b20=-1.2; b21= 0.18; b22=0.00; END;
OTHERWISE DO; b10=-1.5; b11= 0.12; b12= 0.00;  b20=-1.3; b21= 0.12; b22=0.01; END;
END;

DO qtr = 1 TO &T;
obs   = (qtr >= first_t AND qtr <= last_t);
cap_t = &max_const;

y1 = .; y2 = .;

IF obs THEN DO;
mu1 = b10 + b11*qtr + b12*(qtr*qtr);
mu2 = b20 + b21*qtr + b22*(qtr*qtr);

z1 = mu1 + RAND('normal', 0, &sigma1);
z2 = mu2 + RAND('normal', 0, &sigma2);

y1 = MAX(0, MIN(cap_t, z1));
y2 = MAX(0, MIN(cap_t, z2));
END;

OUTPUT;
END;

END;
KEEP ID class qtr y1 y2 obs;
RUN;

PROC SORT DATA=sim_long; BY ID qtr; RUN;

PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=_y1_w PREFIX=Y1_;
BY ID; ID qtr; VAR y1;
RUN;

PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=_y2_w PREFIX=Y2_;
BY ID; ID qtr; VAR y2;
RUN;

DATA sim_wide;
MERGE _y1_w _y2_w;
BY ID;
RUN;
%MEND sim_data_cont;


/*##########################################################################################################################
*STEP 2: BUILD / LOAD BASE_FILE_SRS
  + adds CAP1-CAP&T as real variables (important fix)
##########################################################################################################################*/
%MACRO add_caps(ds=BASE_FILE_SRS);
%LOCAL _i;
DATA &ds;
SET &ds;
ARRAY CAP[&T] CAP1-CAP&T;
    /* load from macro list */
    /* _i below is a MACRO index, not a DATA step variable, so there is
       nothing to DROP. The original "drop _i;" referred to a variable
       that is never created in this step. */
%DO _i=1 %TO &T;
CAP[&_i] = %SCAN(&max_values, &_i);
%END;
RUN;
%MEND add_caps;

%MACRO build_base_file_srs;
%IF &USE_SIM = 1 %THEN %DO;

%sim_data_cont(class=&class, N=500, T=&T, seed=1, miss_pattern=balanced, max_const=10);

DATA BASE_FILE_SRS;
SET sim_wide;
RENAME
ID   = &idvar
Y1_1 = QHH1   Y1_2 = QHH2   Y1_3 = QHH3   Y1_4 = QHH4   Y1_5 = QHH5   Y1_6 = QHH6
Y1_7 = QHH7   Y1_8 = QHH8   Y1_9 = QHH9   Y1_10= QHH10  Y1_11= QHH11  Y1_12= QHH12
Y2_1 = QINP1  Y2_2 = QINP2  Y2_3 = QINP3  Y2_4 = QINP4  Y2_5 = QINP5  Y2_6 = QINP6
Y2_7 = QINP7  Y2_8 = QINP8  Y2_9 = QINP9  Y2_10= QINP10 Y2_11= QINP11 Y2_12= QINP12
;
RUN;

DATA BASE_FILE_SRS;
SET BASE_FILE_SRS;
ARRAY quar[&T] &tvars;
DO _i=1 TO &T;
quar[_i] = _i;
END;
DROP _i;
RUN;

%add_caps(ds=BASE_FILE_SRS);

%END;
%ELSE %DO;
%PUT NOTE: USE_SIM=0, expecting BASE_FILE_SRS already exists with ID=&idvar, time=&tvars, y1=&y1vars, y2=&y2vars.;
%add_caps(ds=BASE_FILE_SRS);
%END;
%MEND build_base_file_srs;

/* wrapped in a macro so it runs on SAS 9.4 releases before M5, where
   %IF is not allowed in open code */
%MACRO _traj2_demo_build;
%IF &TRAJ2_RUN_DEMO = 1 %THEN %build_base_file_srs;
%MEND _traj2_demo_build;
%_traj2_demo_build;


/*##########################################################################################################################
*STEP 3: MACRO LIBRARY
##########################################################################################################################*/
%LET class_names = A B C D E F G H I J K L M N O P Q R S T;

%MACRO starting_value_alpha(class);
%LOCAL s class_;
%DO s=2 %TO &class;
%LET class_ = %SCAN(&class_names, &s);
alpha0_&class_.=0
%END;
%MEND starting_value_alpha;

/* Starting values for the polynomial coefficients and residual scale.

   int_lo=, int_hi= and sigma0= are OPTIONAL and backward compatible. If
   int_lo= and int_hi= are left blank the macro reproduces the original
   behaviour exactly: every class intercept 0, every sigma 30.

   WARNING: leaving them blank gives every latent class an IDENTICAL
   starting point. That is an exchangeable, symmetric start for a finite
   mixture, and the classes have no gradient information to separate them.
   Supplying int_lo= and int_hi= spreads the class intercepts evenly over
   the plausible outcome range, which is what the validation runs use.
   Example for an outcome censored to [0,100] with residual SD near 6:
     %starting_value_beta_sigma(class=3, outcome=1, order=2, equal_sigma=T,
                                int_lo=10, int_hi=45, sigma0=6)           */
%MACRO starting_value_beta_sigma(class,outcome,ORDER,equal_sigma,
int_lo=,int_hi=,sigma0=30);
%LOCAL s i class_ v_;

%IF %UPCASE(&equal_sigma)=T %THEN %DO;
sigma&outcome._ = &sigma0
%END;
%ELSE %DO;
sigma&outcome._A = &sigma0
%DO s=2 %TO &class;
%LET class_=%SCAN(&class_names, &s);
sigma&outcome._&class_. = &sigma0
%END;
%END;

%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names, &s);

%IF %LENGTH(%SUPERQ(int_lo)) = 0 OR %LENGTH(%SUPERQ(int_hi)) = 0 %THEN %LET v_ = 0;
%ELSE %IF &class = 1 %THEN %LET v_ = %SYSEVALF((&int_lo + &int_hi)/2);
%ELSE %LET v_ = %SYSEVALF(&int_lo + (&int_hi - &int_lo)*(&s - 1)/(&class - 1));

beta&outcome._&class_.0 = &v_
%DO i=1 %TO &order;
beta&outcome._&class_.&i = 0
%END;
%END;
%MEND starting_value_beta_sigma;

%MACRO bounds_alpha(bounds_alpha,class);
%LOCAL s class_;
%DO s=2 %TO &class;
%LET class_=%SCAN(&class_names, &s);
-&bounds_alpha.<alpha0_&class_.<&bounds_alpha.
%IF &s < &class %THEN ,;
%END;
%MEND bounds_alpha;

%MACRO bounds_sigma(bounds_sigma,class,outcome,equal_sigma);
%LOCAL s class_;
%IF %UPCASE(&equal_sigma)=T %THEN %DO;
sigma&outcome._ > &bounds_sigma.
%END;
%ELSE %DO;
sigma&outcome._A > &bounds_sigma.
%DO s=2 %TO &class;
%LET class_=%SCAN(&class_names, &s);
, sigma&outcome._&class_. > &bounds_sigma.
%END;
%END;
%MEND bounds_sigma;

%MACRO initiation_universal(T);
ARRAY X[&T] &tvars;
ARRAY Y1[&T] &y1vars;
ARRAY Y2[&T] &y2vars;
%MEND initiation_universal;

%MACRO initiation(T,class,outcome);
%LOCAL s class_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
ARRAY PI&class_.&outcome.[&T] PI&class_.&outcome._1-PI&class_.&outcome._&T;
ARRAY mu&class_.&outcome.[&T] mu&class_.&outcome._1-mu&class_.&outcome._&T;
PROD&class_.&outcome.=0;
%END;
%MEND initiation;

%MACRO MODEL(ORDER,class,outcome);
%LOCAL s i class_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
mu&class_.&outcome.[I] = beta&outcome._&class_.0
%DO i=1 %TO &order;
+ beta&outcome._&class_.&i * (X[I]**&i)
%END;
;
%END;
%MEND MODEL;

%MACRO residual(class,outcome);
%LOCAL s class_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
e&outcome._&class_. = Y&outcome.[I] - mu&class_.&outcome.[I];
%END;
%MEND residual;

/* Numerical guard. The raw residual is clipped to +/- &float residual SDs
   BEFORE the density is evaluated. This is what stops exp(PROD) from
   underflowing to zero for classes that sit far from an observation, which
   would make l_latclass = 0 and log(0) missing. It is a deviation from the
   exact likelihood written in Section 2.3 of the paper and should be stated
   there. A log-sum-exp evaluation of l_latclass would remove the need for
   it; see the ZIP prototype for that pattern.                              */
%MACRO float_control(float,class,outcome,equal_sigma);
%LOCAL s class_ class2_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
%IF %UPCASE(&equal_sigma)=T %THEN %LET class2_=;
%ELSE %LET class2_=&class_;
e&outcome._&class_. = MIN(MAX(e&outcome._&class_.,
-&float.*sigma&outcome._&class2_.),
&float.*sigma&outcome._&class2_.);
%END;
%MEND float_control;

/* censored-normal contribution uses CAP[I] (CAP array points to CAP1-CAP&T) */
%MACRO Prob_cnorm(class,outcome,equal_sigma);
%LOCAL s class_ class2_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
%IF %UPCASE(&equal_sigma)=T %THEN %LET class2_=;
%ELSE %LET class2_=&class_;

    /* A missing time point must contribute nothing to the likelihood.
       Without this guard PROD goes missing for any subject with a gap,
       which silently removes that subject from the fit. Table 1 of the
       paper promises the contribution is skipped, and miss_pattern=
       unbalanced in %sim_data_cont produces exactly this situation.
       For complete data the result is bit-identical to the old code.   */
IF Y&outcome.[I] = . THEN PI&class_.&outcome.[I] = 0;
ELSE DO;
PI&class_.&outcome.[I] = LOGPDF('NORMAL', e&outcome._&class_. / sigma&outcome._&class2_.)
- LOG(sigma&outcome._&class2_.);

IF Y&outcome.[I] = 0 THEN
PI&class_.&outcome.[I] = LOGCDF('NORMAL', e&outcome._&class_. / sigma&outcome._&class2_.);
ELSE IF Y&outcome.[I] = CAP[I] THEN
PI&class_.&outcome.[I] = LOGCDF('NORMAL', (-e&outcome._&class_.) / sigma&outcome._&class2_.);
END;

PROD&class_.&outcome. = PROD&class_.&outcome. + PI&class_.&outcome.[I];
%END;
%MEND Prob_cnorm;

%MACRO Class_Membership(class);
%LOCAL s class_;
alpha0_A = 0;

%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
pinumer_&class_. = EXP(alpha0_&class_.);
%END;

pideno = 0
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
+ pinumer_&class_.
%END;
;

%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
pie_&class_. = pinumer_&class_. / pideno;
%END;
%MEND Class_Membership;

%MACRO LogLike(class,outcome);
%LOCAL s class_;
l_latclass =
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
%IF &s>1 %THEN + ;
pie_&class_. * EXP(PROD&class_.&outcome.)
%END;
;
%MEND LogLike;

%MACRO LogLike_multi(class);
%LOCAL s class_;
( %DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
%IF &s>1 %THEN + ;
pie_&class_. * EXP(PROD&class_.1) * EXP(PROD&class_.2)
%END;
)
%MEND LogLike_multi;

%MACRO Posterior(class);
%LOCAL s class_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
post_&class_. = pie_&class_. * EXP(PROD&class_.1) * EXP(PROD&class_.2) / %LogLike_multi(class=&class);
%END;

KEEP
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
post_&class_.
%END;
;
%MEND Posterior;

%MACRO initiation_pred(T,class,outcome);
%LOCAL s class_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
ARRAY Pred&class_.&outcome.[&T] Pred&class_.&outcome._1-Pred&class_.&outcome._&T;
ARRAY mu&class_.&outcome.[&T]   mu&class_.&outcome._1-mu&class_.&outcome._&T;
%END;
%MEND initiation_pred;

%MACRO Pred_cnorm(class,outcome,equal_sigma);
%LOCAL s class_ class2_;
%DO s=1 %TO &class;
%LET class_=%SCAN(&class_names,&s);
%IF %UPCASE(&equal_sigma)=T %THEN %LET class2_=;
%ELSE %LET class2_=&class_;

temp = (EXP(LOGCDF('normal',(CAP[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)) -
EXP(LOGCDF('normal',(0      -mu&class_.&outcome.[I])/sigma&outcome._&class2_.)));
IF temp = 0 THEN temp = 0.000001;

Pred&class_.&outcome.[I] =
0
+ temp * (
mu&class_.&outcome.[I]
+ sigma&outcome._&class2_. *
(EXP(LOGPDF('normal', -mu&class_.&outcome.[I]/sigma&outcome._&class2_.)) -
EXP(LOGPDF('normal', (CAP[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.))) / temp
)
+ CAP[I] * EXP(LOGCDF('normal', (-CAP[I]+mu&class_.&outcome.[I])/sigma&outcome._&class2_.))
;
%END;
%MEND Pred_cnorm;


/*##########################################################################################################################
*STEP 4: PLOTTING
##########################################################################################################################*/
%MACRO plot_prep(T,LC,result,ORDER,equal_sigma);

DATA parameter; SET &result; KEEP parameter ESTIMATE; RUN;
PROC TRANSPOSE DATA=parameter OUT=parameter; ID parameter; RUN;

PROC SQL;
CREATE TABLE data_pred AS
SELECT * FROM BASE_FILE_SRS, parameter;
QUIT;

DATA pred_membership_y;
SET data_pred;

%initiation_universal(T=&T);
ARRAY CAP[&T] CAP1-CAP&T;

%initiation(T=&T, class=&LC, outcome=1);
%initiation(T=&T, class=&LC, outcome=2);

DO I=1 TO &T;
%MODEL(ORDER=&order, class=&LC, outcome=1);
%residual(class=&LC, outcome=1);
%float_control(float=8, class=&LC, outcome=1, equal_sigma=&equal_sigma);
%Prob_cnorm(class=&LC, outcome=1, equal_sigma=&equal_sigma);

%MODEL(ORDER=&order, class=&LC, outcome=2);
%residual(class=&LC, outcome=2);
%float_control(float=8, class=&LC, outcome=2, equal_sigma=&equal_sigma);
%Prob_cnorm(class=&LC, outcome=2, equal_sigma=&equal_sigma);
END;

%Class_Membership(class=&LC);
%Posterior(class=&LC);
RUN;

DATA y_1; SET BASE_FILE_SRS(KEEP=&y1vars); RUN;
PROC iml;
start colsum(m); RETURN(m[+,]); finish;
use y_1; read all VAR _num_ INTO y;
use pred_membership_y; read all VAR _num_ INTO p;
SUM = colsum(p);
invsum = 1/SUM;
avg_y  = y`*p;
avg_y2 = avg_y#invsum;
CREATE avg_y1 FROM avg_y2; append FROM avg_y2; CLOSE avg_y1;
QUIT;

DATA y_2; SET BASE_FILE_SRS(KEEP=&y2vars); RUN;
PROC iml;
start colsum(m); RETURN(m[+,]); finish;
use y_2; read all VAR _num_ INTO y;
use pred_membership_y; read all VAR _num_ INTO p;
SUM = colsum(p);
invsum = 1/SUM;
avg_y  = y`*p;
avg_y2 = avg_y#invsum;
CREATE avg_y2 FROM avg_y2; append FROM avg_y2; CLOSE avg_y2;
QUIT;

  /* build single-row dataset for prediction curves */
DATA wide;
SET parameter;
%DO j=1 %TO &T;
quar&j = &j;
CAP&j  = %SCAN(&max_values,&j);
%END;
RUN;

DATA wide;
SET wide;
ARRAY X[&T] &tvars;
ARRAY CAP[&T] CAP1-CAP&T;

%initiation_pred(T=&T, class=&LC, outcome=1);
%initiation_pred(T=&T, class=&LC, outcome=2);

DO I=1 TO &T;
%MODEL(ORDER=&order, class=&LC, outcome=1);
%MODEL(ORDER=&order, class=&LC, outcome=2);
%Pred_cnorm(class=&LC, outcome=1, equal_sigma=&equal_sigma);
%Pred_cnorm(class=&LC, outcome=2, equal_sigma=&equal_sigma);
END;
RUN;

PROC TRANSPOSE DATA=wide OUT=pred_temp name=_NAME_;
VAR _numeric_;
RUN;

DATA pred_temp;
SET pred_temp;
RENAME COL1 = ESTIMATE;
RUN;

DATA pred;
SET pred_temp;
LENGTH outcome $3 type $4 class $1;
WHERE _NAME_ contains 'Pred';
class  = UPCASE(SUBSTR(_NAME_,5,1));
y      = SUBSTR(_NAME_,6,1);
IF y='1' THEN outcome='HH';
IF y='2' THEN outcome='INP';
type   = 'pred';
quar   = INPUT(SUBSTR(_NAME_,8), best.);
KEEP outcome type class quar ESTIMATE;
RUN;

TITLE 'Predicted Trajectory by Latent Class';
PROC SGPANEL DATA=pred;
PANELBY class / columns=4;
SERIES x=quar y=ESTIMATE / GROUP=outcome;
RUN;
TITLE;

DATA avg_y1_plot;
SET avg_y1;
LENGTH outcome $3 type $3 class $1;
ARRAY col col1-col&LC;
quar=_n_;
DO i=1 TO &LC;
ESTIMATE=col[i];
class=UPCASE(byte(i+64));
outcome="HH"; type="Avg";
OUTPUT;
END;
KEEP ESTIMATE class quar outcome type;
RUN;

DATA avg_y2_plot;
SET avg_y2;
LENGTH outcome $3 type $3 class $1;
ARRAY col col1-col&LC;
quar=_n_;
DO i=1 TO &LC;
ESTIMATE=col[i];
class=UPCASE(byte(i+64));
outcome="INP"; type="Avg";
OUTPUT;
END;
KEEP ESTIMATE class quar outcome type;
RUN;

DATA avg; SET avg_y1_plot avg_y2_plot; RUN;

PROC SORT DATA=avg;  BY outcome quar class; RUN;
PROC SORT DATA=pred; BY outcome quar class; RUN;

DATA pred2;
SET pred;
pred = ESTIMATE;
KEEP outcome quar class pred;
RUN;

DATA data_plot;
MERGE avg(RENAME=(ESTIMATE=avg)) pred2;
BY outcome quar class;
RUN;

TITLE 'Predicted vs Averaged Observed by Latent Class';
PROC SGPANEL DATA=data_plot;
PANELBY outcome;
SERIES  x=quar y=pred / GROUP=class name="pred";
SCATTER x=quar y=avg  / GROUP=class name="obs";
KEYLEGEND "pred" / TITLE="Predicted";
KEYLEGEND "obs"  / TITLE="Averaged Observed";
RUN;
TITLE;

TITLE 'Averaged Posterior Class Membership';
PROC MEANS DATA=pred_membership_y MEAN MISSING;
VAR
%DO s=1 %TO &LC;
%LET class_=%SCAN(&class_names,&s);
post_&class_.
%END;
;
OUTPUT OUT=avg_membership MEAN=Avg STD=SD;
RUN;
TITLE;

%MEND plot_prep;


/*##########################################################################################################################
*STEP 5: NLMIXED FITTING MACROS
##########################################################################################################################*/
%MACRO nlmixed_1(T,LC,Y,starting,OUTPUT,ORDER,equal_sigma);

PROC NLMIXED DATA=BASE_FILE_SRS itdetails qpoints=40 noad maxiter=1000 tech=dbldog;

BOUNDS %bounds_alpha(bounds_alpha=3, class=&LC),
%bounds_sigma(bounds_sigma=0, class=&LC, outcome=&Y, equal_sigma=&equal_sigma);

PARMS &starting;

%initiation_universal(T=&T);
ARRAY CAP[&T] CAP1-CAP&T;

%initiation(T=&T, class=&LC, outcome=&Y);

DO I=1 TO &T;
%MODEL(ORDER=&order, class=&LC, outcome=&Y);
%residual(class=&LC, outcome=&Y);
%float_control(float=8, class=&LC, outcome=&Y, equal_sigma=&equal_sigma);
%Prob_cnorm(class=&LC, outcome=&Y, equal_sigma=&equal_sigma);
END;

%Class_Membership(class=&LC);
%LogLike(class=&LC, outcome=&Y);

ll_latclass = LOG(l_latclass);
MODEL ll_latclass ~ GENERAL(ll_latclass);

ODS OUTPUT ParameterEstimates=work.&output;
RUN;

%MEND nlmixed_1;

%MACRO nlmixed_MultiTraj(T,LC,starting,OUTPUT,ORDER,equal_sigma);

PROC NLMIXED DATA=BASE_FILE_SRS itdetails qpoints=40 noad maxiter=1000 tech=dbldog;

BOUNDS %bounds_alpha(bounds_alpha=3, class=&LC),
%bounds_sigma(bounds_sigma=0, class=&LC, outcome=1, equal_sigma=&equal_sigma),
%bounds_sigma(bounds_sigma=0, class=&LC, outcome=2, equal_sigma=&equal_sigma);

PARMS &starting;

%initiation_universal(T=&T);
ARRAY CAP[&T] CAP1-CAP&T;

%initiation(T=&T, class=&LC, outcome=1);
%initiation(T=&T, class=&LC, outcome=2);

DO I=1 TO &T;
%MODEL(ORDER=&order, class=&LC, outcome=1);
%residual(class=&LC, outcome=1);
%float_control(float=8, class=&LC, outcome=1, equal_sigma=&equal_sigma);
%Prob_cnorm(class=&LC, outcome=1, equal_sigma=&equal_sigma);

%MODEL(ORDER=&order, class=&LC, outcome=2);
%residual(class=&LC, outcome=2);
%float_control(float=8, class=&LC, outcome=2, equal_sigma=&equal_sigma);
%Prob_cnorm(class=&LC, outcome=2, equal_sigma=&equal_sigma);
END;

%Class_Membership(class=&LC);

l_latclass  = %LogLike_multi(class=&LC);
ll_latclass = LOG(l_latclass);
MODEL ll_latclass ~ GENERAL(ll_latclass);

ODS OUTPUT ParameterEstimates=work.&output;
RUN;

%MEND nlmixed_MultiTraj;


/*##########################################################################################################################
*STEP 6: DEMO PIPELINE, SKIPPED UNLESS TRAJ2_RUN_DEMO=1
##########################################################################################################################*/
%MACRO _traj2_demo_pipeline;
%IF &TRAJ2_RUN_DEMO NE 1 %THEN %DO;
%PUT NOTE: Traj2 continuous macros loaded. Demo pipeline skipped (TRAJ2_RUN_DEMO=&TRAJ2_RUN_DEMO).;
%RETURN;
%END;

%nlmixed_1(
T=&T, LC=&class, Y=1,
starting=%starting_value_alpha(class=&class)
%starting_value_beta_sigma(class=&class,outcome=1,ORDER=&order_model,equal_sigma=&equal),
OUTPUT=nlm_fix_T1&class,
ORDER=&order_model,
equal_sigma=&equal
);

%nlmixed_1(
T=&T, LC=&class, Y=2,
starting=%starting_value_alpha(class=&class)
%starting_value_beta_sigma(class=&class,outcome=2,ORDER=&order_model,equal_sigma=&equal),
OUTPUT=nlm_fix_T2&class,
ORDER=&order_model,
equal_sigma=&equal
);

DATA work.nlm_2y_starting;
SET nlm_fix_T1&class nlm_fix_T2&class;
IF parameter =: 'alpha' THEN DELETE;
RUN;

%nlmixed_MultiTraj(
T=&T, LC=&class,
starting=%starting_value_alpha(class=&class) / DATA=work.nlm_2y_starting,
OUTPUT=nlm_fix_T1_T2&class,
ORDER=&order_model,
equal_sigma=&equal
);

%plot_prep(
T=&T, LC=&class,
result=nlm_fix_T1_T2&class,
ORDER=&order_model,
equal_sigma=&equal
);

%MEND _traj2_demo_pipeline;
%_traj2_demo_pipeline;

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT DATASETS
sim_wide                wide simulated outcomes, one row per subject, from SIM_DATA_CONT
BASE_FILE_SRS           the fitting contract. Outcomes, quar1 to quarT and CAP1 to CAPT
work.<output>           ParameterEstimates from NLMIXED_1 or NLMIXED_MULTITRAJ, named by the OUTPUT argument
pred_membership_y       one row per subject with predicted values and posterior class probabilities
avg_membership          averaged posterior class membership by assigned class

TYPICAL WORKFLOW
  %let TRAJ2_RUN_DEMO = 0;
  %include "continuous_2outcomes.sas";

  %nlmixed_1(T=12, LC=3, Y=1, starting=<parms>, output=fit_y1, order=2, equal_sigma=T);
  %nlmixed_1(T=12, LC=3, Y=2, starting=<parms>, output=fit_y2, order=2, equal_sigma=T);
  %nlmixed_MultiTraj(T=12, LC=3, starting=<alphas> / data=<stacked>, output=fit_joint,
                     order=2, equal_sigma=T);
  %plot_prep(T=12, LC=3, result=fit_joint, order=2, equal_sigma=T);

NOTE ON STARTING VALUES
  Leaving int_lo= and int_hi= blank in STARTING_VALUE_BETA_SIGMA gives every latent class an identical
  starting point, which is an exchangeable start for a finite mixture. Supply int_lo= and int_hi= to
  spread the class intercepts over the plausible outcome range, as the validation runs do.
##########################################################################################################################*/
