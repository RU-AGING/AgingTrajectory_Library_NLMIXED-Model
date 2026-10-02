/*Copyright (c) 2026 Community Health and Aging Outcomes (CHAO) Lab, Rutgers University.
Released under the MIT License. Full text in LICENSE at
https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model;
*PROJECT NAME: Traj2 Ordinal-Probit Group-Based Trajectory Macros
LAST UPDATED DATE: 01 OCT 2026
DATA SOURCES: NONE. This file defines macros only. Input is either the optional simulator or a user table
PURPOSE: Single-outcome ordinal-probit latent class trajectory modelling in base SAS. Provides
(a)an optional data simulator for worked examples
(b)a wide-format data-prep helper
(c)one PROC NLMIXED fit at a fixed number of latent classes
(d)posterior class membership, modal assignment, class sizes with APP and OCC
(e)class-proportion and predicted-trajectory plots
IDENTIFICATION: the first threshold is fixed at 0 and the class intercept beta_*0 is freely estimated.
The free threshold parameters are the increments i*2, i*3 and onward. There is no i*1. The increments are
held positive by a BOUNDS statement, so the thresholds stay strictly ordered within each class.
AUTHOR: Haiqun Lin, Weiyi Xia, Anum Zafar
##########################################################################################################################
*Execution Environment: SAS 9.4 or later, PROC NLMIXED from SAS/STAT. No compiled components         *
*This file defines macros only. Including it runs nothing                                            *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 1:Optional simulator. SIM_DATA builds SIM_WIDE
#STEP 2:Data prep. BUILD_BASE_FROM_SIMWIDE builds BASE_FILE_SRS with quar1 to quarT holding 1 to T
#STEP 3:Helper macros. Class letters, override indexing, default PARMS construction
#STEP 4:Model fit. ORDPROB_MIX_FIT_ONE runs one NLMIXED fit at a fixed class count
#STEP 5:Post-estimation. ORDPROB_MIX_POST_ONE computes posteriors, modal class, class sizes, APP and OCC
#STEP 6:Plots. ORDPROB_MIX_PLOT_ONE draws class proportions and predicted mean trajectories
##########################################################################################################################;
*/
/*##########################################################################################################################
*STEP 1: OPTIONAL SIMULATOR
##########################################################################################################################
Builds SIM_WIDE for the worked examples. DIST selects the generating process. BIN4 and TPOIS4 give each
class a distinct mean trajectory, CAT5 draws categories uniformly. SIGMA_RE adds a subject-level random
effect, so set SIGMA_RE=0 when the fitted model is to assume conditional independence given class.

SAMPLE DATA OUTPUT DATASET (SIM_WIDE):
id   Y1_1   Y1_2   ...   Y1_12   Y2_1   Y2_2   ...   Y2_12
1    0      1            3       1      1            2
##########################################################################################################################*/
%MACRO sim_data(
  dist=BIN4, class=4, n=2000, T=12, seed=1,
  miss_pattern=balanced, p_obs_min=0.6, sigma_re=0.4, cap=4
);
%LOCAL _DIST;
%LET _DIST=%SYSFUNC(UPCASE(%SYSFUNC(STRIP(&dist.))));
DATA sim_long;
CALL STREAMINIT(&seed.);
DO id = 1 TO &n.;
class = CEIL(RAND('uniform') * &class.);
IF "&miss_pattern."="balanced" THEN DO; first_t=1; last_t=&T.; END;
ELSE DO; frac=RAND('uniform')*(1-&p_obs_min.) + &p_obs_min.; n_obs=CEIL(&T.*frac); first_t=1; last_t=n_obs; END;
bi = RAND('normal', 0, &sigma_re.);
DO t = 1 TO &T.;
qtr=t;
obs=(t>=first_t AND t<=last_t);
y1=.;
y2=.;
IF obs THEN DO;
SELECT (class);
WHEN (1) DO; mu1=0.6+0.06*t;          mu2=0.8+0.03*t; END;
WHEN (2) DO; mu1=0.5+0.30*t;          mu2=0.7+0.10*t; END;
WHEN (3) DO; mu1=1.0+0.10*t+0.01*t*t; mu2=0.9+0.08*t; END;
OTHERWISE DO; mu1=0.5+0.20*t;         mu2=0.5+0.15*t; END;
END;
%IF "&_DIST."="BIN4" %THEN %DO;
m1=MAX(0,mu1*EXP(bi)); m2=MAX(0,mu2*EXP(bi));
p1=MIN(MAX(m1/&cap.,0),1); p2=MIN(MAX(m2/&cap.,0),1);
y1=RAND('binomial',p1,&cap.); y2=RAND('binomial',p2,&cap.);
%END;
%ELSE %IF "&_DIST."="TPOIS4" %THEN %DO;
m1=MAX(1e-6,mu1*EXP(bi)); m2=MAX(1e-6,mu2*EXP(bi));
ARRAY p1[%EVAL(&cap.+1)]; ARRAY p2[%EVAL(&cap.+1)];
s1=0; s2=0;
DO k=0 TO &cap.; p1[k+1]=PDF('poisson',k,m1); s1 + p1[k+1];
                 p2[k+1]=PDF('poisson',k,m2); s2 + p2[k+1]; END;
DO k=0 TO &cap.; p1[k+1]=p1[k+1]/s1; p2[k+1]=p2[k+1]/s2; END;
y1=RAND('table', OF p1[*]) - 1; y2=RAND('table', OF p2[*]) - 1;
%END;
%ELSE %IF "&_DIST."="CAT5" %THEN %DO;
y1=RAND('integer',%EVAL(&cap.+1))-1; y2=RAND('integer',%EVAL(&cap.+1))-1;
%END;
%ELSE %DO;
IF id=1 AND t=1 THEN PUT "WARNING: dist=&dist not recognized. Using BIN4.";
m1=MAX(0,mu1*EXP(bi)); m2=MAX(0,mu2*EXP(bi));
p1=MIN(MAX(m1/&cap.,0),1); p2=MIN(MAX(m2/&cap.,0),1);
y1=RAND('binomial',p1,&cap.); y2=RAND('binomial',p2,&cap.);
%END;
END;
OUTPUT;
END;
END;
KEEP id class qtr y1 y2 obs;
RUN;

PROC SORT DATA=sim_long; BY id qtr; RUN;
PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=y1_w PREFIX=Y1_; BY id; ID qtr; VAR y1; RUN;
PROC TRANSPOSE DATA=sim_long(WHERE=(obs=1)) OUT=y2_w PREFIX=Y2_; BY id; ID qtr; VAR y2; RUN;
DATA sim_wide; MERGE y1_w y2_w; BY id; RUN;
%MEND sim_data;

/*##########################################################################################################################
*STEP 2: DATA PREP
##########################################################################################################################
Turns SIM_WIDE into the wide-format contract the fitting macro expects. The time index columns quar1 to
quarT hold the numeric values 1 to T. To use your own data instead, build a table with the same shape and
pass it through data= on the fitting macro.

SAMPLE DATA OUTPUT DATASET (BASE_FILE_SRS):
BENE_ID   Y1_1   ...   Y1_12   quar1   quar2   ...   quar12
1         0            3       1       2             12
##########################################################################################################################*/
%MACRO build_base_from_simwide(T=12);
%IF %SYSFUNC(EXIST(work.sim_wide)) %THEN %DO;
DATA BASE_FILE_SRS;
SET sim_wide;
ARRAY quar[&T] quar1-quar&T;
DO i=1 TO &T; quar[i]=i; END;
DROP i;
RENAME id = BENE_ID;
RUN;
%END;
%ELSE %PUT NOTE: SIM_WIDE not found. Point DATA= in ordprob_mix_fit_one() to your own table.;
%MEND build_base_from_simwide;

/*##########################################################################################################################
*STEP 3: HELPER MACROS
##########################################################################################################################
CL maps a class index to its letter. _INDEX_OVERRIDE_PAIRS flags parameters the user supplied through
start_values so the defaults do not emit them twice. _MAKE_DEFAULT_PARMS builds the PARMS list.

The PARMS list holds the mixing intercepts alpha0_* for classes 2 upward, the polynomial coefficients
beta_*0 to beta_*deg, and the threshold increments i*2 to i*(m-1). The first threshold is fixed at 0 and
is therefore NOT a parameter, so i*1 is never emitted.

The increments default to 1 rather than 0 because they are bounded below by zero in STEP 4 and a start of
zero would sit on the boundary. Starting every increment at 1 places the first threshold at 0, the second
at 1, the third at 2 and so on, which is the same starting point in threshold space that the earlier
exponential-increment version reached from a start of zero.
##########################################################################################################################*/
*Class index to letter, A through O;
%MACRO CL(c); %SCAN(A B C D E F G H I J K L M N O, &c., %STR( )) %MEND CL;

*Mark user-specified overrides so the defaults skip them;
%MACRO _index_override_pairs(pairs);
%LOCAL i token eqpos name;
%DO i=1 %TO %SYSFUNC(COUNTW(&pairs,%STR( )));
%LET token=%SCAN(&pairs,&i,%STR( ));
%LET eqpos=%INDEX(&token,=);
%IF &eqpos>1 %THEN %DO;
%LET name=%SUBSTR(&token,1,%EVAL(&eqpos-1));
%GLOBAL OV_&name;
%LET OV_&name=1;
%END;
%END;
%MEND _index_override_pairs;

*Build the default PARMS list for (k, m, deg), skipping anything the user overrode;
%MACRO _make_default_parms(k, m, deg);
%LOCAL c L l j s m1;
%LET s=;
%LET m1=%EVAL(&m-1);
%DO c=1 %TO &k;
%LET L=%CL(&c); %LET l=%LOWCASE(&L);
%*mixing intercepts. Class A is fixed at 0 so it is skipped;
%IF &c>1 %THEN %IF NOT %SYMEXIST(OV_alpha0_&L) %THEN %LET s=&s alpha0_&L.=0;
%*polynomial coefficients. beta_*0 is free and carries the class location;
%IF NOT %SYMEXIST(OV_beta_&l.0) %THEN %LET s=&s beta_&l.0=0;
%IF &deg>=1 %THEN %DO; %IF NOT %SYMEXIST(OV_beta_&l.1) %THEN %LET s=&s beta_&l.1=0; %END;
%IF &deg>=2 %THEN %DO j=2 %TO &deg; %IF NOT %SYMEXIST(OV_beta_&l.&j) %THEN %LET s=&s beta_&l.&j.=0; %END;
%*threshold increments i*2 to i*(m-1). The first threshold is fixed at 0. Started at 1, not 0,
 because the increments are bounded below by zero in the fitting macro;
%IF &m1>=2 %THEN %DO j=2 %TO &m1; %IF NOT %SYMEXIST(OV_i&l.&j) %THEN %LET s=&s i&l.&j.=1; %END;
%END;
&s
%MEND _make_default_parms;

/*##########################################################################################################################
*STEP 4: MODEL FIT
##########################################################################################################################
One PROC NLMIXED run at a fixed class count k. Class A mixing intercept is fixed at 0, class weights come
from a softmax over alpha0_*, and each class contributes a cumulative-probit likelihood with thresholds
built as cumulative sums of increments, th1 = 0 and thj = th(j-1) + i*j, where every increment is held
positive by a generated BOUNDS statement so the thresholds stay strictly ordered within each class.

Time points whose outcome value does not match any entry in ycodes, including missing values, contribute
nothing to the likelihood. The subject still contributes every time point that does match.

SAMPLE DATA OUTPUT DATASETS (prefix and k as supplied):
ordprob_single_EST_K3    Parameter   Estimate   StandardError   Probt
ordprob_single_FIT_K3    Descr       Value
ordprob_single_ESTS_K3   Label       Estimate   StandardError
##########################################################################################################################*/
%MACRO ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=BENE_ID,
  yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
  m=4, ycodes=0 1 2 3, deg=2,
  k=4, qpoints=40, maxiter=1000, tech=dbldog,
  start_values=,                     /* e.g., %str(alpha0_B=-0.5 beta_a0=0.1 ia2=1 ia3=1) */
  bounds=,                           /* e.g., %str(-6 < beta_a0 beta_b0 beta_c0 beta_d0 < 6)    */
  outestlib=work, prefix=ordprob_single
);
%LOCAL m1 d j _ibnd _bc _bj _bl;
%LET m1=%EVAL(&m-1);
%LET d=&deg;

*Hard stop if the input table is missing;
%IF NOT %SYSFUNC(EXIST(&data)) %THEN %DO;
%PUT ERROR: DATA=&data not found. Build it with SIM_DATA and BUILD_BASE_FROM_SIMWIDE, or point DATA= at your own table.;
%RETURN;
%END;

*Clear override flags left over from an earlier call in this session, so a previous run start_values
 cannot silently drop parameters from this run PARMS statement;
%LOCAL _ovlist _ovi;
%LET _ovlist=;
PROC SQL NOPRINT;
SELECT name INTO :_ovlist SEPARATED BY ' '
FROM dictionary.macros
WHERE scope='GLOBAL' AND name LIKE 'OV_%';
QUIT;
%IF %LENGTH(&_ovlist) %THEN %DO _ovi=1 %TO %SYSFUNC(COUNTW(&_ovlist,%STR( )));
%SYMDEL %SCAN(&_ovlist,&_ovi,%STR( )) / NOWARN;
%END;

*Category codes to macro variables;
%DO j=1 %TO &m; %GLOBAL ycode&j; %LET ycode&j=%SCAN(&ycodes,&j,%STR( )); %END;

*Index the user overrides so the defaults skip them;
%IF %LENGTH(&start_values) %THEN %_index_override_pairs(%SUPERQ(start_values));

*Build the increment bounds i*2>0 to i*(m-1)>0 for every class;
%LET _ibnd=;
%DO _bc=1 %TO &k;
%LET _bl=%LOWCASE(%CL(&_bc));
%IF &m1>=2 %THEN %DO _bj=2 %TO &m1;
%IF %LENGTH(&_ibnd) %THEN %LET _ibnd=&_ibnd,;
%LET _ibnd=&_ibnd i&_bl.&_bj.>0;
%END;
%END;

ODS LISTING;
ODS OUTPUT
ParameterEstimates  =&outestlib..&prefix._EST_K&k
FitStatistics       =&outestlib..&prefix._FIT_K&k
AdditionalEstimates =&outestlib..&prefix._ESTS_K&k;

PROC NLMIXED DATA=&data QPOINTS=&qpoints NOAD MAXITER=&maxiter TECH=&tech;
ARRAY yv[&ttotal] &yvars;
ARRAY tv[&ttotal] &tvars;

*ONE PARMS statement. Defaults with duplicates removed, then the user overrides;
PARMS %_make_default_parms(&k, &m, &deg)
      %SUPERQ(start_values);

*Threshold increments must stay positive so the thresholds remain ordered. Generated here and combined
 with anything the user supplies through bounds=;
%IF %LENGTH(&_ibnd) %THEN %DO;
BOUNDS &_ibnd
%IF %LENGTH(&bounds) %THEN %DO; , &bounds %END;
;
%END;
%ELSE %IF %LENGTH(&bounds) %THEN %DO;
BOUNDS &bounds;
%END;

*Fix the class A mixing intercept;
alpha0_A = 0;

*Softmax mixing weights;
denom=0;
%DO c=1 %TO &k; %LET U=%CL(&c); pin&c = EXP(alpha0_&U); denom + pin&c; %END;
%DO c=1 %TO &k; pie&c = pin&c/denom; %END;

*Class log-likelihoods;
%DO c=1 %TO &k; llik&c=0; %END;

DO t=1 TO &ttotal;
tval = tv[t];
cat = .; %DO j=1 %TO &m; IF yv[t] = &&ycode&j THEN cat=&j; %END;

*cat stays missing when the value is blank or is not one of the ycodes, so that time point is skipped;
IF cat>0 THEN DO;
%DO c=1 %TO &k;
%LET U=%CL(&c); %LET L=%LOWCASE(&U);

*eta(t). Polynomial in t;
eta_&U = beta_&L.0;
%IF &d>=1 %THEN %DO; eta_&U = eta_&U + beta_&L.1*(tval); %END;
%IF &d>=2 %THEN %DO j=2 %TO &d; eta_&U = eta_&U + beta_&L.&j.*((tval)**&j); %END;

*Ordered thresholds. th1 is fixed at 0, thj = th(j-1) + iLj with iLj held positive by BOUNDS;
th1_&U = 0;
%IF &m1>=2 %THEN %DO j=2 %TO &m1;
%IF &j=2 %THEN %DO; th&j._&U = th1_&U + i&L.&j; %END;
%ELSE %DO;          th&j._&U = th%EVAL(&j-1)_&U + i&L.&j; %END;
%END;

*Category probabilities;
%DO j=1 %TO &m;
%IF &j=1 %THEN %DO;
z&j._&U = PROBNORM(th1_&U - eta_&U);
%END; %ELSE %IF &j<&m %THEN %DO;
z&j._&U = PROBNORM(th&j._&U - eta_&U) - PROBNORM(th%EVAL(&j-1)_&U - eta_&U);
%END; %ELSE %DO;
z&j._&U = 1 - PROBNORM(th%EVAL(&m1)_&U - eta_&U);
%END;
%END;

ARRAY z&U[&m] %DO j=1 %TO &m; z&j._&U %END;;
p_obs&c = MAX(1e-12, z&U[cat]);
llik&c  = llik&c + LOG(p_obs&c);
%END;
END;
END;

mix=0; %DO c=1 %TO &k; mix + pie&c*EXP(llik&c); %END;
mix=MAX(mix,1e-300);
dummy=0;
MODEL dummy ~ GENERAL(LOG(mix));

*Thresholds on the natural scale. threshold1 is fixed at 0 and is not estimated, so report
 threshold2 to threshold(m-1) where threshold_jj = iL2 + ... + iL_jj;
%DO c=1 %TO &k; %LET U=%CL(&c); %LET L=%LOWCASE(&U);
%IF &m1>=2 %THEN %DO _jj=2 %TO &m1;
%LOCAL _expr; %LET _expr = i&L.2;
%DO j=3 %TO &_jj; %LET _expr = &_expr + i&L.&j; %END;
ESTIMATE "threshold&_jj._&L" &_expr;
%END;
%END;
RUN;

ODS OUTPUT CLOSE;
%MEND ordprob_mix_fit_one;

/*##########################################################################################################################
*STEP 5: POST-ESTIMATION
##########################################################################################################################
Posterior class-membership probabilities for every subject, modal class assignment, and the class table
with average posterior probability (APP) and odds of correct classification (OCC). Reads the parameter
estimates written by STEP 4 and recomputes each class log-likelihood on the same data with the same
likelihood and the same skipping rule for blanks, so the posteriors are consistent with the fit.

APP_j is the mean posterior probability of class j over the subjects assigned to class j. Nagin recommends
APP of at least 0.7. OCC_j = [APP_j/(1-APP_j)] / [pie_j/(1-pie_j)] and is recommended to be at least 5.

SAMPLE DATA OUTPUT DATASETS (prefix and k as supplied):
ordprob_single_POST_K3    BENE_ID   post1   post2   post3   maxpost   classnum   class
ordprob_single_CLASS_K3   class     n       percent   pie     app       occ
##########################################################################################################################*/
%MACRO ordprob_mix_post_one(
  data=BASE_FILE_SRS, id=BENE_ID,
  yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
  m=4, ycodes=0 1 2 3, deg=2, k=4,
  outestlib=work, prefix=ordprob_single
);
%LOCAL m1 d j;
%LET m1=%EVAL(&m-1);
%LET d=&deg;

*Hard stop if either input is missing;
%IF NOT %SYSFUNC(EXIST(&data)) %THEN %DO;
%PUT ERROR: DATA=&data not found.;
%RETURN;
%END;
%IF NOT %SYSFUNC(EXIST(&outestlib..&prefix._EST_K&k)) %THEN %DO;
%PUT ERROR: &outestlib..&prefix._EST_K&k not found. Run ORDPROB_MIX_FIT_ONE first with the same prefix and k.;
%RETURN;
%END;

*Category codes to macro variables;
%DO j=1 %TO &m; %GLOBAL ycode&j; %LET ycode&j=%SCAN(&ycodes,&j,%STR( )); %END;

*--- Parameter estimates to one wide row -------------------------------------;
DATA _pest;
SET &outestlib..&prefix._EST_K&k;
KEEP Parameter Estimate;
RUN;

PROC TRANSPOSE DATA=_pest OUT=_pwide(DROP=_NAME_);
ID Parameter;
VAR Estimate;
RUN;

*Class A has no alpha0_ in ParameterEstimates because it is fixed at 0, so create it and set it to 0;
DATA _pwide;
SET _pwide;
%DO c=1 %TO &k;
%LET U=%CL(&c);
IF 0 THEN alpha0_&U=0;
IF MISSING(alpha0_&U) THEN alpha0_&U=0;
%END;
RUN;

*--- Estimated mixing proportions, one row per class -------------------------;
DATA _pi;
SET _pwide;
LENGTH class $1;
denom=0;
%DO c=1 %TO &k; %LET U=%CL(&c); denom = denom + EXP(alpha0_&U); %END;
%DO c=1 %TO &k; %LET U=%CL(&c);
class="&U"; pie=EXP(alpha0_&U)/denom; OUTPUT;
%END;
KEEP class pie;
RUN;

*--- Class log-likelihoods, posteriors, modal assignment ---------------------;
DATA &outestlib..&prefix._POST_K&k;
IF _N_=1 THEN SET _pwide;
SET &data;
ARRAY yv[&ttotal] &yvars;
ARRAY tv[&ttotal] &tvars;
%DO c=1 %TO &k; %LET U=%CL(&c);
ARRAY z&U[&m] %DO j=1 %TO &m; z&j._&U %END;;
%END;
ARRAY pst[&k] %DO c=1 %TO &k; post&c %END;;

*Class log-likelihoods;
%DO c=1 %TO &k; llik&c=0; %END;

DO t=1 TO &ttotal;
tval = tv[t];
cat = .; %DO j=1 %TO &m; IF yv[t] = &&ycode&j THEN cat=&j; %END;

*cat stays missing when the value is blank or is not one of the ycodes, so that time point is skipped;
IF cat>0 THEN DO;
%DO c=1 %TO &k;
%LET U=%CL(&c); %LET L=%LOWCASE(&U);

eta_&U = beta_&L.0;
%IF &d>=1 %THEN %DO; eta_&U = eta_&U + beta_&L.1*(tval); %END;
%IF &d>=2 %THEN %DO j=2 %TO &d; eta_&U = eta_&U + beta_&L.&j.*((tval)**&j); %END;

th1_&U = 0;
%IF &m1>=2 %THEN %DO j=2 %TO &m1;
%IF &j=2 %THEN %DO; th&j._&U = th1_&U + i&L.&j; %END;
%ELSE %DO;          th&j._&U = th%EVAL(&j-1)_&U + i&L.&j; %END;
%END;

%DO j=1 %TO &m;
%IF &j=1 %THEN %DO;
z&j._&U = PROBNORM(th1_&U - eta_&U);
%END; %ELSE %IF &j<&m %THEN %DO;
z&j._&U = PROBNORM(th&j._&U - eta_&U) - PROBNORM(th%EVAL(&j-1)_&U - eta_&U);
%END; %ELSE %DO;
z&j._&U = 1 - PROBNORM(th%EVAL(&m1)_&U - eta_&U);
%END;
%END;

llik&c = llik&c + LOG(MAX(1e-12, z&U[cat]));
%END;
END;
END;

*Softmax mixing weights, matching the fit;
denom=0;
%DO c=1 %TO &k; %LET U=%CL(&c); pin&c = EXP(alpha0_&U); denom + pin&c; %END;
%DO c=1 %TO &k; pie&c = pin&c/denom; %END;

*Posterior class membership;
mix=0; %DO c=1 %TO &k; mix + pie&c*EXP(llik&c); %END;
mix=MAX(mix,1e-300);
%DO c=1 %TO &k; post&c = pie&c*EXP(llik&c)/mix; %END;

*Modal assignment;
maxpost  = LARGEST(1, OF pst[*]);
classnum = WHICHN(maxpost, OF pst[*]);
LENGTH class $1;
class = BYTE(64+classnum);

KEEP &id %DO c=1 %TO &k; post&c %END; maxpost classnum class;
RUN;

*--- Class sizes, APP and OCC ------------------------------------------------;
PROC SQL;
CREATE TABLE _csize AS
SELECT class,
       COUNT(*) AS n,
       100*COUNT(*)/(SELECT COUNT(*) FROM &outestlib..&prefix._POST_K&k) AS percent FORMAT=6.2,
       MEAN(maxpost) AS app FORMAT=6.3
FROM &outestlib..&prefix._POST_K&k
GROUP BY class;
QUIT;

PROC SORT DATA=_csize; BY class; RUN;
PROC SORT DATA=_pi;    BY class; RUN;

DATA &outestlib..&prefix._CLASS_K&k;
MERGE _csize _pi;
BY class;
FORMAT pie 6.3 occ 8.1;
IF 0<app<1 AND 0<pie<1 THEN occ = (app/(1-app)) / (pie/(1-pie));
ELSE occ = .;
RUN;

*--- QC ----------------------------------------------------------------------;
PROC PRINT DATA=&outestlib..&prefix._CLASS_K&k NOOBS;
VAR class n percent pie app occ;
TITLE "Class sizes, average posterior probability and odds of correct classification";
TITLE2 "&k. classes. APP should be at least 0.7 and OCC at least 5";
RUN;
TITLE;
TITLE2;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE _pest _pwide _pi _csize;
QUIT;
%MEND ordprob_mix_post_one;

/*##########################################################################################################################
*STEP 6: PLOTS
##########################################################################################################################
Class proportions recovered from the alpha0_* estimates, and predicted mean outcome by class and time.
Class A is absent from ParameterEstimates because its mixing intercept is fixed, so it is added back with
alpha=0 before the softmax. The first threshold is set to 0 here to match the fit.

SAMPLE DATA OUTPUT DATASETS:
mix_props    class   alpha    expa     pie
             A       0        1        0.344
pred_long    class   tval     EY       p1   p2   p3   p4
             A       1        1.4712   0.11 0.35 0.42 0.12
##########################################################################################################################*/
%MACRO ordprob_mix_plot_one(
  outestlib=work, prefix=ordprob_single,
  k=4, ttotal=12, m=4, ycodes=0 1 2 3
);
%LOCAL m1 j;
%LET m1=%EVAL(&m-1);
%DO j=1 %TO &m; %GLOBAL ycode&j; %LET ycode&j=%SCAN(&ycodes,&j,%STR( )); %END;

*Mixing proportions from alpha0_*. Add class A back because it is fixed at 0 and not returned;
DATA _mix;
SET &outestlib..&prefix._EST_K&k(KEEP=Parameter Estimate);
LENGTH class $1 alpha 8;
IF UPCASE(SUBSTR(Parameter,1,7))='ALPHA0_' THEN DO;
class = SUBSTR(Parameter,8,1);
alpha = Estimate;
OUTPUT;
END;
RUN;

PROC SQL NOPRINT; SELECT COUNT(*) INTO :_hasA FROM _mix WHERE UPCASE(class)='A'; QUIT;
%IF %SYSEVALF(&_hasA=0) %THEN %DO;
DATA _addA; LENGTH class $1 alpha 8; class='A'; alpha=0; RUN;
PROC APPEND BASE=_mix DATA=_addA FORCE; RUN;
PROC DATASETS LIB=WORK NOLIST; DELETE _addA; QUIT;
%END;

DATA _mix; SET _mix; class=UPCASE(class); RUN;
PROC SORT DATA=_mix; BY class; RUN;

DATA mix_props;
SET _mix END=last;
RETAIN sumexp 0;
expa=EXP(alpha);
sumexp + expa;
IF last THEN CALL SYMPUTX('_sumexp',sumexp);
RUN;

DATA mix_props; SET mix_props; pie = expa / &_sumexp; RUN;

PROC SGPLOT DATA=mix_props;
FORMAT pie 5.3;
VBAR class / RESPONSE=pie DATALABEL BARWIDTH=0.6;
YAXIS GRID LABEL='Mixing proportion' MIN=0;
XAXIS LABEL='Class';
RUN;

*Predicted trajectories;
PROC SQL;
CREATE TABLE _betas AS
SELECT UPCASE(SUBSTR(Parameter,6,1)) AS CL LENGTH=1,
INPUT(SUBSTR(Parameter,7), BEST.) AS IDX,
Estimate
FROM &outestlib..&prefix._EST_K&k
WHERE UPCASE(SUBSTR(Parameter,1,5))='BETA_';
QUIT;
PROC SORT DATA=_betas; BY CL IDX; RUN;
PROC TRANSPOSE DATA=_betas OUT=betas_w PREFIX=beta_; BY CL; ID IDX; VAR Estimate; RUN;

*Thresholds. ESTS holds threshold2 to threshold(m-1). threshold1 is fixed at 0;
PROC SQL;
CREATE TABLE _ths AS
SELECT UPCASE(SCAN(Label,2,'_')) AS CL LENGTH=1,
INPUT(COMPRESS(SUBSTR(LOWCASE(Label),10),,'kd'), BEST.) AS IDX,
Estimate
FROM &outestlib..&prefix._ESTS_K&k
WHERE UPCASE(SUBSTR(Label,1,9))='THRESHOLD';
QUIT;
PROC SORT DATA=_ths; BY CL IDX; RUN;
PROC TRANSPOSE DATA=_ths OUT=ths_w PREFIX=th_; BY CL; ID IDX; VAR Estimate; RUN;

PROC SORT DATA=betas_w; BY CL; RUN;
PROC SORT DATA=ths_w; BY CL; RUN;
DATA coeffs; MERGE betas_w ths_w; BY CL; RUN;

DATA pred_long;
SET coeffs;
LENGTH class $1;
class=CL;
ARRAY b  beta_0-beta_99;*only the existing indices are used;
ARRAY th th_1-th_99;
ARRAY p  p1-p&m;

th_1 = 0;*first threshold fixed at 0, matching the fit;

DO tval=1 TO &ttotal;
eta=0;
DO j=1 TO DIM(b);
IF NOT MISSING(b[j]) THEN eta + b[j]*(tval**(j-1));
END;

p[1] = PROBNORM(th[1] - eta);
%IF &m1>=2 %THEN %DO jj=2 %TO &m1;
p[&jj] = PROBNORM(th[&jj] - eta) - PROBNORM(th[%EVAL(&jj-1)] - eta);
%END;
p[&m] = 1 - PROBNORM(th[&m1] - eta);

EY=0; %DO j=1 %TO &m; EY + %SCAN(&ycodes,&j,%STR( ))*p[&j]; %END;
OUTPUT;
END;
KEEP class tval EY p1-p&m;
RUN;

PROC SGPLOT DATA=pred_long;
SERIES X=tval Y=EY / GROUP=class MARKERS
LINEATTRS=(THICKNESS=2) MARKERATTRS=(SIZE=8);
XAXIS VALUES=(1 TO &ttotal BY 1) LABEL='Quarter (t)';
YAXIS LABEL='E[Y_t | class]';
KEYLEGEND / TITLE='Class';
RUN;
%MEND ordprob_mix_plot_one;

/*##########################################################################################################################
*END
##########################################################################################################################
This file defines macros only. Including it runs nothing.

TYPICAL WORKFLOW
  %include "ordinal_probit.sas";

  %sim_data(dist=BIN4, class=3, n=500, T=12, seed=2026, cap=3);
  %build_base_from_simwide(T=12);

  %ordprob_mix_fit_one(
    data=BASE_FILE_SRS, id=BENE_ID,
    yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
    m=4, ycodes=0 1 2 3, deg=2, k=3,
    prefix=ordprob_T1C3
  );

  %ordprob_mix_post_one(
    data=BASE_FILE_SRS, id=BENE_ID,
    yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
    m=4, ycodes=0 1 2 3, deg=2, k=3,
    outestlib=work, prefix=ordprob_T1C3
  );

  %ordprob_mix_plot_one(
    outestlib=work, prefix=ordprob_T1C3,
    k=3, ttotal=12, m=4, ycodes=0 1 2 3
  );

TO USE YOUR OWN DATA INSTEAD OF THE SIMULATOR
  yvars = your wide-format outcome columns, for example Y1_1-Y1_12
  tvars = quar1-quar12 holding the numeric time values 1 to T. Create them if needed
  id    = your subject identifier, or rename it to BENE_ID

MISSING AND OUT-OF-RANGE VALUES
  A time point whose outcome value is blank, or is not one of the values listed in ycodes, contributes
  nothing to the likelihood. The subject still contributes every time point that does match. Values that
  are present but outside ycodes are dropped silently, so check the coding of the outcome before fitting.

IDENTIFICATION
  The first threshold is fixed at 0 and the class intercept beta_*0 is free. The free threshold
  parameters are the INCREMENTS i*2, i*3 and onward. There is no i*1. Each increment is held positive by
  a BOUNDS statement that the macro generates, so the thresholds stay strictly ordered within a class.
  If you supply start values for thresholds, start from i*2 and use positive values, for example

    start_values=%str(alpha0_B=-0.5 beta_a0=0.1 ib2=1 ib3=1)

  A start value of 0 for an increment sits on the boundary and should not be used. When start_values= is
  left blank the increments default to 1, which places the thresholds at 0, 1, 2 and so on.

  Never supply alpha0_A. It is fixed at 0 for identification.

OUTPUT DATASETS
  SIM_WIDE                  wide simulated outcomes, one row per subject
  BASE_FILE_SRS             the fitting-macro input contract
  <prefix>_EST_K<k>         parameter estimates
  <prefix>_FIT_K<k>         fit statistics including BIC
  <prefix>_ESTS_K<k>        thresholds on the natural scale
  <prefix>_POST_K<k>        posterior class-membership probabilities and modal class, one row per subject
  <prefix>_CLASS_K<k>       class sizes, estimated proportions, APP and OCC
  mix_props                 estimated class proportions
  pred_long                 predicted mean outcome by class and time
##########################################################################################################################*/
