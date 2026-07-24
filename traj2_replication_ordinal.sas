*************************************************************************************************************************
*************************************************************************************************************************
COMMUNITY HEALTH AND AGING OUTCOMES (CHAO) LAB - INSTITUTE FOR HEALTH, HEALTH CARE POLICY & AGING RESEARCH - 
RUTGERS, THE STATE UNIVERSITY OF NEW JERSEY                  
*************************************************************************************************************************
*************************************************************************************************************************
*PROJECT NAME: Traj2 Replication, Ordinal-Probit Outcome Family
LAST UPDATED DATE: 23 JUL 2026
DATA SOURCES: NONE. All input is simulated in STEP 2 from the data generating process declared in STEP 0
PURPOSE: Single standalone script reproducing the ordinal-outcome tables and figures in
Traj2: A Native Macro Library for Single and Multi-Outcome Group-Based Trajectory Modeling in SAS
Zafar, Xia, Lin, Jarrin, Journal of Statistical Software. Produces
(a)STEP 3, Table 14. Ordinal DGP recovery, 200 replications, N=2,500, SEED_MC=20260000
(b)STEP 4, Table 15 and Figure 2. Cross-implementation check against the R package lcmm, SEED_CC=20260601
(c)STEP 5, Figures 5 and 6. Worked example of paper Section 5.2, N=500, seed 2026
COMPANION: traj2_lcmm_crosscheck.R supplies the lcmm half of STEP 4. Run this script first so the two CSV
files exist, then run the R script, then read its output into the Table 15 shell printed at the end of STEP 4.
NOT REPRODUCIBLE: Figure 1 and paper Section 3.6, the 2019 Medicare decedent cohort of 181,463 subjects.
Those claims data live inside the CMS VRDC under a data use agreement and cannot be distributed.
AUTHOR: Anum Zafar, Weiyi Xia, Haiqun Lin, Olga F. Jarrin
##########################################################################################################################
*Execution Environment: SAS 9.4 or later with SAS/STAT PROC NLMIXED. No compiled components           *
*ordinal_probit.sas must sit in SRCPATH. OUTPATH is built automatically if its parent folder exists    *
*Runtime: STEP 3 is the slow one. It fits 6 models per replication, so 200 replications is 1200 fits   *
*Smoke test at NREP=5 before committing to the full run                                                *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 0:Parameters. Edit only this block
#STEP 1:Load the macro library
#STEP 2:Shared ordinal simulator and the shape-recovery truth
#STEP 3:Table 14. Ordinal-probit DGP recovery
#STEP 4:Table 15 and Figure 2. Cross-implementation check against lcmm
#STEP 5:Figures 5 and 6. Worked example
#STEP 6:Check output
#STEP 7:Clean up the WORK library
##########################################################################################################################;


/*##########################################################################################################################
*STEP 0: PARAMETERS (EDIT THESE ONLY)
##########################################################################################################################*/
%LET SRCPATH=T:\datasets_t;                *folder holding ordinal_probit.sas;
%LET OUTPATH=T:\datasets_t\replication;    *folder for the PDF and CSV output. Must already exist;

%LET RUN_RECOVERY=1;   *1=run STEP 3, Table 14;
%LET RUN_CROSSCHECK=1; *1=run STEP 4, the SAS half of Table 15 and Figure 2;
%LET RUN_EXAMPLE=1;    *1=run STEP 5, Figures 5 and 6;

%LET NREP=200;         *replications for Table 14. Use 5 for a smoke test first;
%LET NSIM=2500;        *subjects per replication;
%LET SEED_MC=20260000; *replication r uses seed SEED_MC+r, so 20260001 through 20260200;
%LET SEED_CC=20260601; *single seed for the lcmm cross-check;

%LET KTRUE=3;          *true number of latent classes;
%LET KMIN=2;           *BIC sweep low end;
%LET KMAX=6;           *BIC sweep high end;
%LET T=12;             *time points;
%LET MCAT=4;           *ordinal categories, coded 0 1 2 3;
%LET DEG=2;            *polynomial degree, quadratic;
%LET IDVAR=BENE_ID;    *subject identifier;

*Cumulative probit DGP, described in paper Section 4.4. Class-specific
 quadratic mean trajectories on the latent scale, class-common thresholds,
 no subject-level random effect, and non-uniform class proportions so the
 run also tests recovery of the mixing weights;
%LET A_B0=-1.0;  %LET A_B1=0.00;  %LET A_B2=0.000;   *low and flat;
%LET B_B0=-0.8;  %LET B_B1=0.18;  %LET B_B2=0.000;   *moderate rise;
%LET C_B0=-0.3;  %LET C_B1=0.10;  %LET C_B2=0.012;   *starts mid, accelerates;

*class-common threshold increments: tau1=I1, tau2=tau1+EXP(I2), tau3=tau2+EXP(I3);
%LET I1=0.4;  %LET I2=-0.1;  %LET I3=0.0;

*true class proportions, must sum to 1;
%LET PA=0.45;  %LET PB=0.35;  %LET PC=0.20;

OPTIONS NOFMTERR MPRINT NOMLOGIC NOSYMBOLGEN;

*Create OUTPATH if it is not there yet. Without this an ODS PDF FILE=
 opens against a missing folder and only fails later, at ODS PDF CLOSE,
 with "Physical file does not exist". The parent folder must already exist;
OPTIONS DLCREATEDIR;
LIBNAME _mkdir "&OUTPATH.";
LIBNAME _mkdir CLEAR;
OPTIONS NODLCREATEDIR;


/*##########################################################################################################################
*STEP 1: LOAD THE MACRO LIBRARY

  ordinal_probit.sas defines macros only and runs nothing on
  %include, so no demo switch is needed here.

  The released version fixes the first threshold at 0 and estimates
  the class intercept beta_*0 freely. The free threshold parameters
  are therefore the increments i*2, i*3 and so on. There is no i*1
  and there is no alpha0_A, which is fixed at 0 for identification.
##########################################################################################################################*/
%INCLUDE "&SRCPATH.\ordinal_probit.sas";


/*##########################################################################################################################
*STEP 2: SHARED ORDINAL SIMULATOR AND TRUTH

  %sim_ord_dgp generates one replication of the Section 4.4 DGP in
  the wide contract of Table 1. It draws one uniform for the class
  and then one uniform per time point, so a given seed always yields
  the same data set. STEP 3 and STEP 4 differ only in the seed.

  true_traj holds E[Y_t | class] on the 0-3 outcome scale under the
  generating parameters. Shape recovery is measured against this.
##########################################################################################################################*/
%MACRO sim_ord_dgp(out=BASE_FILE_SRS, seed=, n=&NSIM., tpts=&T., id=&IDVAR.);
DATA &out.;
CALL STREAMINIT(&seed.);
ARRAY B0[3] _TEMPORARY_ (&A_B0. &B_B0. &C_B0.);
ARRAY B1[3] _TEMPORARY_ (&A_B1. &B_B1. &C_B1.);
ARRAY B2[3] _TEMPORARY_ (&A_B2. &B_B2. &C_B2.);
ARRAY Y[&tpts.] Y1_1-Y1_&tpts.;
ARRAY quar[&tpts.] quar1-quar&tpts.;
tau1=&I1.;
tau2=tau1+EXP(&I2.);
tau3=tau2+EXP(&I3.);
cpA=&PA.;
cpB=&PA.+&PB.;
DO &id.=1 TO &n.;
u=RAND('UNIFORM');
IF u<cpA THEN trueclass=1;
ELSE IF u<cpB THEN trueclass=2;
ELSE trueclass=3;
DO t=1 TO &tpts.;
quar[t]=t;
eta=B0[trueclass]+B1[trueclass]*t+B2[trueclass]*t*t;
p0=PROBNORM(tau1-eta);
p1=PROBNORM(tau2-eta)-PROBNORM(tau1-eta);
p2=PROBNORM(tau3-eta)-PROBNORM(tau2-eta);
p3=1-PROBNORM(tau3-eta);
uu=RAND('UNIFORM');
IF uu<p0 THEN Y[t]=0;
ELSE IF uu<p0+p1 THEN Y[t]=1;
ELSE IF uu<p0+p1+p2 THEN Y[t]=2;
ELSE Y[t]=3;
END;
OUTPUT;
END;
KEEP &id. trueclass Y1_1-Y1_&tpts. quar1-quar&tpts.;
RUN;
%MEND sim_ord_dgp;

/*##########################################################################################################################
  SAMPLE DATA OUTPUT DATASET: BASE_FILE_SRS (from %sim_ord_dgp)

  BENE_ID trueclass Y1_1 Y1_2 Y1_3 ... Y1_12 quar1 quar2 ... quar12
        1         1    3    0    1 ...     0     1     2 ...     12
        2         3    1    2    2 ...     3     1     2 ...     12

  One row per subject. Y1_t takes the values 0 1 2 3. trueclass is
  carried for validation only and no fitting macro reads it.
##########################################################################################################################*/

DATA true_traj;
ARRAY B0[3] _TEMPORARY_ (&A_B0. &B_B0. &C_B0.);
ARRAY B1[3] _TEMPORARY_ (&A_B1. &B_B1. &C_B1.);
ARRAY B2[3] _TEMPORARY_ (&A_B2. &B_B2. &C_B2.);
t1=&I1.;
t2=t1+EXP(&I2.);
t3=t2+EXP(&I3.);
DO trueclass=1 TO 3;
DO t=1 TO &T.;
eta=B0[trueclass]+B1[trueclass]*t+B2[trueclass]*t*t;
p0=PROBNORM(t1-eta);
p1=PROBNORM(t2-eta)-PROBNORM(t1-eta);
p2=PROBNORM(t3-eta)-PROBNORM(t2-eta);
p3=1-PROBNORM(t3-eta);
EY_true=0*p0+1*p1+2*p2+3*p3;
OUTPUT;
END;
END;
KEEP trueclass t EY_true;
RUN;

*%ord_coef pulls the betas, thresholds and mixing proportions of one fit out of the
 three ODS data sets the fitting macro writes, and lands them one row per class.
 Used by STEP 3 and STEP 4, which is why it lives here;
%MACRO ord_coef(prefix=, k=, out=_coef);
PROC SQL NOPRINT;
CREATE TABLE _b AS
SELECT UPCASE(SUBSTR(Parameter,6,1)) AS CL LENGTH=1,
       INPUT(SUBSTR(Parameter,7),BEST.) AS IDX, Estimate
FROM work.&prefix._EST_K&k.
WHERE UPCASE(SUBSTR(Parameter,1,5))='BETA_';
QUIT;
PROC TRANSPOSE DATA=_b OUT=_bw PREFIX=b_; BY CL; ID IDX; VAR Estimate; RUN;

*ESTS holds threshold2 through threshold(m-1). Threshold 1 is fixed at 0
 in the released macro and is supplied below when _coef is built;
PROC SQL NOPRINT;
CREATE TABLE _th AS
SELECT UPCASE(SCAN(Label,2,'_')) AS CL LENGTH=1,
       INPUT(COMPRESS(SUBSTR(LOWCASE(Label),10),,'kd'),BEST.) AS IDX, Estimate
FROM work.&prefix._ESTS_K&k.
WHERE UPCASE(SUBSTR(Label,1,9))='THRESHOLD';
QUIT;
PROC TRANSPOSE DATA=_th OUT=_thw PREFIX=th_; BY CL; ID IDX; VAR Estimate; RUN;

*class proportions from the alpha0_* estimates of this fit. Class A is fixed
 at 0 and never appears in ParameterEstimates, so it is added by hand;
DATA _alpha;
LENGTH CL $1 alpha 8;
CL='A';
alpha=0;
RUN;

DATA _alphaB;
SET work.&prefix._EST_K&k.(KEEP=Parameter Estimate);
LENGTH CL $1 alpha 8;
IF UPCASE(SUBSTR(Parameter,1,7))='ALPHA0_';
CL=UPCASE(SUBSTR(Parameter,8,1));
alpha=Estimate;
KEEP CL alpha;
RUN;

PROC APPEND BASE=_alpha DATA=_alphaB FORCE; RUN;

PROC SQL NOPRINT;
SELECT SUM(EXP(alpha)) INTO :_sumexp TRIMMED FROM _alpha;
QUIT;

DATA _mixp;
SET _alpha;
pie=EXP(alpha)/&_sumexp.;
KEEP CL pie;
RUN;

PROC SORT DATA=_mixp; BY CL; RUN;
PROC SORT DATA=_bw; BY CL; RUN;
PROC SORT DATA=_thw; BY CL; RUN;

DATA &out.;
MERGE _bw _thw _mixp;
BY CL;
th_1=0;   *first threshold fixed at 0 in the released macro;
classnum=(CL='A')*1+(CL='B')*2+(CL='C')*3;
RUN;
%MEND ord_coef;


/*##########################################################################################################################
*STEP 3: TABLE 14
  Ordinal probit DGP recovery over NREP replications.

  Three quantities are recorded per replication:

    Shape recovery      worst-class RMSE between fitted and true
                        E[Y_t | class] on the 0-3 scale, after
                        permuting fitted classes to minimise the
                        overall trajectory distance. The worst class
                        is reported rather than the average, so the
                        number is not carried by the classes that are
                        easiest to recover.

    Classification      proportion of subjects whose maximum
                        posterior class matches the simulated class.
                        The majority-class baseline is 45 percent.

    Model selection     the C that minimises BIC over C=KMIN..KMAX.

  A replication counts only if the estimated thresholds are correctly
  ordered, the betas are finite, and the fitted-to-true class map is
  one to one. The clean-Hessian column records whether every free
  parameter came back with a usable standard error.
##########################################################################################################################*/
%MACRO run_recovery;
%IF &RUN_RECOVERY. NE 1 %THEN %RETURN;

%LOCAL r kk seed;
%GLOBAL nrep_seen ndrop nclean nse;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE mc_metrics mc_bicpick;
QUIT;

ODS GRAPHICS OFF;
ODS EXCLUDE ALL;
ODS NORESULTS;
OPTIONS NONOTES NOSOURCE NOSOURCE2 NOMPRINT;

%DO r=1 %TO &NREP.;
%LET seed=%EVAL(&SEED_MC. + &r.);
%sim_ord_dgp(out=BASE_FILE_SRS, seed=&seed.);

*fit at the true number of classes;
%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=&IDVAR.,
  yvars=Y1_1-Y1_&T., tvars=quar1-quar&T., ttotal=&T.,
  m=&MCAT., ycodes=0 1 2 3, deg=&DEG., k=&KTRUE.,
  tech=quanew, maxiter=2000, prefix=oa
);

%IF %SYSFUNC(EXIST(work.oa_EST_K&KTRUE.)) %THEN %DO;

%ord_coef(prefix=oa, k=&KTRUE., out=_coef);

*clean Hessian check: are all standard errors present and positive;
%LET _hasSE=0;
PROC SQL NOPRINT;
SELECT COUNT(*) INTO :_hasSE TRIMMED
FROM dictionary.columns
WHERE libname='WORK' AND UPCASE(memname)=UPCASE("oa_EST_K&KTRUE.")
  AND UPCASE(name)='STANDARDERROR';
QUIT;
%LET nbadse=.;
%IF &_hasSE. %THEN %DO;
PROC SQL NOPRINT;
SELECT SUM(MISSING(StandardError) OR StandardError<=0) INTO :nbadse TRIMMED
FROM work.oa_EST_K&KTRUE.;
QUIT;
%END;

*(1) shape recovery;
DATA _fit_traj;
SET _coef;
ARRAY b b_0 b_1 b_2;
DO t=1 TO &T.;
eta=0;
DO j=1 TO DIM(b);
IF NOT MISSING(b[j]) THEN eta+b[j]*(t**(j-1));
END;
p0=PROBNORM(th_1-eta);
p1=PROBNORM(th_2-eta)-PROBNORM(th_1-eta);
p2=PROBNORM(th_3-eta)-PROBNORM(th_2-eta);
p3=1-PROBNORM(th_3-eta);
EY_fit=0*p0+1*p1+2*p2+3*p3;
OUTPUT;
END;
KEEP classnum t EY_fit;
RUN;

PROC SQL;
CREATE TABLE _pair AS
SELECT f.classnum AS fclass, tr.trueclass, SUM((f.EY_fit-tr.EY_true)**2) AS sse
FROM _fit_traj f CROSS JOIN true_traj tr
WHERE f.t=tr.t
GROUP BY f.classnum, tr.trueclass;
CREATE TABLE _match AS
SELECT fclass, trueclass AS matched_true, sse
FROM _pair p
WHERE sse=(SELECT MIN(sse) FROM _pair p2 WHERE p2.fclass=p.fclass);
QUIT;

PROC SORT DATA=_match NODUPKEY; BY fclass; RUN;
PROC SQL NOPRINT;
SELECT COUNT(DISTINCT matched_true) INTO :ndist TRIMMED FROM _match;
QUIT;

PROC SQL;
CREATE TABLE _shape AS
SELECT &r. AS rep, (&ndist.=&KTRUE.) AS bijective, MAX(SQRT(sse/&T.)) AS max_rmse
FROM _match;
QUIT;

*(2) classification by maximum posterior membership, prior weighted argmax;
DATA _assign;
IF _N_=1 THEN DO;
IF 0 THEN SET _coef;
DECLARE HASH h(dataset:'_coef');
_rc=h.defineKey('classnum');
_rc=h.defineData('b_0','b_1','b_2','th_1','th_2','th_3','pie');
_rc=h.defineDone();
END;
SET BASE_FILE_SRS;
ARRAY Y Y1_1-Y1_&T.;
ARRAY bb[3,0:2] _TEMPORARY_;
ARRAY tt[3,3] _TEMPORARY_;
ARRAY pw[3] _TEMPORARY_;
DO cc=1 TO 3;
CALL MISSING(b_0,b_1,b_2,th_1,th_2,th_3,pie);
_r=h.find(key:cc);
bb[cc,0]=COALESCE(b_0,0);
bb[cc,1]=COALESCE(b_1,0);
bb[cc,2]=COALESCE(b_2,0);
tt[cc,1]=th_1;
tt[cc,2]=th_2;
tt[cc,3]=th_3;
pw[cc]=pie;
END;
bestll=.;
assigned=.;
DO cc=1 TO 3;
ll=0;
DO t=1 TO &T.;
eta=bb[cc,0]+bb[cc,1]*t+bb[cc,2]*t*t;
IF Y[t]=0 THEN pr=PROBNORM(tt[cc,1]-eta);
ELSE IF Y[t]=1 THEN pr=PROBNORM(tt[cc,2]-eta)-PROBNORM(tt[cc,1]-eta);
ELSE IF Y[t]=2 THEN pr=PROBNORM(tt[cc,3]-eta)-PROBNORM(tt[cc,2]-eta);
ELSE pr=1-PROBNORM(tt[cc,3]-eta);
ll + LOG(MAX(pr,1e-12));
END;
IF NOT MISSING(pw[cc]) THEN ll + LOG(MAX(pw[cc],1e-12));
IF bestll=. OR ll>bestll THEN DO;
bestll=ll;
assigned=cc;
END;
END;
KEEP &IDVAR. trueclass assigned;
RUN;

PROC SQL;
CREATE TABLE _asg2 AS
SELECT a.*, m.matched_true AS assigned_true
FROM _assign a LEFT JOIN _match m ON a.assigned=m.fclass;
CREATE TABLE _acc AS
SELECT &r. AS rep, MEAN(assigned_true=trueclass) AS agree FROM _asg2;
QUIT;

PROC SQL NOPRINT;
SELECT SUM( (th_1<th_2 AND th_2<th_3)
            AND (b_0 IS NOT NULL AND b_1 IS NOT NULL AND b_2 IS NOT NULL) )
INTO :nok TRIMMED FROM _coef;
QUIT;

DATA _row;
MERGE _shape _acc;
BY rep;
valid=(&nok.=&KTRUE.) AND bijective;
%IF &_hasSE. %THEN %DO;
clean_hess=(&nbadse.=0);
%END;
%ELSE %DO;
clean_hess=.;
%END;
RUN;

PROC APPEND BASE=mc_metrics DATA=_row FORCE; RUN;
%END;

*(3) BIC model selection over C=KMIN..KMAX;
%DO kk=&KMIN. %TO &KMAX.;
%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=&IDVAR.,
  yvars=Y1_1-Y1_&T., tvars=quar1-quar&T., ttotal=&T.,
  m=&MCAT., ycodes=0 1 2 3, deg=&DEG., k=&kk.,
  tech=quanew, maxiter=2000, prefix=oak&kk.
);
%IF %SYSFUNC(EXIST(work.oak&kk._FIT_K&kk.)) %THEN %DO;
DATA _bicrow;
SET work.oak&kk._FIT_K&kk.;
WHERE UPCASE(Descr) LIKE '%BIC%';
rep=&r.;
kfit=&kk.;
bicval=Value;
KEEP rep kfit bicval;
RUN;
PROC APPEND BASE=mc_bicpick DATA=_bicrow FORCE; RUN;
%END;
%END;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE oa_EST_K&KTRUE. oa_ESTS_K&KTRUE. oa_FIT_K&KTRUE.
       _b _bw _th _thw _alpha _alphaB _mixp _coef _fit_traj _pair _match _shape
       _assign _asg2 _acc _row _bicrow
  %DO kk=&KMIN. %TO &KMAX.; oak&kk._EST_K&kk. oak&kk._ESTS_K&kk. oak&kk._FIT_K&kk. %END; ;
QUIT;

%IF %SYSFUNC(MOD(&r.,25))=0 %THEN %PUT NOTE: ===== completed replication &r. of &NREP. =====;
%END;

OPTIONS NOTES SOURCE SOURCE2 MPRINT;
ODS EXCLUDE NONE;
ODS RESULTS;
ODS GRAPHICS ON;

*summarise;
PROC SORT DATA=mc_bicpick; BY rep bicval; RUN;
DATA _bestk;
SET mc_bicpick;
BY rep;
IF first.rep;
RUN;

PROC SQL NOPRINT;
SELECT COUNT(DISTINCT rep) INTO :nrep_seen TRIMMED FROM mc_bicpick;
SELECT SUM(valid=0) INTO :ndrop TRIMMED FROM mc_metrics;
SELECT SUM(clean_hess=1), COUNT(clean_hess) INTO :nclean TRIMMED, :nse TRIMMED FROM mc_metrics;
QUIT;

ODS PDF FILE="&OUTPATH.\traj2_ordinal_validation.pdf" STARTPAGE=NO;
TITLE "Table 14. Traj2 ordinal probit DGP recovery (&NREP. reps, N=&NSIM., T=&T., true C=&KTRUE.)";
TITLE2 'worst-class RMSE is on the 0-3 outcome scale; agreement baseline is the 45 percent majority class';
PROC MEANS DATA=mc_metrics(WHERE=(valid=1)) N MEAN STD MIN MAX MAXDEC=4;
VAR max_rmse agree;
RUN;

TITLE "Model selection: C chosen by minimum BIC over C=&KMIN. to &KMAX. (target C=&KTRUE.)";
TITLE2 "replications with a BIC sweep: &nrep_seen.";
PROC FREQ DATA=_bestk;
TABLES kfit / NOCUM;
RUN;

TITLE 'Estimability: replications with fully estimable standard errors';
TITLE2 "clean-Hessian fits: &nclean. of &nse.; replications dropped on the validity filter: &ndrop.";
PROC FREQ DATA=mc_metrics;
TABLES valid*clean_hess / LIST MISSING;
RUN;
TITLE;
TITLE2;
ODS PDF CLOSE;

PROC EXPORT DATA=mc_metrics
  OUTFILE="&OUTPATH.\ordinal_recovery_metrics.csv"
  DBMS=CSV REPLACE;
RUN;

%PUT NOTE: Dropped &ndrop. rep(s) on the validity and bijective filter. Clean-Hessian fits: &nclean. of &nse.;
%MEND run_recovery;
%run_recovery;


/*##########################################################################################################################
*STEP 4: TABLE 15 AND FIGURE 2
  Cross-implementation check against the R package lcmm, SAS half.

  One data set from the same DGP as STEP 3 at a single fixed seed.
  It is fitted here with %ordprob_mix_fit_one and exported twice:

    ordinal_data_for_R.csv   the raw long-format data, for lcmm
    traj2_EY_by_class.csv    Traj2 fitted E[Y_t | class] on the 0-3 scale

  The two programs are compared on E[Y | class, t], class proportions
  and the log likelihood, all of which are invariant to how each
  program parameterises the model. Raw thresholds and betas are NOT
  compared, because the two use different conventions.

  Run traj2_lcmm_crosscheck.R next, then fill the Table 15 shell.
##########################################################################################################################*/
%MACRO run_crosscheck;
%IF &RUN_CROSSCHECK. NE 1 %THEN %RETURN;

%sim_ord_dgp(out=BASE_FILE_SRS, seed=&SEED_CC.);

*long format for lcmm, one row per subject and time;
DATA ord_long;
SET BASE_FILE_SRS;
ARRAY Yarr[&T.] Y1_1-Y1_&T.;
DO t=1 TO &T.;
yval=Yarr[t];
OUTPUT;
END;
KEEP &IDVAR. t yval trueclass;
RUN;

PROC EXPORT DATA=ord_long
  OUTFILE="&OUTPATH.\ordinal_data_for_R.csv"
  DBMS=CSV REPLACE;
RUN;

ODS EXCLUDE ALL;
OPTIONS NONOTES;
%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=&IDVAR.,
  yvars=Y1_1-Y1_&T., tvars=quar1-quar&T., ttotal=&T.,
  m=&MCAT., ycodes=0 1 2 3, deg=&DEG., k=&KTRUE.,
  tech=quanew, maxiter=2000, prefix=cc
);
OPTIONS NOTES;
ODS EXCLUDE NONE;

%ord_coef(prefix=cc, k=&KTRUE., out=_coef);

DATA traj2_EY;
SET _coef;
ARRAY b b_0 b_1 b_2;
DO t=1 TO &T.;
eta=0;
DO j=1 TO DIM(b);
IF NOT MISSING(b[j]) THEN eta+b[j]*(t**(j-1));
END;
p0=PROBNORM(th_1-eta);
p1=PROBNORM(th_2-eta)-PROBNORM(th_1-eta);
p2=PROBNORM(th_3-eta)-PROBNORM(th_2-eta);
p3=1-PROBNORM(th_3-eta);
EY_traj2=0*p0+1*p1+2*p2+3*p3;
intercept_b0=b_0;   *for rank matching the classes across the two programs;
pi_traj2=pie;
OUTPUT;
END;
KEEP classnum t EY_traj2 intercept_b0 pi_traj2;
RUN;

PROC EXPORT DATA=traj2_EY
  OUTFILE="&OUTPATH.\traj2_EY_by_class.csv"
  DBMS=CSV REPLACE;
RUN;

TITLE 'Table 15, Traj2 column. Fit statistics for the cross-check';
TITLE2 'read -2 log L from here and the lcmm value from the R output';
PROC PRINT DATA=work.cc_FIT_K&KTRUE. NOOBS;
RUN;

TITLE 'Table 15, Traj2 column. Class proportions';
PROC SQL;
SELECT classnum, MEAN(pi_traj2) AS pi_traj2 FORMAT=6.3
FROM traj2_EY GROUP BY classnum;
QUIT;
TITLE;
TITLE2;
%MEND run_crosscheck;
%run_crosscheck;


/*##########################################################################################################################
*STEP 5: FIGURES 5 AND 6
  The worked example of paper Section 5.2, N=500.

  cap=3 makes the outcome take the four values 0 1 2 3, which
  matches m=4 and ycodes=0 1 2 3. The BIN4 process gives each class
  a distinct mean trajectory. Leaving start_values= blank starts
  every parameter at zero, which is what Section 5.2 describes.
##########################################################################################################################*/
%MACRO run_example;
%IF &RUN_EXAMPLE. NE 1 %THEN %RETURN;

%sim_data(dist=BIN4, class=3, n=500, T=12, seed=2026, cap=3);
%build_base_from_simwide(T=12);

%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=BENE_ID,
  yvars=Y1_1-Y1_12, tvars=quar1-quar12, ttotal=12,
  m=4, ycodes=0 1 2 3, deg=2, k=3,
  prefix=ordprob_T1C3
);

*Figure 5 is the predicted mean trajectory plot and Figure 6 the
 estimated class proportions;
%ordprob_mix_plot_one(
  outestlib=work, prefix=ordprob_T1C3,
  k=3, ttotal=12, m=4, ycodes=0 1 2 3
);
%MEND run_example;
%run_example;


/*##########################################################################################################################
*STEP 6: QC
##########################################################################################################################*/
%MACRO run_qc;

%IF &RUN_RECOVERY. = 1 %THEN %DO;

TITLE 'QC 1. Table 14. Shape recovery and classification agreement';
TITLE2 'expected: worst-class RMSE mean near 0.029, agreement mean near 0.910';
PROC MEANS DATA=mc_metrics(WHERE=(valid=1)) N MEAN STD MIN MAX MAXDEC=4;
VAR max_rmse agree;
RUN;

TITLE 'QC 2. Table 14. Replications that failed the validity filter';
TITLE2 'this table must be empty for the paper claim that all 200 simulations satisfied the constraints';
PROC PRINT DATA=mc_metrics NOOBS;
WHERE valid NE 1;
VAR rep bijective max_rmse agree valid clean_hess;
RUN;

TITLE 'QC 3. Table 14. Replications with an unusable standard error';
TITLE2 'this table must be empty for the paper claim about fully estimable standard errors';
PROC PRINT DATA=mc_metrics NOOBS;
WHERE clean_hess NE 1;
VAR rep valid clean_hess;
RUN;

TITLE 'QC 4. Table 14. BIC model selection';
TITLE2 "every replication should choose C=&KTRUE.";
PROC FREQ DATA=_bestk;
TABLES kfit / NOCUM;
RUN;

%END;

%IF &RUN_CROSSCHECK. = 1 %THEN %DO;

TITLE 'QC 5. Table 15. Traj2 fitted E[Y_t | class] on the 0-3 scale';
TITLE2 'compare against the lcmm column produced by traj2_lcmm_crosscheck.R';
PROC PRINT DATA=traj2_EY NOOBS;
VAR classnum t EY_traj2;
FORMAT EY_traj2 8.4;
RUN;

TITLE 'QC 6. Cross-check data: simulated class sizes';
TITLE2 "compare with the generating proportions &PA., &PB., &PC.";
PROC FREQ DATA=BASE_FILE_SRS;
TABLES trueclass / NOCUM;
RUN;

%END;

TITLE;
TITLE2;
%MEND run_qc;
%run_qc;


/*##########################################################################################################################
*STEP 7: CLEAN UP WORK LIBRARY
##########################################################################################################################*/
*DELETE INTERMEDIATE FILES;
PROC DATASETS LIB=work NOLIST NOWARN;
DELETE _b _bw _th _thw _alpha _alphaB _mixp _coef _bicrow
       sim_long sim_wide y1_w y2_w;
QUIT;

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT FILES
traj2_ordinal_validation.pdf   STEP 3. Table 14 recovery, model selection and estimability, written to OUTPATH
ordinal_recovery_metrics.csv   STEP 3. One row per replication, the per-replication metrics behind Table 14
ordinal_data_for_R.csv         STEP 4. Long-format simulated data, input to traj2_lcmm_crosscheck.R
traj2_EY_by_class.csv          STEP 4. Traj2 fitted E[Y_t | class] on the 0 to 3 scale, for the Table 15 comparison

RUN ORDER FOR STEP 4
  1. run this script with RUN_CROSSCHECK=1
  2. run traj2_lcmm_crosscheck.R against ordinal_data_for_R.csv
  3. read the lcmm column back into the Table 15 shell printed at the end of STEP 4
##########################################################################################################################*/
