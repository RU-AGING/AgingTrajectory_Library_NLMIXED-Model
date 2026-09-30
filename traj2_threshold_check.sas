/*Copyright (c) 2026 Community Health and Aging Outcomes (CHAO) Lab, Rutgers University.
Released under the MIT License. Full text in LICENSE at
https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model;
/*************************************************************************************************************************
*************************************************************************************************************************
COMMUNITY HEALTH AND AGING OUTCOMES (CHAO) LAB - INSTITUTE FOR HEALTH, HEALTH CARE POLICY & AGING RESEARCH -
RUTGERS, THE STATE UNIVERSITY OF NEW JERSEY
*************************************************************************************************************************
*************************************************************************************************************************
*PROJECT NAME: Traj2 Class-Specific Threshold Check, Ordinal-Probit Outcome Family
LAST UPDATED DATE: 30 SEP 2026 (threshold arrays renamed TT2/TT3, TAU2/TAU3 clashed with variables tau2/tau3)
DATA SOURCES: NONE. All input is simulated in STEP 2 from the data generating process declared in STEP 0
PURPOSE: Reproduces Section 3.3 (Tables 15 and 16, Figure 6) of Traj2: Censored-normal and ordinal-probit
group-based trajectory modeling in SAS for restricted data environments, Zafar, Xia, Lin, Jarrin, MethodsX.
The lcmm cross-check of Section 3.2 uses a DGP with class-common thresholds. This script uses a DGP in
which the threshold SPACING differs by class, and produces
(a)STEP 3, Table 16 and Figure 6, single data set, SEED_CC=20260930. Traj2 fit exported for traj2_threshold_lcmm.R, which fits the
   same data in lcmm (class-common thresholds) and compares both programs against the truth
(b)STEP 4, Table 15, threshold recovery over NREP replications, SEED_MC=20261000. Mean, bias, RMSE and 95 percent
   interval coverage for the class-specific thresholds and the class proportions
COMPANION: traj2_threshold_lcmm.R. Run this script first so the CSV files exist, then run the R script.
IDENTIFICATION: Traj2 fixes the first threshold at 0 in every class and estimates the class intercept freely,
so what is identified per class is the SPACING of the thresholds. The DGP below therefore sets the first
threshold to 0 in every class and varies the spacing: narrow in class B, unit in class A, wide in class C.
AUTHOR: Anum Zafar, Weiyi Xia, Haiqun Lin, Olga F. Jarrin
##########################################################################################################################
*Execution Environment: SAS 9.4 or later with SAS/STAT PROC NLMIXED. No compiled components           *
*ordinal_probit.sas must sit in SRCPATH. OUTPATH is built automatically if its parent folder exists    *
*Runtime: STEP 4 fits one model per replication (no BIC sweep), about 1/6 of traj2_replication_ordinal *
*Smoke test at NREP=5 before committing to the full run                                                *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 0:Parameters. Edit only this block
#STEP 1:Load the macro library
#STEP 2:Simulator with class-specific thresholds, and the truth
#STEP 3:Single data set, SAS half of the Traj2 versus lcmm comparison
#STEP 4:Threshold recovery over NREP replications
#STEP 5:QC
#STEP 6:Clean up the WORK library
##########################################################################################################################;


/*##########################################################################################################################
*STEP 0: PARAMETERS (EDIT THESE ONLY)
##########################################################################################################################*/
%LET SRCPATH=T:\datasets_t;                *folder holding ordinal_probit.sas;
%LET OUTPATH=T:\datasets_t\replication_1;    *folder for the PDF and CSV output;

%LET RUN_SINGLE=1;     *1=run STEP 3, the data set for the lcmm comparison;
%LET RUN_MC=1;         *1=run STEP 4, threshold recovery;

%LET NREP=200;           *replications for STEP 4. Use 5 for a smoke test first;
%LET NSIM=2500;        *subjects per data set, same as Section 3.2;
%LET SEED_CC=20260930; *single seed for STEP 3;
%LET SEED_MC=20261000; *replication r uses seed SEED_MC+r, so 20261001 through 20261200;

%LET KTRUE=3;          *true number of latent classes;
%LET T=12;             *time points;
%LET MCAT=4;           *ordinal categories, coded 0 1 2 3;
%LET DEG=2;            *polynomial degree, quadratic;
%LET IDVAR=BENE_ID;    *subject identifier;

*latent-scale quadratic trajectories. Chosen so every class uses all four categories;
%LET A_B0=0.2;   %LET A_B1=0.08;  %LET A_B2=0.000;   *slow rise;
%LET B_B0=-0.2;  %LET B_B1=0.25;  %LET B_B2=0.000;   *steady rise;
%LET C_B0=0.3;   %LET C_B1=0.10;  %LET C_B2=0.012;   *starts mid, accelerates;

*class-specific thresholds. tau1=0 in every class, only the spacing differs;
%LET A_TAU2=1.0;  %LET A_TAU3=2.0;   *unit spacing;
%LET B_TAU2=0.5;  %LET B_TAU3=1.0;   *narrow spacing, responses pushed to the ends of the scale;
%LET C_TAU2=2.0;  %LET C_TAU3=4.0;   *wide spacing, responses held in the middle of the scale;

*true class proportions, must sum to 1. Same as Section 3.2;
%LET PA=0.45;  %LET PB=0.35;  %LET PC=0.20;

OPTIONS NOFMTERR MPRINT NOMLOGIC NOSYMBOLGEN;

OPTIONS DLCREATEDIR;
LIBNAME _mkdir "&OUTPATH.";
LIBNAME _mkdir CLEAR;
OPTIONS NODLCREATEDIR;


/*##########################################################################################################################
*STEP 1: LOAD THE MACRO LIBRARY
##########################################################################################################################*/
%INCLUDE "&SRCPATH.\ordinal_probit.sas";


/*##########################################################################################################################
*STEP 2: SIMULATOR WITH CLASS-SPECIFIC THRESHOLDS, AND THE TRUTH

  %sim_ord_cst is %sim_ord_dgp from traj2_replication_ordinal.sas
  with one change: the thresholds are looked up by class instead of
  being shared. The draw order is the same, one uniform for the class
  and then one uniform per time point.

  true_par holds the generating thresholds and proportions, one row
  per class. true_traj holds E[Y_t | class] on the 0-3 scale.
##########################################################################################################################*/
%MACRO sim_ord_cst(out=BASE_FILE_SRS, seed=, n=&NSIM., tpts=&T., id=&IDVAR.);
DATA &out.;
CALL STREAMINIT(&seed.);
ARRAY B0[3] _TEMPORARY_ (&A_B0. &B_B0. &C_B0.);
ARRAY B1[3] _TEMPORARY_ (&A_B1. &B_B1. &C_B1.);
ARRAY B2[3] _TEMPORARY_ (&A_B2. &B_B2. &C_B2.);
ARRAY TT2[3] _TEMPORARY_ (&A_TAU2. &B_TAU2. &C_TAU2.);
ARRAY TT3[3] _TEMPORARY_ (&A_TAU3. &B_TAU3. &C_TAU3.);
ARRAY Y[&tpts.] Y1_1-Y1_&tpts.;
ARRAY quar[&tpts.] quar1-quar&tpts.;
cpA=&PA.;
cpB=&PA.+&PB.;
DO &id.=1 TO &n.;
u=RAND('UNIFORM');
IF u<cpA THEN trueclass=1;
ELSE IF u<cpB THEN trueclass=2;
ELSE trueclass=3;
tau1=0;
tau2=TT2[trueclass];
tau3=TT3[trueclass];
DO t=1 TO &tpts.;
quar[t]=t;
eta=B0[trueclass]+B1[trueclass]*t+B2[trueclass]*t*t;
p0=PROBNORM(tau1-eta);
p1=PROBNORM(tau2-eta)-PROBNORM(tau1-eta);
p2=PROBNORM(tau3-eta)-PROBNORM(tau2-eta);
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
%MEND sim_ord_cst;

/*##########################################################################################################################
  SAMPLE DATA OUTPUT DATASET: BASE_FILE_SRS (from %sim_ord_cst)

  BENE_ID trueclass Y1_1 Y1_2 Y1_3 ... Y1_12 quar1 quar2 ... quar12
        1         2    0    3    1 ...     3     1     2 ...     12
        2         3    1    1    2 ...     2     1     2 ...     12

  Same layout as %sim_ord_dgp. trueclass is carried for validation
  only and no fitting macro reads it.
##########################################################################################################################*/

DATA true_par;
ARRAY TT2[3] _TEMPORARY_ (&A_TAU2. &B_TAU2. &C_TAU2.);
ARRAY TT3[3] _TEMPORARY_ (&A_TAU3. &B_TAU3. &C_TAU3.);
ARRAY PP[3] _TEMPORARY_ (&PA. &PB. &PC.);
DO trueclass=1 TO 3;
tau2_true=TT2[trueclass];
tau3_true=TT3[trueclass];
pi_true=PP[trueclass];
OUTPUT;
END;
RUN;

DATA true_traj;
ARRAY B0[3] _TEMPORARY_ (&A_B0. &B_B0. &C_B0.);
ARRAY B1[3] _TEMPORARY_ (&A_B1. &B_B1. &C_B1.);
ARRAY B2[3] _TEMPORARY_ (&A_B2. &B_B2. &C_B2.);
ARRAY TT2[3] _TEMPORARY_ (&A_TAU2. &B_TAU2. &C_TAU2.);
ARRAY TT3[3] _TEMPORARY_ (&A_TAU3. &B_TAU3. &C_TAU3.);
DO trueclass=1 TO 3;
DO t=1 TO &T.;
eta=B0[trueclass]+B1[trueclass]*t+B2[trueclass]*t*t;
p0=PROBNORM(0-eta);
p1=PROBNORM(TT2[trueclass]-eta)-PROBNORM(0-eta);
p2=PROBNORM(TT3[trueclass]-eta)-PROBNORM(TT2[trueclass]-eta);
p3=1-PROBNORM(TT3[trueclass]-eta);
EY_true=0*p0+1*p1+2*p2+3*p3;
OUTPUT;
END;
END;
KEEP trueclass t EY_true;
RUN;

PROC EXPORT DATA=true_traj
  OUTFILE="&OUTPATH.\cst_true_EY.csv"
  DBMS=CSV REPLACE;
RUN;

*%cst_coef pulls betas, thresholds (with their 95 percent limits) and mixing proportions of
 one fit into one row per fitted class, then matches fitted classes to true classes by the
 smallest squared distance between fitted and true E[Y_t | class];
%MACRO cst_coef(prefix=, k=&KTRUE., out=_coef);
%GLOBAL ndist _sumexp;   *read by run_mc after this macro returns;
PROC SQL NOPRINT;
CREATE TABLE _b AS
SELECT UPCASE(SUBSTR(Parameter,6,1)) AS CL LENGTH=1,
       INPUT(SUBSTR(Parameter,7),BEST.) AS IDX, Estimate
FROM work.&prefix._EST_K&k.
WHERE UPCASE(SUBSTR(Parameter,1,5))='BETA_';
QUIT;
PROC SORT DATA=_b; BY CL IDX; RUN;
PROC TRANSPOSE DATA=_b OUT=_bw PREFIX=b_; BY CL; ID IDX; VAR Estimate; RUN;

*thresholds on the natural scale with delta-method 95 percent limits from NLMIXED;
PROC SQL NOPRINT;
CREATE TABLE _th AS
SELECT UPCASE(SCAN(Label,2,'_')) AS CL LENGTH=1,
       INPUT(COMPRESS(SUBSTR(LOWCASE(Label),10),,'kd'),BEST.) AS IDX,
       Estimate, Lower, Upper
FROM work.&prefix._ESTS_K&k.
WHERE UPCASE(SUBSTR(Label,1,9))='THRESHOLD';
QUIT;
PROC SORT DATA=_th; BY CL IDX; RUN;
PROC TRANSPOSE DATA=_th OUT=_thw PREFIX=th_; BY CL; ID IDX; VAR Estimate; RUN;
PROC TRANSPOSE DATA=_th OUT=_thlo PREFIX=lo_; BY CL; ID IDX; VAR Lower; RUN;
PROC TRANSPOSE DATA=_th OUT=_thhi PREFIX=hi_; BY CL; ID IDX; VAR Upper; RUN;

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

DATA _coef0;
MERGE _bw(DROP=_NAME_) _thw(DROP=_NAME_) _thlo(DROP=_NAME_) _thhi(DROP=_NAME_) _mixp;
BY CL;
th_1=0;   *first threshold fixed at 0;
classnum=(CL='A')*1+(CL='B')*2+(CL='C')*3;
RUN;

*fitted E[Y_t | class] and the class match;
DATA _fit_traj;
SET _coef0;
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
FROM _fit_traj f, true_traj tr
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
CREATE TABLE &out. AS
SELECT c.*, m.matched_true AS trueclass, SQRT(m.sse/&T.) AS rmse_EY
FROM _coef0 c LEFT JOIN _match m ON c.classnum=m.fclass
ORDER BY trueclass;
QUIT;
%MEND cst_coef;


/*##########################################################################################################################
*STEP 3: SINGLE DATA SET, SAS HALF OF THE TRAJ2 VERSUS LCMM COMPARISON

  One data set at SEED_CC, fitted with %ordprob_mix_fit_one at C=3
  and exported for traj2_threshold_lcmm.R:

    cst_data_for_R.csv     long-format data, input to lcmm
    cst_traj2_EY.csv       Traj2 fitted E[Y_t | class], matched to the true class
    cst_traj2_fit.csv      Traj2 -2 log L, so the R script can form the likelihood ratio test
    cst_true_EY.csv        true E[Y_t | class], written in STEP 2

  lcmm fixes one threshold set for all classes, which is the
  special case iA2=iB2=iC2 and iA3=iB3=iC3 of the Traj2 model. The two
  models are nested, so -2 log L(lcmm) minus -2 log L(Traj2) is a
  likelihood ratio statistic on (C-1)(M-2)=4 degrees of freedom.
##########################################################################################################################*/
%MACRO run_single;
%IF &RUN_SINGLE. NE 1 %THEN %RETURN;

%sim_ord_cst(out=BASE_FILE_SRS, seed=&SEED_CC.);

DATA cst_long;
SET BASE_FILE_SRS;
ARRAY Yarr[&T.] Y1_1-Y1_&T.;
DO t=1 TO &T.;
yval=Yarr[t];
OUTPUT;
END;
KEEP &IDVAR. t yval trueclass;
RUN;

PROC EXPORT DATA=cst_long
  OUTFILE="&OUTPATH.\cst_data_for_R.csv"
  DBMS=CSV REPLACE;
RUN;

ODS EXCLUDE ALL;
OPTIONS NONOTES;
%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=&IDVAR.,
  yvars=Y1_1-Y1_&T., tvars=quar1-quar&T., ttotal=&T.,
  m=&MCAT., ycodes=0 1 2 3, deg=&DEG., k=&KTRUE.,
  tech=quanew, maxiter=2000, prefix=cs
);
OPTIONS NOTES;
ODS EXCLUDE NONE;

%cst_coef(prefix=cs, out=cs_coef);

DATA cst_traj2_EY;
SET cs_coef;
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
pi_traj2=pie;
OUTPUT;
END;
KEEP trueclass t EY_traj2 pi_traj2;
RUN;

PROC EXPORT DATA=cst_traj2_EY
  OUTFILE="&OUTPATH.\cst_traj2_EY.csv"
  DBMS=CSV REPLACE;
RUN;

DATA cst_traj2_fit;
SET work.cs_FIT_K&KTRUE.;
WHERE UPCASE(Descr) LIKE '-2%';
m2ll_traj2=Value;
KEEP m2ll_traj2;
RUN;

PROC EXPORT DATA=cst_traj2_fit
  OUTFILE="&OUTPATH.\cst_traj2_fit.csv"
  DBMS=CSV REPLACE;
RUN;

*Traj2 thresholds and proportions for the single data set;
DATA cst_table;
MERGE true_par cs_coef(KEEP=trueclass th_2 th_3 lo_2 hi_2 lo_3 hi_3 pie rmse_EY);
BY trueclass;
RUN;
%MEND run_single;
%run_single;

/*##########################################################################################################################
  SAMPLE DATA OUTPUT DATASET: cst_table (one row per true class)

  trueclass tau2_true tau3_true pi_true  th_2  lo_2  hi_2  th_3  lo_3  hi_3   pie  rmse_EY
          1       1.0       2.0    0.45  1.02  0.97  1.07  2.01  1.93  2.09  0.45    0.010
          2       0.5       1.0    0.35  0.49  0.46  0.52  0.99  0.94  1.04  0.35    0.012
          3       2.0       4.0    0.20  1.98  1.85  2.11  3.95  3.66  4.24  0.20    0.015
##########################################################################################################################*/


/*##########################################################################################################################
*STEP 4: THRESHOLD RECOVERY OVER NREP REPLICATIONS

  One fit per replication at C=3. For each true class the script
  records the estimated tau2 and tau3, whether their 95 percent
  intervals cover the truth, the estimated class proportion and the
  RMSE of fitted E[Y_t | class] against the truth. A replication
  counts only if the fitted-to-true class map is one to one.
##########################################################################################################################*/
%MACRO run_mc;
%IF &RUN_MC. NE 1 %THEN %RETURN;

%LOCAL r seed;
%GLOBAL ndrop_cst;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE cst_mc;
QUIT;

ODS GRAPHICS OFF;
ODS EXCLUDE ALL;
ODS NORESULTS;
OPTIONS NONOTES NOSOURCE NOSOURCE2 NOMPRINT;

%DO r=1 %TO &NREP.;
%LET seed=%EVAL(&SEED_MC. + &r.);
%sim_ord_cst(out=BASE_FILE_SRS, seed=&seed.);

%ordprob_mix_fit_one(
  data=BASE_FILE_SRS, id=&IDVAR.,
  yvars=Y1_1-Y1_&T., tvars=quar1-quar&T., ttotal=&T.,
  m=&MCAT., ycodes=0 1 2 3, deg=&DEG., k=&KTRUE.,
  tech=quanew, maxiter=2000, prefix=cm
);

%IF %SYSFUNC(EXIST(work.cm_EST_K&KTRUE.)) %THEN %DO;
%cst_coef(prefix=cm, out=_mcrow0);
DATA _mcrow;
MERGE true_par _mcrow0(KEEP=trueclass th_2 th_3 lo_2 hi_2 lo_3 hi_3 pie rmse_EY);
BY trueclass;
rep=&r.;
bijective=(&ndist.=&KTRUE.);
cov_2=(lo_2<=tau2_true<=hi_2);
cov_3=(lo_3<=tau3_true<=hi_3);
RUN;
PROC APPEND BASE=cst_mc DATA=_mcrow FORCE; RUN;
%END;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE cm_EST_K&KTRUE. cm_ESTS_K&KTRUE. cm_FIT_K&KTRUE. _mcrow _mcrow0;
QUIT;

%IF %SYSFUNC(MOD(&r.,25))=0 %THEN %PUT ===== completed replication &r. of &NREP. =====;
%END;

OPTIONS NOTES SOURCE SOURCE2 MPRINT;
ODS EXCLUDE NONE;
ODS RESULTS;
ODS GRAPHICS ON;

PROC SQL NOPRINT;
SELECT COUNT(DISTINCT rep) INTO :ndrop_cst TRIMMED FROM cst_mc WHERE bijective=0;
QUIT;

*one row per true class: truth, mean, bias, RMSE and coverage;
PROC SQL;
CREATE TABLE cst_mc_summary AS
SELECT trueclass,
       MEAN(tau2_true) AS tau2_true, MEAN(th_2) AS tau2_mean, MEAN(th_2-tau2_true) AS tau2_bias,
       SQRT(MEAN((th_2-tau2_true)**2)) AS tau2_rmse, MEAN(cov_2) AS tau2_cov,
       MEAN(tau3_true) AS tau3_true, MEAN(th_3) AS tau3_mean, MEAN(th_3-tau3_true) AS tau3_bias,
       SQRT(MEAN((th_3-tau3_true)**2)) AS tau3_rmse, MEAN(cov_3) AS tau3_cov,
       MEAN(pi_true) AS pi_true, MEAN(pie) AS pi_mean, MEAN(rmse_EY) AS rmse_EY,
       COUNT(*) AS nrep
FROM cst_mc
WHERE bijective=1
GROUP BY trueclass;
QUIT;

ODS PDF FILE="&OUTPATH.\traj2_threshold_recovery.pdf" STARTPAGE=NO;
TITLE "Class-specific threshold recovery (&NREP. reps, N=&NSIM., T=&T., C=&KTRUE.)";
TITLE2 "tau1=0 in every class; coverage is for Wald 95 percent intervals; reps dropped as non-bijective: &ndrop_cst.";
PROC PRINT DATA=cst_mc_summary NOOBS;
FORMAT tau2_true tau2_mean tau3_true tau3_mean pi_true pi_mean 6.3
       tau2_bias tau2_rmse tau3_bias tau3_rmse rmse_EY 7.4 tau2_cov tau3_cov 6.3;
RUN;
TITLE;
TITLE2;
ODS PDF CLOSE;

PROC EXPORT DATA=cst_mc
  OUTFILE="&OUTPATH.\cst_recovery_metrics.csv"
  DBMS=CSV REPLACE;
RUN;
%MEND run_mc;
%run_mc;

/*##########################################################################################################################
  SAMPLE DATA OUTPUT DATASET: cst_mc_summary (one row per true class)

  trueclass tau2_true tau2_mean tau2_bias tau2_rmse tau2_cov tau3_true ... pi_true pi_mean rmse_EY nrep
          1     1.000     1.001    0.0012    0.0290    0.950     2.000 ...   0.450   0.449  0.0110  200
##########################################################################################################################*/


/*##########################################################################################################################
*STEP 5: QC
##########################################################################################################################*/
%MACRO run_qc;

%IF &RUN_SINGLE. = 1 %THEN %DO;

TITLE 'QC 1. Single data set. Traj2 thresholds and proportions against the truth';
TITLE2 'each true value should sit inside its interval; the three fitted classes must map to three different true classes';
PROC PRINT DATA=cst_table NOOBS;
FORMAT tau2_true tau3_true pi_true th_2 lo_2 hi_2 th_3 lo_3 hi_3 pie 6.3 rmse_EY 7.4;
RUN;

TITLE 'QC 2. Single data set. Traj2 -2 log L';
TITLE2 'carried into traj2_threshold_lcmm.R for the likelihood ratio test';
PROC PRINT DATA=cst_traj2_fit NOOBS;
RUN;

%END;

%IF &RUN_MC. = 1 %THEN %DO;

TITLE 'QC 3. Recovery. Replications that failed the one to one class map';
TITLE2 'this table should be empty';
PROC PRINT DATA=cst_mc NOOBS;
WHERE bijective NE 1;
VAR rep trueclass bijective rmse_EY;
RUN;

TITLE 'QC 4. Recovery. Coverage by class';
TITLE2 'expected: every coverage rate within about 0.92 to 0.98 (binomial SE 0.015 at 200 reps)';
PROC PRINT DATA=cst_mc_summary NOOBS;
VAR trueclass nrep tau2_cov tau3_cov;
RUN;

%END;

TITLE;
TITLE2;
%MEND run_qc;
%run_qc;


/*##########################################################################################################################
*STEP 6: CLEAN UP WORK LIBRARY
##########################################################################################################################*/
*DELETE INTERMEDIATE FILES;
PROC DATASETS LIB=work NOLIST NOWARN;
DELETE _b _bw _th _thw _thlo _thhi _alpha _alphaB _mixp _coef0 _fit_traj _pair _match
       cst_long;
QUIT;

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT FILES, all written to OUTPATH
cst_true_EY.csv                  STEP 2. True E[Y_t | class] on the 0-3 scale
cst_data_for_R.csv               STEP 3. Long-format simulated data, input to traj2_threshold_lcmm.R
cst_traj2_EY.csv                 STEP 3. Traj2 fitted E[Y_t | class] and class proportions, by true class
cst_traj2_fit.csv                STEP 3. Traj2 -2 log L
traj2_threshold_recovery.pdf     STEP 4. Recovery table
cst_recovery_metrics.csv         STEP 4. One row per replication and true class

RUN ORDER
  1. run this script with NREP=5 as a smoke test, then NREP=200
  2. run traj2_threshold_lcmm.R against the STEP 3 CSV files
##########################################################################################################################*/

