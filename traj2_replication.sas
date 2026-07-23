*PROJECT NAME: Traj2 Replication, Continuous Censored-Normal Outcome Family
LAST UPDATED DATE: 23 JUL 2026
DATA SOURCES: NONE. All input is simulated in STEP 2 from the data generating process declared in STEP 0
PURPOSE: Single standalone script reproducing the continuous-outcome tables and figures in
Traj2: A Native Macro Library for Single and Multi-Outcome Group-Based Trajectory Modeling in SAS
Zafar, Xia, Lin, Jarrin, Journal of Statistical Software. Produces
(a)STEP 3, Tables 9, 11, 12, 13 and Figures 3 and 4. Validation run, N=10,000, SEED_VAL=20260609
(b)STEP 4, Table 10. Monte Carlo bias and coverage, 200 replications, N=2,500, SEED_MC=770000
(c)STEP 5, Table 17 and Figure 7. Worked example of paper Sections 5.2 and 5.3, N=500
NOT COVERED: the ordinal-probit results, paper Sections 3.6, 4.4, 4.5, Tables 14 and 15, Figures 1, 2, 5, 6.
Those belong to the companion script traj2_replication_ordinal.sas
AUTHOR: Anum Zafar, Weiyi Xia, Haiqun Lin, Olga F. Jarrin
##########################################################################################################################
*Execution Environment: SAS 9.4 or later with SAS/STAT PROC NLMIXED and SAS/IML. No compiled components*
*continuous_2outcomes.sas must sit in SRCPATH. OUTPATH and FIGPATH are built if their parent exists    *
*The PROC TRAJ blocks in STEP 3 need the compiled PROC TRAJ binary and will not run inside the CMS VRDC*
*Runtime: STEP 3 about 5 minutes, STEP 4 about 60 minutes at NREP=200, STEP 5 under 1 minute           *
*Smoke test at NREP=5 before committing to the full run                                                *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 0:Parameters. Edit only this block
#STEP 1:Load the macro library and define the figure style
#STEP 2:Shared simulator and the generating values
#STEP 3:Tables 9, 11, 12, 13 and Figures 3 and 4. Validation run against PROC TRAJ
#STEP 4:Table 10. Monte Carlo bias and coverage
#STEP 5:Table 17 and Figure 7. Worked example
#STEP 6:Check output
#STEP 7:Clean up the WORK library
##########################################################################################################################;


/*##########################################################################################################################
*STEP 0: PARAMETERS (EDIT THESE ONLY)
##########################################################################################################################*/
%LET SRCPATH=T:\datasets_t;                *folder holding continuous_2outcomes.sas;
%LET OUTPATH=T:\datasets_t\replication;    *folder for the PDF and CSV output;
%LET FIGPATH=T:\datasets_t\figs;           *folder for the standalone vector figures;

%LET RUN_VALIDATION=1;   *1=run STEP 3, Tables 9 11 12 13 and Figures 3 and 4;
%LET RUN_TRAJ=1;         *1=also run the PROC TRAJ half of STEP 3. Set to 0 inside the VRDC;
%LET RUN_MC=1;           *1=run STEP 4, Table 10. This is the slow one;
%LET RUN_EXAMPLE=1;      *1=run STEP 5, Table 17 and Figure 7;

%LET SEED_VAL=20260609;  *seed behind Tables 9 11 12 13 and Figures 3 and 4;
%LET NVAL=10000;         *subjects in the validation data set;

%LET NREP=200;           *Monte Carlo replications. Use 5 for a smoke test first;
%LET NSUB=2500;          *subjects per Monte Carlo replication;
%LET SEED_MC=770000;     *replication r uses seed SEED_MC+r, so 770001 through 770200;

*Data generating process for STEP 3 and STEP 4, described in paper Section 4.1.
 Coefficients are listed class by class as b0 b1 b2;
%LET DGP_LC=3;                                                    *latent classes;
%LET DGP_DEG=2;                                                   *polynomial degree, quadratic;
%LET DGP_PI=0.40 0.35 0.25;                                       *mixing proportions;
%LET DGP_Y1=8.0 -0.30 0.000  20.0 4.00 -0.100  40.0 6.00 0.000;   *HH, skilled home based care;
%LET DGP_Y2=5.0 0.50 0.000  30.0 -1.00 0.050  60.0 2.00 0.000;    *INP, inpatient care;
%LET DGP_SIG=6;                                                   *residual SD, both outcomes;
%LET DGP_CAP=100;                                                 *upper censoring bound, lower bound is 0;
%LET CAPLIST=100 100 100 100 100 100 100 100 100 100 100 100;     *the same cap repeated T times;

*Starting values. Class intercepts are spread across the outcome range and the
 slope and curvature terms start at zero. Written out in full rather than
 generated, so the exact numbers behind the published fits are on the page;
*The mixing intercepts are held separately from the outcome blocks. The joint
 fit needs both outcome blocks but only ONE copy of the alphas: listing them
 twice makes PROC NLMIXED stop with "Parameter alpha0_B has a duplicate
 specification" and no ParameterEstimates table is written;
%LET ALPHA0=alpha0_B=0 alpha0_C=0;

%LET BETA_Y1=beta1_A0=10 beta1_A1=0 beta1_A2=0
             beta1_B0=25 beta1_B1=0 beta1_B2=0
             beta1_C0=45 beta1_C1=0 beta1_C2=0
             sigma1_=6;

%LET BETA_Y2=beta2_A0=5  beta2_A1=0 beta2_A2=0
             beta2_B0=30 beta2_B1=0 beta2_B2=0
             beta2_C0=60 beta2_C1=0 beta2_C2=0
             sigma2_=6;

%LET START_Y1=&ALPHA0. &BETA_Y1.;              *single outcome fit, HH;
%LET START_Y2=&ALPHA0. &BETA_Y2.;              *single outcome fit, INP;
%LET MC_START=&ALPHA0. &BETA_Y1. &BETA_Y2.;    *joint fit, alphas listed once;

OPTIONS NOFMTERR MPRINT NOMLOGIC NOSYMBOLGEN;

*Create OUTPATH and FIGPATH if they are not there yet. Without this an
 ODS PDF FILE= opens against a missing folder and only fails later, at
 ODS PDF CLOSE, with "Physical file does not exist". The parent folder
 must already exist;
OPTIONS DLCREATEDIR;
LIBNAME _mkdir "&OUTPATH.";
LIBNAME _mkdir CLEAR;
LIBNAME _mkdir "&FIGPATH.";
LIBNAME _mkdir CLEAR;
OPTIONS NODLCREATEDIR;


/*##########################################################################################################################
*STEP 1: LOAD THE LIBRARY AND DEFINE THE FIGURE STYLE

  TRAJ2_RUN_DEMO=0 loads the macro definitions and runs nothing,
  which is the behaviour paper Section 5.1 describes.

  The WORK template store is session temporary, so defining
  styles.jss here keeps the script self contained. For colour
  figures change the parent to styles.statistical.
##########################################################################################################################*/
%LET TRAJ2_RUN_DEMO=0;
%INCLUDE "&SRCPATH.\continuous_2outcomes.sas";

PROC TEMPLATE;
DEFINE STYLE styles.jss;
PARENT=styles.journal;
CLASS GraphFonts /
  'GraphTitleText'   = ("Times New Roman", 11pt, bold)
  'GraphLabelText'   = ("Times New Roman", 10pt)
  'GraphValueText'   = ("Times New Roman",  9pt)
  'GraphDataText'    = ("Times New Roman",  9pt)
  'GraphFootnoteText'= ("Times New Roman",  8pt);
CLASS GraphData1 / linethickness=2px markersize=7px;
CLASS GraphData2 / linethickness=2px markersize=7px;
CLASS GraphData3 / linethickness=2px markersize=7px;
END;
RUN;


/*##########################################################################################################################
*STEP 2: SHARED SIMULATOR

  %sim_cont_dgp builds the paper Section 4.1 validation DGP in the
  wide contract of Table 1, CAP columns included. One generator now
  feeds both the validation run and the Monte Carlo study, instead
  of the DGP being written out again in each script.

  It draws random numbers in the same order as the original inline
  DATA steps (one uniform per subject, then Y1 and Y2 normals
  alternating within the time loop), so the streams are identical
  and the published numbers reproduce.
##########################################################################################################################*/
%MACRO sim_cont_dgp(out=BASE_FILE_SRS, n=2500, T=12, seed=20260609,
                    lc=3, deg=2, pilist=%STR(0.40 0.35 0.25),
                    y1coef=, y2coef=, sigma=6, cap=100,
                    y1pre=QHH, y2pre=QINP, tpre=quar, idvar=BENE_ID);
DATA &out.;
CALL STREAMINIT(&seed.);
ARRAY B1[&lc.,%EVAL(&deg.+1)] _TEMPORARY_ (&y1coef.);
ARRAY B2[&lc.,%EVAL(&deg.+1)] _TEMPORARY_ (&y2coef.);
ARRAY PP[&lc.] _TEMPORARY_ (&pilist.);
ARRAY Y1[&T.] &y1pre.1-&y1pre.&T.;
ARRAY Y2[&T.] &y2pre.1-&y2pre.&T.;
ARRAY TT[&T.] &tpre.1-&tpre.&T.;
ARRAY CP[&T.] CAP1-CAP&T.;
DO &idvar.=1 TO &n.;
u=RAND('UNIFORM');
cum=0;
TRUECLASS=.;
DO c=1 TO &lc.;
cum=cum+PP[c];
IF TRUECLASS=. AND u<cum THEN TRUECLASS=c;
END;
IF TRUECLASS=. THEN TRUECLASS=&lc.;
DO j=1 TO &T.;
TT[j]=j;
CP[j]=&cap.;
m1=0;
m2=0;
DO d=0 TO &deg.;
m1=m1+B1[TRUECLASS,d+1]*(j**d);
m2=m2+B2[TRUECLASS,d+1]*(j**d);
END;
Y1[j]=MIN(MAX(m1+RAND('NORMAL')*&sigma.,0),&cap.);
Y2[j]=MIN(MAX(m2+RAND('NORMAL')*&sigma.,0),&cap.);
END;
OUTPUT;
END;
KEEP &idvar. TRUECLASS &y1pre.1-&y1pre.&T. &y2pre.1-&y2pre.&T. &tpre.1-&tpre.&T. CAP1-CAP&T.;
RUN;
%MEND sim_cont_dgp;

/*##########################################################################################################################
  SAMPLE DATA OUTPUT DATASET: BASE_FILE_SRS (from %sim_cont_dgp)

  BENE_ID TRUECLASS QHH1 QHH2 ... QHH12 QINP1 ... QINP12 quar1 ... quar12 CAP1 ... CAP12
        1         2 21.4 26.9 ...  85.1  27.3 ...   25.0     1 ...     12  100 ...    100
        2         1  8.9  4.1 ...   0.0   6.2 ...   10.4     1 ...     12  100 ...    100

  One row per subject, one column per outcome per time point.
  TRUECLASS is carried for validation only and no fitting macro reads it.
##########################################################################################################################*/

*The generating values, used by STEP 3 and STEP 4. The true mixing
 logits are the log odds against class A:
   alpha0_B = LOG(0.35/0.40) = -0.133531
   alpha0_C = LOG(0.25/0.40) = -0.470004;
DATA dgp_truth;
LENGTH Parameter $20;
INPUT Parameter $ true_val ord;
DATALINES;
alpha0_B -0.133531 1
alpha0_C -0.470004 2
beta1_A0 8.0 3
beta1_A1 -0.30 4
beta1_A2 0.00 5
beta1_B0 20.0 6
beta1_B1 4.00 7
beta1_B2 -0.10 8
beta1_C0 40.0 9
beta1_C1 6.00 10
beta1_C2 0.00 11
sigma1_ 6.0 12
beta2_A0 5.0 13
beta2_A1 0.50 14
beta2_A2 0.00 15
beta2_B0 30.0 16
beta2_B1 -1.00 17
beta2_B2 0.05 18
beta2_C0 60.0 19
beta2_C1 2.00 20
beta2_C2 0.00 21
sigma2_ 6.0 22
;
RUN;

*truth by intercept rank, for the class matched comparison tables;
DATA truth_y1;
INPUT cls_true $ b0_t b1_t b2_t;
DATALINES;
c1 8.0 -0.30 0.000
c2 20.0 4.00 -0.100
c3 40.0 6.00 0.000
;
RUN;

DATA truth_pi;
rk=1; pi_true=0.25; OUTPUT;
rk=2; pi_true=0.35; OUTPUT;
rk=3; pi_true=0.40; OUTPUT;
RUN;


/*##########################################################################################################################
*STEP 3: TABLES 9, 11, 12, 13 AND FIGURES 3 AND 4
  The validation run, N=10,000, seed 20260609.

  PART A   Traj2 single outcome fit on HH, recovery against truth
  PART B   PROC TRAJ on the same outcome, coefficient and proportion agreement
  PART C   joint two-outcome fit via the warm start pipeline, Table 9, Figure 3
  PART D   PROC TRAJ multi-trajectory, Table 12
  PART E   twin predicted versus observed panel, supplementary

  Classes are matched between engines by intercept rank throughout,
  so "rank" means the same thing in every table and figure. In this
  DGP the lowest intercept class carries the largest proportion, so
  intercept rank and proportion rank run in opposite directions.
  Keying everything to intercept avoids the column flip that causes.
##########################################################################################################################*/
%MACRO run_validation;
%IF &RUN_VALIDATION. NE 1 %THEN %RETURN;

%LET T=12;
%LET class=&DGP_LC.;
%LET order_model=&DGP_DEG.;
%LET equal=T;
%LET max_values=&CAPLIST.;
%LET y1vars=QHH1-QHH12;
%LET y2vars=QINP1-QINP12;
%LET idvar=BENE_ID;
%LET tvars=quar1-quar12;

%sim_cont_dgp(out=BASE_FILE_SRS, n=&NVAL., T=&T., seed=&SEED_VAL.,
              lc=&DGP_LC., deg=&DGP_DEG., pilist=%STR(&DGP_PI.),
              y1coef=%STR(&DGP_Y1.), y2coef=%STR(&DGP_Y2.),
              sigma=&DGP_SIG., cap=&DGP_CAP.);

*keep a copy for the QC section, since STEP 5 overwrites BASE_FILE_SRS;
DATA val_data;
SET BASE_FILE_SRS;
RUN;

ODS GRAPHICS ON / RESET WIDTH=6.5IN HEIGHT=4IN IMAGEFMT=PDF ANTIALIAS=ON;
ODS PDF FILE="&OUTPATH.\traj2_comparison_results.pdf" STYLE=styles.jss STARTPAGE=YES;
ODS LISTING GPATH="&FIGPATH." STYLE=styles.jss;
ODS EXCLUDE IterHistory;   *drop the long NLMIXED iteration table, keep estimates and fit statistics;

/*##########################################################################################################################
  PART A: Traj2 single outcome fit on HH and parameter recovery
##########################################################################################################################*/
ODS OUTPUT FitStatistics=t2_fit_y1;
%nlmixed_1(T=&T., LC=&class., Y=1, starting=&START_Y1.,
           output=fit_y1, order=&order_model., equal_sigma=&equal.);
ODS OUTPUT CLOSE;

*estimated betas, long to one row per class holding b0 b1 b2;
DATA t2_beta;
SET fit_y1;
WHERE INDEX(parameter,'beta1_')=1;
cls=SUBSTR(parameter,7,1);
ord=INPUT(SUBSTR(parameter,8,2),2.);
est=estimate;
KEEP cls ord est;
RUN;
PROC SORT DATA=t2_beta; BY cls ord; RUN;
PROC TRANSPOSE DATA=t2_beta OUT=t2_w(DROP=_name_) PREFIX=b;
BY cls;
ID ord;
VAR est;
RUN;

PROC RANK DATA=t2_w OUT=t2_rank DESCENDING; VAR b0; RANKS rk; RUN;
PROC RANK DATA=truth_y1 OUT=truth_rank DESCENDING; VAR b0_t; RANKS rk; RUN;
PROC SORT DATA=t2_rank; BY rk; RUN;
PROC SORT DATA=truth_rank; BY rk; RUN;

DATA recovery_y1;
MERGE truth_rank t2_rank;
BY rk;
bias_b0=b0-b0_t;
bias_b1=b1-b1_t;
bias_b2=b2-b2_t;
RUN;

TITLE 'Table A. Parameter recovery, outcome 1 (true vs Traj2, classes matched by intercept)';
PROC PRINT DATA=recovery_y1 NOOBS LABEL;
VAR cls_true cls b0_t b0 bias_b0 b1_t b1 bias_b1 b2_t b2 bias_b2;
LABEL cls_true='TrueCls' cls='FitCls'
      b0_t='b0(true)' b0='b0(est)' bias_b0='bias'
      b1_t='b1(true)' b1='b1(est)' bias_b1='bias'
      b2_t='b2(true)' b2='b2(est)' bias_b2='bias';
FORMAT b0_t b0 bias_b0 b1_t b1 bias_b1 b2_t b2 bias_b2 8.4;
RUN;
TITLE;

*proportion recovery, keyed to intercept rank;
PROC SQL NOPRINT;
SELECT estimate INTO :aB TRIMMED FROM fit_y1 WHERE parameter='alpha0_B';
SELECT estimate INTO :aC TRIMMED FROM fit_y1 WHERE parameter='alpha0_C';
QUIT;

DATA t2_pi;
LENGTH cls $1;
eA=1;
eB=EXP(&aB.);
eC=EXP(&aC.);
s=eA+eB+eC;
cls='A'; pi=eA/s; OUTPUT;
cls='B'; pi=eB/s; OUTPUT;
cls='C'; pi=eC/s; OUTPUT;
KEEP cls pi;
RUN;

PROC SQL;
CREATE TABLE prop_cmp AS
SELECT r.rk, t.pi AS pi_fit
FROM t2_rank r JOIN t2_pi t ON r.cls=t.cls;
QUIT;

PROC SORT DATA=prop_cmp; BY rk; RUN;
DATA prop_cmp;
MERGE truth_pi prop_cmp;
BY rk;
bias=pi_fit-pi_true;
RUN;

TITLE 'Table A2. Class-proportion recovery (by intercept rank): Traj2 vs truth';
PROC PRINT DATA=prop_cmp NOOBS LABEL;
VAR rk pi_true pi_fit bias;
LABEL rk='Rank' pi_true='pi(true)' pi_fit='pi(Traj2)' bias='bias';
FORMAT pi_true pi_fit bias 6.3;
RUN;
TITLE;

/*##########################################################################################################################
  PART B: PROC TRAJ on the same outcome, then agreement.
  This is paper Table 11 (coefficients), Table 13 (proportions)
  and Figure 4 (fitted class trajectories).

  PROC TRAJ output schema, confirmed:
    traj_out  : BENE_ID, QHH1-12, quar1-12, GRP1PRB-GRP3PRB, GROUP
    traj_stat : per group polynomial coefficients BETA0-BETA5 plus PI
    traj_est  : parameter estimates
  PROC TRAJ does not write BIC to a data set. It prints BIC and the
  log likelihood in the Results window, and that output does not
  route to ODS PDF, so read it off the screen.
##########################################################################################################################*/
%IF &RUN_TRAJ. = 1 %THEN %DO;

PROC TRAJ DATA=BASE_FILE_SRS OUT=traj_out OUTEST=traj_est OUTSTAT=traj_stat;
ID &idvar.;
VAR QHH1-QHH12;
INDEP quar1-quar12;
MODEL CNORM;
MIN 0;
MAX &DGP_CAP.;
NGROUPS &class.;
ORDER 2 2 2;
RUN;

PROC RANK DATA=traj_stat OUT=traj_coef DESCENDING; VAR BETA0; RANKS rk; RUN;
PROC RANK DATA=t2_w OUT=t2_coef DESCENDING; VAR b0; RANKS rk; RUN;
PROC SORT DATA=traj_coef; BY rk; RUN;
PROC SORT DATA=t2_coef; BY rk; RUN;

DATA coef_cmp;
MERGE t2_coef(KEEP=rk b0 b1 b2 RENAME=(b0=b0_t2 b1=b1_t2 b2=b2_t2))
      traj_coef(KEEP=rk BETA0 BETA1 BETA2 RENAME=(BETA0=b0_tr BETA1=b1_tr BETA2=b2_tr));
BY rk;
d0=b0_t2-b0_tr;
d1=b1_t2-b1_tr;
d2=b2_t2-b2_tr;
RUN;

TITLE 'Table 11. Coefficient agreement: Traj2 vs PROC TRAJ (classes matched by intercept)';
PROC PRINT DATA=coef_cmp NOOBS LABEL;
VAR rk b0_t2 b0_tr d0 b1_t2 b1_tr d1 b2_t2 b2_tr d2;
LABEL rk='Class' b0_t2='b0 Traj2' b0_tr='b0 TRAJ' d0='diff'
      b1_t2='b1 Traj2' b1_tr='b1 TRAJ' d1='diff'
      b2_t2='b2 Traj2' b2_tr='b2 TRAJ' d2='diff';
FORMAT b0_t2 b0_tr b1_t2 b1_tr b2_t2 b2_tr 9.5 d0 d1 d2 e10.;
RUN;
TITLE;

DATA traj_props;
SET traj_coef;
pi_traj=PI/100;
KEEP rk pi_traj;
RUN;
PROC SORT DATA=traj_props; BY rk; RUN;

DATA prop_all;
MERGE prop_cmp traj_props;
BY rk;
RUN;

TITLE 'Table 13. Class proportions (by intercept rank): truth vs Traj2 vs PROC TRAJ';
PROC PRINT DATA=prop_all NOOBS LABEL;
VAR rk pi_true pi_fit pi_traj;
LABEL rk='Rank' pi_true='pi(true)' pi_fit='pi(Traj2)' pi_traj='pi(PROC TRAJ)';
FORMAT pi_true pi_fit pi_traj 6.3;
RUN;
TITLE;

DATA t2_bic;
SET t2_fit_y1;
WHERE INDEX(UPCASE(descr),'BIC')>0;
engine='Traj2';
bic=value;
KEEP engine bic;
RUN;

TITLE 'Traj2 BIC (read the PROC TRAJ BIC off its Results window; equal fit means equal log-L)';
PROC PRINT DATA=t2_bic NOOBS;
FORMAT bic 12.1;
RUN;
TITLE;

*Figure 4, fitted class trajectories on the latent scale.
 Wide form, one row per class and quarter with one column per source,
 because WHERE is not valid on a SERIES or SCATTER statement;
PROC SORT DATA=truth_rank; BY rk; RUN;
DATA curve_src;
MERGE coef_cmp truth_rank(KEEP=rk b0_t b1_t b2_t);
BY rk;
class=rk;
DO quar=1 TO &T.;
mu_true=b0_t+b1_t*quar+b2_t*quar*quar;
mu_t2=b0_t2+b1_t2*quar+b2_t2*quar*quar;
mu_tr=b0_tr+b1_tr*quar+b2_tr*quar*quar;
OUTPUT;
END;
KEEP class quar mu_true mu_t2 mu_tr;
RUN;

TITLE 'Figure 4. Fitted class trajectories on the latent scale';
TITLE2 'lines = truth, filled circles = Traj2, x = PROC TRAJ (markers sit on the lines)';
ODS GRAPHICS ON / IMAGENAME="fig_traj_overlay";
PROC SGPLOT DATA=curve_src;
SERIES X=quar Y=mu_true / GROUP=class LINEATTRS=(THICKNESS=2) NAME="ln";
SCATTER X=quar Y=mu_t2 / GROUP=class MARKERATTRS=(SYMBOL=CIRCLEFILLED SIZE=9);
SCATTER X=quar Y=mu_tr / GROUP=class MARKERATTRS=(SYMBOL=X SIZE=11);
KEYLEGEND "ln" / TITLE="Latent class";
XAXIS LABEL="Quarter (t)";
YAXIS LABEL="Mean trajectory (latent scale)";
RUN;
TITLE;
TITLE2;

%END;

/*##########################################################################################################################
  PART C: joint two-outcome fit via the warm start pipeline.
  This is paper Table 9 and Figure 3.
##########################################################################################################################*/
%nlmixed_1(T=&T., LC=&class., Y=2, starting=&START_Y2.,
           output=fit_y2, order=&order_model., equal_sigma=&equal.);

DATA nlm_2y_starting;
SET fit_y1 fit_y2;
IF parameter=:'alpha' THEN DELETE;
KEEP parameter estimate;
RUN;

%nlmixed_MultiTraj(T=&T., LC=&class.,
                   starting=%starting_value_alpha(class=&class.) / data=nlm_2y_starting,
                   output=fit_joint, order=&order_model., equal_sigma=&equal.);

*Table 9, joint estimates beside the generating values;
PROC SQL;
CREATE TABLE table9 AS
SELECT t.ord, t.Parameter, t.true_val, e.Estimate, e.StandardError,
       (e.Estimate - t.true_val) AS deviation
FROM dgp_truth t LEFT JOIN fit_joint e ON t.Parameter = e.Parameter
ORDER BY t.ord;
QUIT;

*Figure 3, predicted versus posterior weighted observed, both outcomes;
ODS GRAPHICS ON / IMAGENAME="fig_plotprep";
%plot_prep(T=&T., LC=&class., result=fit_joint,
           order=&order_model., equal_sigma=&equal.);

*mixing proportions from the joint fit, for the Table 9 footnote;
PROC SQL NOPRINT;
SELECT estimate INTO :jB TRIMMED FROM fit_joint WHERE parameter='alpha0_B';
SELECT estimate INTO :jC TRIMMED FROM fit_joint WHERE parameter='alpha0_C';
QUIT;

DATA joint_pi;
LENGTH class $1;
eA=1;
eB=EXP(&jB.);
eC=EXP(&jC.);
s=eA+eB+eC;
class='A'; pi=eA/s; OUTPUT;
class='B'; pi=eB/s; OUTPUT;
class='C'; pi=eC/s; OUTPUT;
KEEP class pi;
RUN;

TITLE 'Figure. Estimated class proportions (joint two-outcome fit)';
ODS GRAPHICS ON / IMAGENAME="fig_proportions";
PROC SGPLOT DATA=joint_pi;
VBAR class / RESPONSE=pi DATALABEL;
YAXIS LABEL="Mixing proportion";
XAXIS LABEL="Latent class";
RUN;
TITLE;

/*##########################################################################################################################
  PART D: PROC TRAJ multi-trajectory on both outcomes.
  This is paper Table 12.

  Note the statement pattern: the FIRST outcome uses the unnumbered
  VAR / INDEP / MODEL / MIN / MAX / ORDER statements and the SECOND
  uses the numbered VAR2 / INDEP2 / MODEL2 / MIN2 / MAX2 / ORDER2,
  with MULTGROUPS giving the number of shared latent classes.
##########################################################################################################################*/
%IF &RUN_TRAJ. = 1 %THEN %DO;

PROC TRAJ DATA=BASE_FILE_SRS OUT=trajj_out OUTEST=trajj_est OUTSTAT=trajj_stat
          OUTPLOT=trajj_plot OUTSTAT2=trajj_stat2 OUTPLOT2=trajj_plot2;
ID &idvar.;
VAR QHH1-QHH12;
INDEP quar1-quar12;
MODEL CNORM;
MIN 0;
MAX &DGP_CAP.;
ORDER 2 2 2;
VAR2 QINP1-QINP12;
INDEP2 quar1-quar12;
MODEL2 CNORM;
MIN2 0;
MAX2 &DGP_CAP.;
ORDER2 2 2 2;
MULTGROUPS &class.;
RUN;

*PROC TRAJ side, stacked over the two outcomes and ranked by the
 outcome 1 intercept, so both outcomes carry the same class label;
DATA trajj_o1;
SET trajj_stat;
grp=_N_;
outc=1;
b0=BETA0; b1=BETA1; b2=BETA2;
KEEP grp outc b0 b1 b2;
RUN;

DATA trajj_o2;
SET trajj_stat2;
grp=_N_;
outc=2;
b0=BETA0; b1=BETA1; b2=BETA2;
KEEP grp outc b0 b1 b2;
RUN;

PROC RANK DATA=trajj_o1 OUT=trajj_rk DESCENDING; VAR b0; RANKS rk; RUN;

DATA trajj_both;
SET trajj_o1 trajj_o2;
RUN;

PROC SQL;
CREATE TABLE trajj_cmp AS
SELECT a.outc, m.rk, a.b0 AS b0_tr, a.b1 AS b1_tr, a.b2 AS b2_tr
FROM trajj_both a JOIN (SELECT grp, rk FROM trajj_rk) m ON a.grp=m.grp;
QUIT;

*Traj2 side, same treatment;
DATA t2j_beta;
SET fit_joint;
IF INDEX(parameter,'beta1_')=1 THEN outc=1;
ELSE IF INDEX(parameter,'beta2_')=1 THEN outc=2;
ELSE DELETE;
cls=SUBSTR(parameter,INDEX(parameter,'_')+1,1);
ord=INPUT(SUBSTR(parameter,INDEX(parameter,'_')+2,2),2.);
est=estimate;
KEEP outc cls ord est;
RUN;
PROC SORT DATA=t2j_beta; BY outc cls ord; RUN;
PROC TRANSPOSE DATA=t2j_beta OUT=t2j_w(DROP=_name_) PREFIX=b;
BY outc cls;
ID ord;
VAR est;
RUN;

PROC RANK DATA=t2j_w(WHERE=(outc=1)) OUT=t2j_rk DESCENDING; VAR b0; RANKS rk; RUN;

PROC SQL;
CREATE TABLE t2j_cmp AS
SELECT w.outc, m.rk, w.b0 AS b0_t2, w.b1 AS b1_t2, w.b2 AS b2_t2
FROM t2j_w w JOIN (SELECT cls, rk FROM t2j_rk) m ON w.cls=m.cls;
QUIT;

PROC SORT DATA=t2j_cmp; BY outc rk; RUN;
PROC SORT DATA=trajj_cmp; BY outc rk; RUN;

DATA table12;
MERGE t2j_cmp trajj_cmp;
BY outc rk;
LENGTH outcome $3;
IF outc=1 THEN outcome='HH';
ELSE outcome='INP';
d0=b0_t2-b0_tr;
d1=b1_t2-b1_tr;
d2=b2_t2-b2_tr;
RUN;

TITLE 'Table 12. Joint two-outcome coefficient agreement: Traj2 vs PROC TRAJ';
TITLE2 'classes matched by the outcome 1 intercept rank';
PROC PRINT DATA=table12 NOOBS LABEL;
VAR outcome rk b0_t2 b0_tr d0 b1_t2 b1_tr d1 b2_t2 b2_tr d2;
LABEL outcome='Outcome' rk='Class' b0_t2='b0 Traj2' b0_tr='b0 TRAJ' d0='diff'
      b1_t2='b1 Traj2' b1_tr='b1 TRAJ' d1='diff'
      b2_t2='b2 Traj2' b2_tr='b2 TRAJ' d2='diff';
FORMAT b0_t2 b0_tr b1_t2 b1_tr b2_t2 b2_tr 9.5 d0 d1 d2 e10.;
RUN;
TITLE;
TITLE2;

/*##########################################################################################################################
  PART E: twin predicted versus observed panel, single outcome HH.
  Supplementary figure, not in the paper. Both engines are drawn in
  the %plot_prep style: predicted censored mean curve per class plus
  posterior weighted observed averages, classes aligned by intercept
  rank so the two panels are directly comparable.

  Coefficient level agreement is deliberately not automated for the
  figure. PROC TRAJ may rescale time internally, so its raw betas are
  not always directly comparable. The fitted mean trajectory is
  invariant to that rescaling, so the comparison is done there.
##########################################################################################################################*/
PROC SQL NOPRINT;
SELECT estimate INTO :sig TRIMMED FROM fit_y1 WHERE parameter='sigma1_';
SELECT rk INTO :rkA TRIMMED FROM t2_coef WHERE cls='A';
SELECT rk INTO :rkB TRIMMED FROM t2_coef WHERE cls='B';
SELECT rk INTO :rkC TRIMMED FROM t2_coef WHERE cls='C';
QUIT;

PROC TRANSPOSE DATA=fit_y1 OUT=p_y1(DROP=_name_);
ID parameter;
VAR estimate;
RUN;

DATA t2_post;
IF _N_=1 THEN SET p_y1;
SET BASE_FILE_SRS;
ARRAY Y[12] QHH1-QHH12;
ARRAY X[12] quar1-quar12;
alpha0_A=0;
LA=0; LB=0; LC=0;
DO t=1 TO 12;
mA=beta1_A0+beta1_A1*X[t]+beta1_A2*X[t]**2;
mB=beta1_B0+beta1_B1*X[t]+beta1_B2*X[t]**2;
mC=beta1_C0+beta1_C1*X[t]+beta1_C2*X[t]**2;
eA=Y[t]-mA; eB=Y[t]-mB; eC=Y[t]-mC;
IF Y[t]=0 THEN DO;
LA+LOGCDF('NORMAL',eA/sigma1_);
LB+LOGCDF('NORMAL',eB/sigma1_);
LC+LOGCDF('NORMAL',eC/sigma1_);
END;
ELSE IF Y[t]=&DGP_CAP. THEN DO;
LA+LOGCDF('NORMAL',-eA/sigma1_);
LB+LOGCDF('NORMAL',-eB/sigma1_);
LC+LOGCDF('NORMAL',-eC/sigma1_);
END;
ELSE DO;
LA+(LOGPDF('NORMAL',eA/sigma1_)-LOG(sigma1_));
LB+(LOGPDF('NORMAL',eB/sigma1_)-LOG(sigma1_));
LC+(LOGPDF('NORMAL',eC/sigma1_)-LOG(sigma1_));
END;
END;
pA=EXP(alpha0_A); pB=EXP(alpha0_B); pC=EXP(alpha0_C); pd=pA+pB+pC;
lnA=LOG(pA/pd)+LA; lnB=LOG(pB/pd)+LB; lnC=LOG(pC/pd)+LC;
mx=MAX(lnA,lnB,lnC);
sden=EXP(lnA-mx)+EXP(lnB-mx)+EXP(lnC-mx);
postA=EXP(lnA-mx)/sden;
postB=EXP(lnB-mx)/sden;
postC=EXP(lnC-mx)/sden;
KEEP QHH1-QHH12 postA postB postC;
RUN;

DATA t2_long;
SET t2_post;
ARRAY Q[12] QHH1-QHH12;
ARRAY P[3] postA postB postC;
ARRAY RK[3] _TEMPORARY_ (&rkA. &rkB. &rkC.);
DO c=1 TO 3;
DO t=1 TO 12;
rank=RK[c]; quar=t; y=Q[t]; w=P[c];
OUTPUT;
END;
END;
KEEP rank quar y w;
RUN;

PROC SUMMARY DATA=t2_long NWAY;
CLASS rank quar;
WEIGHT w;
VAR y;
OUTPUT OUT=t2_obs(KEEP=rank quar obs) MEAN=obs;
RUN;

DATA t2_pred;
SET p_y1;
ARRAY B0[3] beta1_A0 beta1_B0 beta1_C0;
ARRAY B1[3] beta1_A1 beta1_B1 beta1_C1;
ARRAY B2[3] beta1_A2 beta1_B2 beta1_C2;
ARRAY RK[3] _TEMPORARY_ (&rkA. &rkB. &rkC.);
DO c=1 TO 3;
DO quar=1 TO 12;
mu=B0[c]+B1[c]*quar+B2[c]*quar*quar;
a=(0-mu)/sigma1_;
bb=(&DGP_CAP.-mu)/sigma1_;
pred=mu*(CDF('normal',bb)-CDF('normal',a))
     +sigma1_*(PDF('normal',a)-PDF('normal',bb))
     +&DGP_CAP.*(1-CDF('normal',bb));
rank=RK[c];
OUTPUT;
END;
END;
KEEP rank quar pred;
RUN;

PROC SORT DATA=t2_obs; BY rank quar; RUN;
PROC SORT DATA=t2_pred; BY rank quar; RUN;
DATA t2_panel;
MERGE t2_pred t2_obs;
BY rank quar;
LENGTH engine $9;
engine='Traj2';
RUN;

DATA traj_stat_g;
SET traj_stat;
grp=_N_;
RUN;
PROC RANK DATA=traj_stat_g OUT=traj_coef2 DESCENDING; VAR BETA0; RANKS rk; RUN;
PROC SQL NOPRINT;
SELECT rk INTO :rg1 TRIMMED FROM traj_coef2 WHERE grp=1;
SELECT rk INTO :rg2 TRIMMED FROM traj_coef2 WHERE grp=2;
SELECT rk INTO :rg3 TRIMMED FROM traj_coef2 WHERE grp=3;
QUIT;

DATA traj_long;
SET traj_out;
ARRAY Q[12] QHH1-QHH12;
ARRAY P[3] GRP1PRB GRP2PRB GRP3PRB;
ARRAY RG[3] _TEMPORARY_ (&rg1. &rg2. &rg3.);
DO g=1 TO 3;
DO t=1 TO 12;
rank=RG[g]; quar=t; y=Q[t]; w=P[g];
OUTPUT;
END;
END;
KEEP rank quar y w;
RUN;

PROC SUMMARY DATA=traj_long NWAY;
CLASS rank quar;
WEIGHT w;
VAR y;
OUTPUT OUT=traj_obs2(KEEP=rank quar obs) MEAN=obs;
RUN;

DATA traj_pred2;
SET traj_coef2;
DO quar=1 TO 12;
mu=BETA0+BETA1*quar+BETA2*quar*quar;
a=(0-mu)/&sig.;
bb=(&DGP_CAP.-mu)/&sig.;
pred=mu*(CDF('normal',bb)-CDF('normal',a))
     +&sig.*(PDF('normal',a)-PDF('normal',bb))
     +&DGP_CAP.*(1-CDF('normal',bb));
rank=rk;
OUTPUT;
END;
KEEP rank quar pred;
RUN;

PROC SORT DATA=traj_obs2; BY rank quar; RUN;
PROC SORT DATA=traj_pred2; BY rank quar; RUN;
DATA traj_panel;
MERGE traj_pred2 traj_obs2;
BY rank quar;
LENGTH engine $9;
engine='PROC TRAJ';
RUN;

DATA twin_panel;
SET t2_panel traj_panel;
RUN;

TITLE 'Figure. Predicted (line) vs posterior-weighted observed (markers) by class';
TITLE2 'QHH single-outcome fit: Traj2 vs PROC TRAJ (classes aligned by intercept rank)';
ODS GRAPHICS ON / IMAGENAME="fig_twin_panel";
PROC SGPANEL DATA=twin_panel;
PANELBY engine / COLUMNS=2 NOVARNAME;
SERIES X=quar Y=pred / GROUP=rank;
SCATTER X=quar Y=obs / GROUP=rank;
COLAXIS LABEL="Quarter (t)";
ROWAXIS LABEL="QHH";
RUN;
TITLE;
TITLE2;

%END;

ODS LISTING CLOSE;
ODS PDF CLOSE;
ODS GRAPHICS OFF;
ODS LISTING;
%MEND run_validation;
%run_validation;


/*##########################################################################################################################
*STEP 4: TABLE 10
  Monte Carlo bias and empirical 95 percent CI coverage.

  NREP independent data sets from the same DGP as STEP 3, each at
  a moderate N. Every replication is fitted from one fixed set of
  scale matched starting values, so the replications are mutually
  independent and none of them uses the warm start pipeline.

  Coverage is deliberately assessed at a moderate N. At N=10,000 the
  standard errors shrink so far that a trivial finite sample bias can
  fall outside the interval and make coverage look broken when the
  estimator is fine.

  With 200 replications the Monte Carlo SE of a coverage rate is
  about SQRT(0.95*0.05/200) = 0.015, so rates roughly between 0.92
  and 0.98 are consistent with the nominal 0.95.
##########################################################################################################################*/
%MACRO run_mc;
%IF &RUN_MC. NE 1 %THEN %RETURN;

%LOCAL r seed;

%LET T=12;
%LET class=&DGP_LC.;
%LET order_model=&DGP_DEG.;
%LET equal=T;
%LET max_values=&CAPLIST.;
%LET y1vars=QHH1-QHH12;
%LET y2vars=QINP1-QINP12;
%LET idvar=BENE_ID;
%LET tvars=quar1-quar12;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE mc_all mc_conv;
QUIT;

ODS GRAPHICS OFF;
ODS EXCLUDE ALL;
ODS NORESULTS;
OPTIONS NONOTES NOSOURCE NOSOURCE2 NOMPRINT;

%DO r=1 %TO &NREP.;
%LET seed=%EVAL(&SEED_MC. + &r.);

%sim_cont_dgp(out=BASE_FILE_SRS, n=&NSUB., T=&T., seed=&seed.,
              lc=&DGP_LC., deg=&DGP_DEG., pilist=%STR(&DGP_PI.),
              y1coef=%STR(&DGP_Y1.), y2coef=%STR(&DGP_Y2.),
              sigma=&DGP_SIG., cap=&DGP_CAP.);

ODS OUTPUT ConvergenceStatus=cs;
%nlmixed_MultiTraj(T=&T., LC=&class., starting=&MC_START.,
                   output=pe_tmp, order=&order_model., equal_sigma=&equal.);
ODS OUTPUT CLOSE;

%IF %SYSFUNC(EXIST(pe_tmp)) %THEN %DO;
DATA pe_tmp;
SET pe_tmp;
rep=&r.;
RUN;
PROC APPEND BASE=mc_all DATA=pe_tmp FORCE;
RUN;
%END;

%IF %SYSFUNC(EXIST(cs)) %THEN %DO;
DATA cs;
SET cs;
rep=&r.;
RUN;
PROC APPEND BASE=mc_conv DATA=cs FORCE;
RUN;
%END;

PROC DATASETS LIB=work NOLIST NOWARN;
DELETE pe_tmp cs;
QUIT;
%END;

OPTIONS NOTES SOURCE SOURCE2 MPRINT;
ODS EXCLUDE NONE;
ODS RESULTS;
ODS GRAPHICS ON;

*Keep only the replications that converged. Status=0 is the machine
 readable converged flag. PROC NLMIXED still writes a
 ParameterEstimates table when it stops early, so without this filter
 a failed fit would enter the bias and coverage averages unnoticed;
PROC SQL;
CREATE TABLE mc_ok AS
SELECT DISTINCT rep FROM mc_conv WHERE Status = 0;
QUIT;

PROC SQL;
CREATE TABLE mc_join AS
SELECT a.rep, a.Parameter, a.Estimate, a.Lower, a.Upper, t.true_val, t.ord
FROM mc_all a
INNER JOIN dgp_truth t ON a.Parameter = t.Parameter
INNER JOIN mc_ok k ON a.rep = k.rep;
QUIT;

*A replication with no usable standard error is counted, not silently
 dropped. Dropping it without counting would push coverage upward;
DATA mc_join;
SET mc_join;
bias = Estimate - true_val;
IF Lower=. OR Upper=. THEN DO;
covered=.;
noSE=1;
END;
ELSE DO;
covered=(Lower <= true_val <= Upper);
noSE=0;
END;
RUN;

PROC SUMMARY DATA=mc_join NWAY;
CLASS Parameter ord;
VAR Estimate bias covered noSE;
OUTPUT OUT=mc_summary(DROP=_type_ _freq_)
       MEAN(Estimate)=mean_est
       MEAN(bias)=mean_bias
       MEAN(covered)=coverage
       N(bias)=n_reps
       SUM(noSE)=n_noSE;
RUN;

PROC SQL;
CREATE TABLE table10 AS
SELECT s.ord, s.Parameter, t.true_val, s.mean_est, s.mean_bias,
       s.coverage, s.n_reps, s.n_noSE
FROM mc_summary s JOIN dgp_truth t ON s.Parameter = t.Parameter
ORDER BY s.ord;
QUIT;

ODS PDF FILE="&OUTPATH.\traj2_mc_coverage.pdf" STARTPAGE=NO;
TITLE 'Table 10. Traj2 Monte Carlo recovery and coverage (joint two-outcome model)';
TITLE2 "&NREP. replications, N=&NSUB. per replication; empirical 95% CI coverage";
PROC PRINT DATA=table10 NOOBS LABEL;
VAR Parameter true_val mean_est mean_bias coverage n_reps n_noSE;
LABEL Parameter='Parameter' true_val='True' mean_est='Mean est'
      mean_bias='Mean bias' coverage='95% coverage' n_reps='n reps'
      n_noSE='n reps with no SE';
FORMAT true_val mean_est mean_bias 8.4 coverage 6.3;
RUN;

TITLE 'Convergence status across replications';
TITLE2 'every replication behind Table 10 must show Status=0';
PROC FREQ DATA=mc_conv;
TABLES Status*Reason / LIST MISSING;
RUN;
TITLE;
TITLE2;
ODS PDF CLOSE;

PROC EXPORT DATA=table10
  OUTFILE="&OUTPATH.\mc_coverage_summary.csv"
  DBMS=CSV REPLACE;
RUN;
%MEND run_mc;
%run_mc;


/*##########################################################################################################################
*STEP 5: TABLE 17 AND FIGURE 7
  The worked example of paper Sections 5.2 and 5.3, N=500.

  This section uses the library simulator %sim_data_cont through
  %build_base_file_srs, which is the code path Table 16 documents,
  and the neutral starting values paper Section 5.2 describes.
  Keeping the published starting values is what lets Table 17
  reproduce.

  %build_base_file_srs hard codes seed=1 in its call to
  %sim_data_cont. That seed is what pins Table 17 and Figure 7.
##########################################################################################################################*/
%MACRO run_example;
%IF &RUN_EXAMPLE. NE 1 %THEN %RETURN;

%LET USE_SIM=1;                                          *1=simulate, 0=use your own BASE_FILE_SRS;
%LET T=12;                                               *time points;
%LET class=3;                                            *latent classes;
%LET order_model=2;                                      *1=linear, 2=quadratic, 3=cubic;
%LET equal=T;                                            *T=equal residual SD across classes;
%LET max_values=10 10 10 10 10 10 10 10 10 10 10 10;     *censoring caps, one per time point;
%LET y1vars=QHH1-QHH12;                                  *outcome 1 repeated measures;
%LET y2vars=QINP1-QINP12;                                *outcome 2 repeated measures;
%LET idvar=BENE_ID;                                      *subject identifier;
%LET tvars=quar1-quar12;                                 *time index variables;

%build_base_file_srs;                                    *calls %sim_data_cont with seed=1;

*STEP 3 closed and reopened ODS LISTING, which drops the GPATH, so
 point it back at FIGPATH before the plotting macro runs;
ODS GRAPHICS ON / RESET WIDTH=6.5IN HEIGHT=4IN IMAGEFMT=PDF ANTIALIAS=ON;
ODS PDF FILE="&OUTPATH.\traj2_example_results.pdf" STYLE=styles.jss STARTPAGE=YES;
ODS LISTING GPATH="&FIGPATH." STYLE=styles.jss;
%nlmixed_1(T=&T., LC=&class., Y=1,
           starting=%starting_value_alpha(class=&class.)
                    %starting_value_beta_sigma(class=&class., outcome=1,
                                               order=&order_model., equal_sigma=&equal.),
           output=nlm_fix_T1C3, order=&order_model., equal_sigma=&equal.);

%nlmixed_1(T=&T., LC=&class., Y=2,
           starting=%starting_value_alpha(class=&class.)
                    %starting_value_beta_sigma(class=&class., outcome=2,
                                               order=&order_model., equal_sigma=&equal.),
           output=nlm_fix_T2C3, order=&order_model., equal_sigma=&equal.);

*If either single outcome fit above stalls at the symmetric starting
 point, add int_lo=0, int_hi=4, sigma0=1 to %starting_value_beta_sigma,
 which is scale matched to this simulator;

DATA work.nlm_2y_starting;
SET nlm_fix_T1C3 nlm_fix_T2C3;
IF parameter=:'alpha' THEN DELETE;
KEEP parameter estimate;
RUN;

*Table 17 is the fit statistics table PROC NLMIXED prints for this call;
%nlmixed_MultiTraj(T=&T., LC=&class.,
                   starting=%starting_value_alpha(class=&class.) / data=work.nlm_2y_starting,
                   output=nlm_fix_T1_T2C3, order=&order_model., equal_sigma=&equal.);
*Figure 7 and the averaged posterior class membership summary;
%plot_prep(T=&T., LC=&class., result=nlm_fix_T1_T2C3,
           order=&order_model., equal_sigma=&equal.);

ODS LISTING CLOSE;
ODS PDF CLOSE;
ODS GRAPHICS OFF;
ODS LISTING;
%MEND run_example;
%run_example;

/*##########################################################################################################################
*STEP 6: QC
##########################################################################################################################*/
%MACRO run_qc;

%IF &RUN_VALIDATION. = 1 %THEN %DO;

TITLE 'QC 1. Table 9. Recovery of the generating values, joint fit, N=10,000';
TITLE2 'expected: intercepts within about 0.12, sigma near 6.02 and 6.03';
PROC PRINT DATA=table9 NOOBS;
VAR Parameter true_val Estimate StandardError deviation;
FORMAT true_val Estimate StandardError deviation 10.4;
RUN;

TITLE 'QC 2. Table 9. Fitted mixing proportions from the joint fit';
TITLE2 'expected: 0.388, 0.362, 0.250 against generating 0.400, 0.350, 0.250';
PROC PRINT DATA=joint_pi NOOBS;
FORMAT pi 6.4;
RUN;

TITLE 'QC 3. Censoring check on the validation data';
TITLE2 'both bounds must carry mass, otherwise the censored normal reduces to a plain normal';
PROC MEANS DATA=val_data N NMISS MIN MAX MEAN;
VAR QHH1 QHH12 QINP1 QINP12;
RUN;

TITLE 'QC 4. Simulated class sizes on the validation data';
TITLE2 'compare with the generating mixing proportions 0.40, 0.35, 0.25';
PROC FREQ DATA=val_data;
TABLES TRUECLASS / NOCUM;
RUN;

%END;

%IF &RUN_MC. = 1 %THEN %DO;

TITLE 'QC 5. Table 10. Monte Carlo bias and coverage';
TITLE2 'n reps should equal NREP and n reps with no SE should be 0 for every parameter';
PROC PRINT DATA=table10 NOOBS;
VAR Parameter true_val mean_est mean_bias coverage n_reps n_noSE;
FORMAT true_val mean_est mean_bias 8.4 coverage 6.3;
RUN;

TITLE 'QC 6. Table 10. Replications that did not converge';
TITLE2 'this table must be empty for the paper claim that all replications met GCONV';
PROC PRINT DATA=mc_conv NOOBS;
WHERE Status NE 0;
VAR rep Status Reason;
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
DELETE sim_long sim_wide _y1_w _y2_w parameter data_pred pred_temp pred pred2
       wide avg avg_y1 avg_y2 avg_y1_plot avg_y2_plot y_1 y_2 data_plot
       mc_join mc_summary mc_ok pe_tmp cs
       t2_beta t2_w t2_rank truth_rank t2_pi t2_coef traj_coef traj_props
       trajj_o1 trajj_o2 trajj_both trajj_rk t2j_beta t2j_w t2j_rk
       t2_post t2_long t2_obs t2_pred t2_panel p_y1
       traj_stat_g traj_coef2 traj_long traj_obs2 traj_pred2 traj_panel;
QUIT;

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT FILES
traj2_comparison_results.pdf   STEP 3. Tables 9, 11, 12, 13 and the PROC TRAJ comparison, written to OUTPATH
traj2_mc_coverage.pdf          STEP 4. Table 10 bias and coverage, written to OUTPATH
mc_coverage_summary.csv        STEP 4. The same coverage summary as a flat file
traj2_example_results.pdf      STEP 5. Table 17 and Figure 7, written to OUTPATH
fig_traj_overlay.png           STEP 3. Figure 4, written to FIGPATH
fig_plotprep.png               STEP 3. Figure 3, written to FIGPATH
fig_proportions.png            STEP 3. Estimated class proportions, written to FIGPATH
fig_twin_panel.png             STEP 4. Monte Carlo twin panel, written to FIGPATH
##########################################################################################################################*/
