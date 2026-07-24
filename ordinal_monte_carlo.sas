*************************************************************************************************************************
*************************************************************************************************************************
COMMUNITY HEALTH AND AGING OUTCOMES (CHAO) LAB - INSTITUTE FOR HEALTH, HEALTH CARE POLICY & AGING RESEARCH - 
RUTGERS, THE STATE UNIVERSITY OF NEW JERSEY                  
*************************************************************************************************************************
*************************************************************************************************************************
*PROJECT NAME: Traj2 Ordinal-Probit Validation
LAST UPDATED DATE: 23 JUL 2026
DATA SOURCES: NONE. All input is simulated in STEP 1 from the data generating process declared in STEP 0
PURPOSE: Reproduces the ordinal-probit recovery results for the Traj2 paper. Produces
(a)worst-class trajectory RMSE against the known data generating process
(b)class-assignment agreement under maximum posterior membership
(c)the number of latent classes chosen by minimum BIC over a C=KMIN to C=KMAX sweep
Requires the released ORDINAL_PROBIT.SAS, in which the first threshold is fixed at 0 and the class
intercept BETA_*0 is freely estimated. The threshold parameters are the increments I*2 onward.
AUTHOR: Anum Zafar
##########################################################################################################################
*Execution Environment: SAS 9.4 or later, PROC NLMIXED from SAS/STAT. No compiled components                         *
*Set NREP to a small value first to confirm timing before running the full replication count                         *
##########################################################################################################################
### CODE OVERVIEW ##
#STEP 0:Macro library, data generating process constants, run controls
#STEP 1:Probit simulator. One replication of BASE_FILE_SRS
#STEP 2:Shape-recovery truth. Expected outcome by class and time under the DGP
#STEP 3:Validation loop. Fit at the true class count, shape recovery, classification, BIC sweep
#STEP 4:Check output
##########################################################################################################################;

/*##########################################################################################################################
#SAS MACRO LIBRARY
##########################################################################################################################*/
OPTIONS FULLSTIMER;
%LET MACLIB = T:\datasets_t;                /* folder holding ORDINAL_PROBIT.SAS. Use a relative path for release */
%INCLUDE "&MACLIB.\ordinal_probit.sas";

/*##########################################################################################################################
GLOBAL PARAMETERS:EDIT ONLY THIS BLOCK
##########################################################################################################################*/
*Class-specific latent means. ETA = B0 + B1*T + B2*T*T over the T time points;
%LET A_B0=-1.0;*class A, low and flat;
%LET A_B1=0.00;
%LET A_B2=0.000;
%LET B_B0=-0.8;*class B, moderate rise;
%LET B_B1=0.18;
%LET B_B2=0.000;
%LET C_B0=-0.3;*class C, starts mid and accelerates;
%LET C_B1=0.10;
%LET C_B2=0.012;
*Class-common threshold increments. TAU1=I1, TAU2=TAU1+EXP(I2), TAU3=TAU2+EXP(I3);
%LET I1=0.4;
%LET I2=-0.1;
%LET I3=0.0;
*True class proportions. Must sum to 1;
%LET PA=0.45;
%LET PB=0.35;
%LET PC=0.20;
*Run controls;
%LET NSIM=2500;*subjects per replication;
%LET NREP=200;*replications;
%LET SEED0=20260000;*base seed. Replication R uses SEED0+R;
%LET KTRUE=3;*true number of latent classes;
%LET KMIN=2;*BIC sweep low;
%LET KMAX=6;*BIC sweep high;
%LET T=12;*time points;
%LET MCAT=4;*ordinal categories, coded 0 to 3;
%LET DEG=2;*polynomial degree for the class mean;
%LET IDVAR=BENE_ID;*subject identifier;
%LET TECH=QUANEW;*NLMIXED optimizer;
%LET MAXITER=2000;*NLMIXED iteration cap;
%LET OUTDIR=%SYSFUNC(PATHNAME(WORK));*WORK is wiped when the session ends. Point this at a permanent folder to keep the PDF;
%PUT NOTE: Ordinal validation. &NREP. replications at N=&NSIM. per replication.;

/*##########################################################################################################################
*STEP 1: PROBIT SIMULATOR
##########################################################################################################################
Builds one replication of BASE_FILE_SRS from the DGP in STEP 0. Category probabilities come from the
cumulative probit with class-common thresholds. TRUECLASS is carried so STEP 3 can score classification.

SAMPLE DATA OUTPUT DATASET:
BENE_ID   TRUECLASS   Y1_1   Y1_2   ...   Y1_12   QUAR1   QUAR2   ...   QUAR12
1         1           0      0            1       1       2             12
##########################################################################################################################*/
%MACRO SIM_ONE(SEED=);
DATA BASE_FILE_SRS;
CALL STREAMINIT(&SEED.);
ARRAY B0[3] _TEMPORARY_ (&A_B0. &B_B0. &C_B0.);
ARRAY B1[3] _TEMPORARY_ (&A_B1. &B_B1. &C_B1.);
ARRAY B2[3] _TEMPORARY_ (&A_B2. &B_B2. &C_B2.);
ARRAY Y[&T.] Y1_1-Y1_&T.;
ARRAY QUAR[&T.] QUAR1-QUAR&T.;
TAU1=&I1.;
TAU2=TAU1+EXP(&I2.);
TAU3=TAU2+EXP(&I3.);
CPA=&PA.;*cumulative class probabilities;
CPB=&PA.+&PB.;
DO &IDVAR. = 1 TO &NSIM.;
U=RAND('UNIFORM');
IF U<CPA THEN TRUECLASS=1;
ELSE IF U<CPB THEN TRUECLASS=2;
ELSE TRUECLASS=3;
DO T=1 TO &T.;
QUAR[T]=T;
ETA=B0[TRUECLASS]+B1[TRUECLASS]*T+B2[TRUECLASS]*T*T;
P0=PROBNORM(TAU1-ETA);
P1=PROBNORM(TAU2-ETA)-PROBNORM(TAU1-ETA);
P2=PROBNORM(TAU3-ETA)-PROBNORM(TAU2-ETA);
P3=1-PROBNORM(TAU3-ETA);
UU=RAND('UNIFORM');
IF UU<P0 THEN Y[T]=0;
ELSE IF UU<P0+P1 THEN Y[T]=1;
ELSE IF UU<P0+P1+P2 THEN Y[T]=2;
ELSE Y[T]=3;
END;
OUTPUT;
END;
KEEP &IDVAR. TRUECLASS Y1_1-Y1_&T. QUAR1-QUAR&T.;
RUN;
%MEND SIM_ONE;

/*##########################################################################################################################
*STEP 2: SHAPE-RECOVERY TRUTH
##########################################################################################################################
Expected outcome on the 0 to 3 scale for each class and time point under the DGP. STEP 3 matches each
fitted class to a true class by minimum squared distance against this table, then reports the worst class.

SAMPLE DATA OUTPUT DATASET:
TRUECLASS   T   EY_TRUE
1           1   0.0918
2           1   0.1827
3           1   0.3523
##########################################################################################################################*/
DATA TRUE_TRAJ;
ARRAY B0[3] _TEMPORARY_ (&A_B0. &B_B0. &C_B0.);
ARRAY B1[3] _TEMPORARY_ (&A_B1. &B_B1. &C_B1.);
ARRAY B2[3] _TEMPORARY_ (&A_B2. &B_B2. &C_B2.);
T1=&I1.;
T2=T1+EXP(&I2.);
T3=T2+EXP(&I3.);
DO TRUECLASS=1 TO 3;
DO T=1 TO &T.;
ETA=B0[TRUECLASS]+B1[TRUECLASS]*T+B2[TRUECLASS]*T*T;
P0=PROBNORM(T1-ETA);
P1=PROBNORM(T2-ETA)-PROBNORM(T1-ETA);
P2=PROBNORM(T3-ETA)-PROBNORM(T2-ETA);
P3=1-PROBNORM(T3-ETA);
EY_TRUE = 0*P0 + 1*P1 + 2*P2 + 3*P3;
OUTPUT;
END;
END;
KEEP TRUECLASS T EY_TRUE;
RUN;

/*##########################################################################################################################
*STEP 3: VALIDATION LOOP
##########################################################################################################################
One pass per replication. Fits at the true class count, scores trajectory recovery and classification,
then sweeps C=KMIN to C=KMAX to record the BIC-selected class count. Intermediate tables are deleted at
the foot of each replication so WORK does not grow across the run.

SAMPLE DATA OUTPUT DATASET (MC_METRICS):
REP   BIJECTIVE   MAX_RMSE   AGREE    VALID   CLEAN_HESS
1     1           0.0317     0.9028   1       1
##########################################################################################################################*/
PROC DATASETS LIB=WORK NOLIST NOWARN; DELETE MC_METRICS MC_BICPICK; QUIT;
ODS GRAPHICS OFF;
ODS EXCLUDE ALL;
OPTIONS NONOTES;

%MACRO RUN_VALIDATION;
%DO R=1 %TO &NREP.;
%LET SEED=%EVAL(&SEED0. + &R.);
%SIM_ONE(SEED=&SEED.);

*Fit at the true class count;
%ORDPROB_MIX_FIT_ONE(
DATA=BASE_FILE_SRS, ID=&IDVAR.,
YVARS=Y1_1-Y1_&T., TVARS=QUAR1-QUAR&T., TTOTAL=&T.,
M=&MCAT., YCODES=0 1 2 3, DEG=&DEG., K=&KTRUE.,
TECH=&TECH., MAXITER=&MAXITER., PREFIX=OA
);

%IF %SYSFUNC(EXIST(WORK.OA_EST_K&KTRUE.)) %THEN %DO;

PROC SQL NOPRINT;
CREATE TABLE _B AS
SELECT UPCASE(SUBSTR(PARAMETER,6,1)) AS CL LENGTH=1,
INPUT(SUBSTR(PARAMETER,7),BEST.) AS IDX, ESTIMATE
FROM WORK.OA_EST_K&KTRUE.
WHERE UPCASE(SUBSTR(PARAMETER,1,5))='BETA_';
QUIT;
PROC TRANSPOSE DATA=_B OUT=_BW PREFIX=B_; BY CL; ID IDX; VAR ESTIMATE; RUN;

*ESTS holds threshold2 to threshold(m-1). Threshold1 is fixed at 0 in the released macro and is set below;
PROC SQL NOPRINT;
CREATE TABLE _TH AS
SELECT UPCASE(SCAN(LABEL,2,'_')) AS CL LENGTH=1,
INPUT(COMPRESS(SUBSTR(LOWCASE(LABEL),10),,'kd'),BEST.) AS IDX, ESTIMATE
FROM WORK.OA_ESTS_K&KTRUE.
WHERE UPCASE(SUBSTR(LABEL,1,9))='THRESHOLD';
QUIT;
PROC TRANSPOSE DATA=_TH OUT=_THW PREFIX=TH_; BY CL; ID IDX; VAR ESTIMATE; RUN;

*Class proportions from the alpha0_* estimates of this fit. Class A is fixed at 0 and is not returned;
DATA _ALPHA;
LENGTH CL $1 ALPHA 8;
CL='A';
ALPHA=0;
RUN;

DATA _ALPHAB;
SET WORK.OA_EST_K&KTRUE.(KEEP=PARAMETER ESTIMATE);
LENGTH CL $1 ALPHA 8;
IF UPCASE(SUBSTR(PARAMETER,1,7))='ALPHA0_';
CL=UPCASE(SUBSTR(PARAMETER,8,1));
ALPHA=ESTIMATE;
KEEP CL ALPHA;
RUN;

PROC APPEND BASE=_ALPHA DATA=_ALPHAB FORCE; RUN;

PROC SQL NOPRINT;
SELECT SUM(EXP(ALPHA)) INTO :_SUMEXP TRIMMED FROM _ALPHA;
QUIT;

DATA _MIXP;
SET _ALPHA;
PIE=EXP(ALPHA)/&_SUMEXP.;
KEEP CL PIE;
RUN;

PROC SORT DATA=_MIXP; BY CL; RUN;
PROC SORT DATA=_BW; BY CL; RUN;
PROC SORT DATA=_THW; BY CL; RUN;

DATA _COEF;
MERGE _BW _THW _MIXP;
BY CL;
TH_1=0;*first threshold fixed at 0 in the released macro;
CLASSNUM=(CL='A')*1+(CL='B')*2+(CL='C')*3;
RUN;

*Standard error check. Are all SEs present and positive;
%LET _HASSE=0;
PROC SQL NOPRINT;
SELECT COUNT(*) INTO :_HASSE TRIMMED
FROM DICTIONARY.COLUMNS
WHERE LIBNAME='WORK' AND UPCASE(MEMNAME)=UPCASE("OA_EST_K&KTRUE.")
  AND UPCASE(NAME)='STANDARDERROR';
QUIT;
%LET NBADSE=.;
%IF &_HASSE. %THEN %DO;
PROC SQL NOPRINT;
SELECT SUM(MISSING(STANDARDERROR) OR STANDARDERROR<=0) INTO :NBADSE TRIMMED
FROM WORK.OA_EST_K&KTRUE.;
QUIT;
%END;

*Metric 1. Shape recovery;
DATA _FIT_TRAJ;
SET _COEF;
ARRAY B B_0 B_1 B_2;
DO T=1 TO &T.;
ETA=0;
DO J=1 TO DIM(B);
IF NOT MISSING(B[J]) THEN ETA+B[J]*(T**(J-1));
END;
P0=PROBNORM(TH_1-ETA);
P1=PROBNORM(TH_2-ETA)-PROBNORM(TH_1-ETA);
P2=PROBNORM(TH_3-ETA)-PROBNORM(TH_2-ETA);
P3=1-PROBNORM(TH_3-ETA);
EY_FIT=0*P0+1*P1+2*P2+3*P3;
OUTPUT;
END;
KEEP CLASSNUM T EY_FIT;
RUN;

PROC SQL;
CREATE TABLE _PAIR AS
SELECT F.CLASSNUM AS FCLASS, TR.TRUECLASS, SUM((F.EY_FIT-TR.EY_TRUE)**2) AS SSE
FROM _FIT_TRAJ F CROSS JOIN TRUE_TRAJ TR
WHERE F.T=TR.T GROUP BY F.CLASSNUM, TR.TRUECLASS;
CREATE TABLE _MATCH AS
SELECT FCLASS, TRUECLASS AS MATCHED_TRUE, SSE
FROM _PAIR P
WHERE SSE=(SELECT MIN(SSE) FROM _PAIR P2 WHERE P2.FCLASS=P.FCLASS);
QUIT;
PROC SORT DATA=_MATCH NODUPKEY; BY FCLASS; RUN;

PROC SQL NOPRINT; SELECT COUNT(DISTINCT MATCHED_TRUE) INTO :NDIST TRIMMED FROM _MATCH; QUIT;
PROC SQL;
CREATE TABLE _SHAPE AS
SELECT &R. AS REP, (&NDIST.=&KTRUE.) AS BIJECTIVE, MAX(SQRT(SSE/&T.)) AS MAX_RMSE
FROM _MATCH;
QUIT;

*Metric 2. Classification by maximum posterior membership;
DATA _ASSIGN;
IF _N_=1 THEN DO;
IF 0 THEN SET _COEF;
DECLARE HASH H(DATASET:'_COEF');
_RC=H.DEFINEKEY('CLASSNUM');
_RC=H.DEFINEDATA('B_0','B_1','B_2','TH_1','TH_2','TH_3','PIE');
_RC=H.DEFINEDONE();
END;
SET BASE_FILE_SRS;
ARRAY Y Y1_1-Y1_&T.;
ARRAY BB[3,0:2] _TEMPORARY_;
ARRAY TT[3,3] _TEMPORARY_;
ARRAY PW[3] _TEMPORARY_;
DO CC=1 TO 3;
CALL MISSING(B_0,B_1,B_2,TH_1,TH_2,TH_3,PIE);
_RC=H.FIND(KEY:CC);
BB[CC,0]=COALESCE(B_0,0);
BB[CC,1]=COALESCE(B_1,0);
BB[CC,2]=COALESCE(B_2,0);
TT[CC,1]=TH_1;
TT[CC,2]=TH_2;
TT[CC,3]=TH_3;
PW[CC]=PIE;
END;
BESTLL=.;
ASSIGNED=.;
DO CC=1 TO 3;
LL=0;
DO T=1 TO &T.;
ETA=BB[CC,0]+BB[CC,1]*T+BB[CC,2]*T*T;
IF Y[T]=0 THEN PR=PROBNORM(TT[CC,1]-ETA);
ELSE IF Y[T]=1 THEN PR=PROBNORM(TT[CC,2]-ETA)-PROBNORM(TT[CC,1]-ETA);
ELSE IF Y[T]=2 THEN PR=PROBNORM(TT[CC,3]-ETA)-PROBNORM(TT[CC,2]-ETA);
ELSE PR=1-PROBNORM(TT[CC,3]-ETA);
LL + LOG(MAX(PR,1E-12));
END;
IF NOT MISSING(PW[CC]) THEN LL + LOG(MAX(PW[CC],1E-12));*class prior;
IF BESTLL=. OR LL>BESTLL THEN DO;
BESTLL=LL;
ASSIGNED=CC;
END;
END;
KEEP &IDVAR. TRUECLASS ASSIGNED;
RUN;

PROC SQL;
CREATE TABLE _ASG2 AS
SELECT A.*, M.MATCHED_TRUE AS ASSIGNED_TRUE
FROM _ASSIGN A LEFT JOIN _MATCH M ON A.ASSIGNED=M.FCLASS;
CREATE TABLE _ACC AS
SELECT &R. AS REP, MEAN(ASSIGNED_TRUE=TRUECLASS) AS AGREE FROM _ASG2;
QUIT;

*Validity filter. Ordered thresholds and finite betas in every class;
PROC SQL NOPRINT;
SELECT SUM( (TH_1<TH_2 AND TH_2<TH_3)
            AND (B_0 IS NOT NULL AND B_1 IS NOT NULL AND B_2 IS NOT NULL) )
INTO :NOK TRIMMED FROM _COEF;
QUIT;

DATA _ROW;
MERGE _SHAPE _ACC;
BY REP;
VALID = (&NOK.=&KTRUE.) AND BIJECTIVE;
%IF &_HASSE. %THEN %DO; CLEAN_HESS=(&NBADSE.=0); %END;
%ELSE %DO; CLEAN_HESS=.; %END;
RUN;
PROC APPEND BASE=MC_METRICS DATA=_ROW FORCE; RUN;
%END;

*Metric 3. BIC model selection over C=KMIN to C=KMAX;
%DO KK=&KMIN. %TO &KMAX.;
%ORDPROB_MIX_FIT_ONE(
DATA=BASE_FILE_SRS, ID=&IDVAR.,
YVARS=Y1_1-Y1_&T., TVARS=QUAR1-QUAR&T., TTOTAL=&T.,
M=&MCAT., YCODES=0 1 2 3, DEG=&DEG., K=&KK.,
TECH=&TECH., MAXITER=&MAXITER., PREFIX=OAK&KK.
);
%IF %SYSFUNC(EXIST(WORK.OAK&KK._FIT_K&KK.)) %THEN %DO;
DATA _BICROW;
SET WORK.OAK&KK._FIT_K&KK.;
WHERE UPCASE(DESCR) LIKE '%BIC%';
REP=&R.;
KFIT=&KK.;
BICVAL=VALUE;
KEEP REP KFIT BICVAL;
RUN;
PROC APPEND BASE=MC_BICPICK DATA=_BICROW FORCE; RUN;
%END;
%END;

*DELETE INTERMEDIATE FILES FOR THIS REPLICATION;
PROC DATASETS LIB=WORK NOLIST NOWARN MEMTYPE=DATA;
DELETE OA_EST_K&KTRUE. OA_ESTS_K&KTRUE. OA_FIT_K&KTRUE.
_B _BW _TH _THW _ALPHA _ALPHAB _MIXP _COEF _FIT_TRAJ _PAIR _MATCH _SHAPE
_ASSIGN _ASG2 _ACC _ROW _BICROW
%DO KK=&KMIN. %TO &KMAX.; OAK&KK._EST_K&KK. OAK&KK._ESTS_K&KK. OAK&KK._FIT_K&KK. %END; ;
QUIT;

%IF %SYSFUNC(MOD(&R.,25))=0 %THEN %PUT NOTE: ===== completed replication &R. of &NREP. =====;
%END;
%MEND RUN_VALIDATION;
%RUN_VALIDATION;

OPTIONS NOTES;
ODS EXCLUDE NONE;

/*##########################################################################################################################
*STEP 4: CHECK OUTPUT
##########################################################################################################################
The BIC-selected class count per replication, then the recovery summary across replications that passed
the validity filter. The same two tables are written to the PDF in OUTDIR.
##########################################################################################################################*/
PROC SORT DATA=MC_BICPICK; BY REP BICVAL; RUN;

DATA _BESTK;
SET MC_BICPICK;
BY REP;
IF FIRST.REP;
RUN;

PROC SQL NOPRINT;
SELECT COUNT(DISTINCT REP) INTO :NREP_SEEN TRIMMED FROM MC_BICPICK;
SELECT SUM(VALID=0) INTO :NDROP TRIMMED FROM MC_METRICS;
SELECT SUM(CLEAN_HESS=1), COUNT(CLEAN_HESS) INTO :NCLEAN TRIMMED, :NSE TRIMMED FROM MC_METRICS;
QUIT;

TITLE "TRAJ2 ORDINAL VALIDATION &NREP. REPLICATIONS AT N=&NSIM.";

PROC FREQ DATA=_BESTK;
TABLES KFIT / NOCUM;
TITLE2 "Class count chosen by minimum BIC. Target C=&KTRUE. over &NREP_SEEN. replications";
RUN;

PROC MEANS DATA=MC_METRICS(WHERE=(VALID=1)) N MEAN STD MIN MAX MAXDEC=4;
VAR MAX_RMSE AGREE;
TITLE2 "Recovery summary. Fully estimable standard errors in &NCLEAN. of &NSE. fits";
RUN;
TITLE;

%PUT NOTE: Dropped &NDROP. rep(s) on the validity filter. Fully estimable SEs: &NCLEAN. / &NSE.;

ODS PDF FILE="&OUTDIR./ordprob_validation_final.pdf" STARTPAGE=NO;
TITLE "TRAJ2 ORDINAL VALIDATION &NREP. REPLICATIONS AT N=&NSIM.";
PROC MEANS DATA=MC_METRICS(WHERE=(VALID=1)) N MEAN STD MIN MAX MAXDEC=4;
VAR MAX_RMSE AGREE;
RUN;
PROC FREQ DATA=_BESTK;
TABLES KFIT / NOCUM;
RUN;
ODS PDF CLOSE;
TITLE;

*DELETE INTERMEDIATE FILES;
PROC DATASETS LIB=WORK NOLIST NOWARN MEMTYPE=DATA;
DELETE BASE_FILE_SRS TRUE_TRAJ;
QUIT;

/*##########################################################################################################################
*END
##########################################################################################################################
OUTPUT FILES
MC_METRICS   one row per replication. Worst-class RMSE, classification agreement, validity and SE flags
MC_BICPICK   one row per replication and fitted class count. BIC value
_BESTK       one row per replication. The class count with the minimum BIC
ordprob_validation_final.pdf   the two summary tables, written to OUTDIR
##########################################################################################################################*/
