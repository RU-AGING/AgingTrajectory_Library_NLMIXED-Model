*************************************************************************************************************************
*************************************************************************************************************************
COMMUNITY HEALTH AND AGING OUTCOMES (CHAO) LAB - INSTITUTE FOR HEALTH, HEALTH CARE POLICY & AGING RESEARCH - 
RUTGERS, THE STATE UNIVERSITY OF NEW JERSEY                  
*************************************************************************************************************************
*************************************************************************************************************************
* PROJECT NAME: Group-Based Trajectory Modeling using NLMIX                                                                            
* LAST UPDATED DATE: 19 Feb 2025                                                                                
*************************************************************************************************************************
* DATA SOURCES: Simulation or real-world longitudinal datasets containing repeated measures for multiple time points 
*************************************************************************************************************************
* PURPOSE:                                                                                                              
* This code implements group-based trajectory modeling with joint trajectories for two outcomes (T1 and T2) using the          
* SAS NLMIXED procedure. It includes single-outcome models, a multi-trajectory model, and a plotting macro for visualization.     		                                     
*************************************************************************************************************************
* Cohort: Individuals with repeated measures data for two outcomes over N time periods.     
*************************************************************************************************************************
* AUTHORS:                                                                                                             
* Weiyi Xia, Haiqun Lin                                                                                   
*************************************************************************************************************************
* SYSTEM RESOURCES & PERFORMANCE:                                                                                       
* - Execution Environment: CMS VIRTUAL RESEARCH DATA CENTER (VRDC), running SAS 8.1 or later versions                   
* - Expected runtime in Standard Analytic Container: 45 minutes using data for 200,000 beneficiaries       
* - Expected runtime in XL Analytic Container: 1.5 hours using data for 200,000 beneficiaries       		
* - Storage Space For Construction of the Datasets: 20 GB
* - File Size of the Final Analytic File: 1 GB
*************************************************************************************************************************;                
* OUTPUT FILE DESCRIPTION:                                                                                               
*   1. Single-outcome trajectory model results: nlm_fix_T1&class. and nlm_fix_T2&class.                                                 
*   2. Multitrajectory model results: nlm_fix_T1_T2&class.                                                                              
*   3. Visualizations of predicted vs. observed trajectories and averaged posterior class memberships.    
*************************************************************************************************************************
* Ref: Jones, B. L., & Nagin, D. S. (2007). Advances in Group-Based Trajectory Modeling and an SAS Procedure for Estimating Them.       
* Sociological Methods & Research, 35(4), 542–571. https://doi.org/10.1177/0049124106292364      
*************************************************************************************************************************;

/*==============================Step 0: Load macros================================*/
%let class_names = A B C D E F G H I J K; /* Define all latent class labels. */
%let max_values = 92 91 91 92 92 91 91 92 91 91 91 91;

/* Macro to initialize starting values for alpha parameters across classes */
%macro starting_value_alpha(class);   
%local s;  
* Loop through latent classes starting from the second (class 2 to class N); 
%do s=2 %to &class.;
%let class_ = %scan(&class_names, &s); /* Extract the class name. */
alpha0_&class_.=0 /* Set the initial value for alpha for each class. */
%end;
%mend starting_value_alpha;


/* Macro to initialize starting values for beta and sigma parameters */
%macro starting_value_beta_sigma(class,outcome,order,equal_sigma);   
%local s;  
/* Initialize sigma parameters. If equal variances, set one global sigma; otherwise, per class. */
%if &equal_sigma. eq T %then %do; sigma&outcome._=30 %let class2_=%str(); %end;
%else %do;
sigma&outcome._A=30 /* Assign sigma for the first class. */
%do s=2 %to &class.; /* Loop for remaining classes. */
%let class_=%scan("&class_all.", &s, ", ");
%let class2_=&class_.;
sigma&outcome._&class2_.=30
%end;
%end;
/* Initialize beta coefficients for the polynomial terms. */
%do s=1 %to &class.;
%let class_=%scan("&class_all.", &s, ", ");
%do i=2 %to &order.;
beta&outcome._&class_.&i.=0
%end;
%end;
%mend starting_value_beta_sigma;

/* Macro to set bounds for alpha parameters */
%macro bounds_alpha(bounds_alpha,class);   
%local s;
-&bounds_alpha.<alpha0_B<&bounds_alpha. /* Set bounds for the second class. */
%do s=3 %to (&class.);
%let class_=%scan("&class_all.", &s, ", ");
,-&bounds_alpha.<alpha0_&class_.<&bounds_alpha. /* Set bounds for remaining classes. */
%end;
%mend bounds_alpha;

/* Macro to set bounds for sigma parameters */
%macro bounds_sigma(bounds_sigma,class,outcome,equal_sigma);   
%local s; %local class2_;
%if &equal_sigma. eq T %then %do; sigma&outcome._>&bounds_sigma.; %end;
%else %do;
sigma&outcome._A>&bounds_sigma. /* Set bounds for the first class. */
%do s=2 %to &class.; /* Loop for remaining classes. */
%let class_=%scan("&class_all.", &s, ", ");
%let class2_=&class_.;
,sigma&outcome._&class2_.>&bounds_sigma.
%end;
%end;
%mend bounds_sigma;

/* Macro to initialize universal arrays */
%macro initiation_universal(T);   
%str(ARRAY) X[&T.] quar1-quar&T.%str(;)
%str(ARRAY) Y1[&T.] /*OUTCOME 1 Variable*/%str(;)
%str(ARRAY) Y2[&T.] /*OUTCOME 2 Variable*/&T.%str(;)
%mend initiation_universal;

/* Macro to initialize arrays specific to classes and outcomes */
%macro initiation(T,class,outcome);   
%local s; 
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%str(ARRAY) PI&class_.&outcome.[&T.] PI&class_.&outcome._1-PI&class_.&outcome._&T.%str(;)
%str(ARRAY) mu&class_.&outcome.[&T.] mu&class_.&outcome._1-mu&class_.&outcome._&T.%str(;)
%end;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
PROD&class_.&outcome.=0%str(;)
%end;
%mend initiation;

/* Macro for constructing the model equation for each class and outcome */
%macro model(pattern,order,class,outcome);   
%local i; 
%local s; 
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
mu&class_.&outcome.[I]=beta&outcome._&class_.0
%do i=1 %to &order.;
+ %sysfunc(tranwrd (%sysfunc(tranwrd(%sysfunc(tranwrd(&pattern,order,&i.)),outcome,&outcome.)),class,&class_.))
%end;%str(;)
%end;
%mend model;

/* Macro to calculate residuals */
%macro residual(class,outcome);   
%local i; 
%local s; 
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
e&outcome._&class_.=Y&outcome.[I]-mu&class_.&outcome.[I] %str(;)
%end;
%mend residual;

/* Macro to apply floating control on residuals */
%macro float_control(float,class,outcome,equal_sigma);   
%local s; 
%local class2_;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%if &equal_sigma. eq T %then %do; %let class2_=%str();%end; %else %do; %let class2_=&class_.; %end;
e&outcome._&class_. = min(max(e&outcome._&class_., -&float.*sigma&outcome._&class2_.), &float.*sigma&outcome._&class2_.);
%end;
%mend float_control;

/* Macro for calculating log-probabilities for normal distributions */
%macro Prob_cnorm(class,outcome,equal_sigma);   
%local s; 
%local class2_;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%if &equal_sigma. eq T %then %do; %let class2_=%str();%end; %else %do; %let class2_=&class_.; %end;
PI&class_.&outcome.[I] = (logpdf('NORMAL',e&outcome._&class_./sigma&outcome._&class2_.))-log(sigma&outcome._&class2_.)%str(;)
if Y&outcome.[I] eq 0  then do%str(;)
PI&class_.&outcome.[I] = logcdf('NORMAL',e&outcome._&class_./sigma&outcome._&class2_.)%str(;)
end%str(;)
   else if Y&outcome.[I] eq MAX[I] then do%str(;)
PI&class_.&outcome.[I] = logcdf('NORMAL',(-e&outcome._&class_.)/sigma&outcome._&class2_.)%str(;)
end%str(;)
PROD&class_.&outcome.=PROD&class_.&outcome.+PI&class_.&outcome.[I]%str(;)
%end;
%mend Prob_cnorm;

/* Macro for calculating log-probabilities for Poisson distributions */
%macro Prob_pois(class,outcome);   
%local s; 
%do s=1 %to &class;
    %let class_=%scan("&class_all.", &s, ", ");
    /* Poisson log PMF per time point */
    PI&class_.&outcome.[I] = logpdf('POISSON', Y&outcome.[I], mu&class_.&outcome.[I]);
    /* accumulate across time */
    PROD&class_.&outcome.=PROD&class_.&outcome.+PI&class_.&outcome.[I];
%end;
%mend Prob_pois;

/* Macro for calculating log-probabilities for Zero-Inflated Poisson */
%macro Prob_zip(class,outcome);   
%local s; 
%do s=1 %to &class;
    %let class_=%scan("&class_all.", &s, ", ");
    
    /* Zero inflation probability for this class */
    p&class_.&outcome. = logistic(gamma&outcome._&class_.);

    if Y&outcome.[I]=0 then do;
        PI&class_.&outcome.[I] = log(
            p&class_.&outcome. + (1 - p&class_.&outcome.)*pdf('POISSON',0,mu&class_.&outcome.[I])
        );
    end;
    else do;
        PI&class_.&outcome.[I] = log(
            (1 - p&class_.&outcome.) * pdf('POISSON',Y&outcome.[I], mu&class_.&outcome.[I])
        );
    end;

    PROD&class_.&outcome.=PROD&class_.&outcome.+PI&class_.&outcome.[I];
%end;
%mend Prob_zip;



/* Macro for calculating class membership*/
%macro Class_Membership(class);   
%local s; 
alpha0_A=0%str(;)
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
pinumer_&class_. = exp(alpha0_&class_.)%str(;)
%end;pideno =1
%do s=2 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
+pinumer_&class_.
%end;%str(;)
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
pie_&class_. = pinumer_&class_./pideno%str(;)
%end;
%mend Class_Membership;

/* Macro for calculating log likelihood */
%macro LogLike(class,outcome);   
%local s; l_latclass = pie_A*exp(PRODA&outcome.)
%do s=2 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
+pie_&class_.*exp(PROD&class_.&outcome.)
%end;%str(;)
%mend LogLike;

%macro LogLike_multi(class);   
%local s;(pie_A*exp(PRODA1)*exp(PRODA2)
%do s=2 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
+pie_&class_.*exp(PROD&class_.1)*exp(PROD&class_.2)
%end;)%str(;)
%mend LogLike_multi;

%macro Posterior(class); 
%local s;

%do s=1 %to &class.;
%let class_=%scan("&class_all.", &s., ", ");
post_&class_.=pie_&class_.*exp(PROD&class_.1)*exp(PROD&class_.2)/%LogLike_multi( class=&class.)
%end;%str(;)
keep  %do s=1 %to &class.;
%let class_=%scan("&class_all.", &s., ", ");
post_&class_.
%end;%str(;)
%mend Posterior;


%macro initiation_pred(T,class,outcome);   
%local s; 
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%str(ARRAY) Pred&class_.&outcome.[&T.] Pred&class_.&outcome._1-Pred&class_.&outcome._&T.%str(;)
%str(ARRAY) mu&class_.&outcome.[&T.] mu&class_.&outcome._1-mu&class_.&outcome._&T.%str(;)
%end;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
PROD&class_.&outcome.=0%str(;)
%end;
%mend initiation_pred;


%macro Pred_cnorm(class,outcome,equal_sigma);   
%local s; 

%local class2_;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%if &equal_sigma. eq T %then %do; %let class2_=%str();%end; %else %do; %let class2_=&class_.; %end;
temp=(exp(logcdf('normal',(MAX[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.))-
            exp(logcdf('normal',(0-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)))%str(;)
if temp = 0 then temp=0.000001 %str(;)

Pred&class_.&outcome.[I] =0+(exp(logcdf('NORMAL',(MAX[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.))-
exp(logcdf('NORMAL',(0-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)))*
(mu&class_.&outcome.[I]+sigma&outcome._&class2_.*(exp(logpdf('normal',-mu&class_.&outcome.[I]/sigma&outcome._&class2_.))-
            exp(logpdf('normal',(MAX[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)))/temp)+
MAX[I]*exp(logcdf('normal',(-MAX[i]+mu&class_.&outcome.[I])/sigma&outcome._&class2_.))

%str(;)

%end;%str(;)
%mend Pred_cnorm;

%macro Pred_pois(class,outcome);   
%local s; 
%do s=1 %to &class;
    %let class_=%scan("&class_all.", &s, ", ");
    /* Predicted mean is just µ for Poisson */
    Pred&class_.&outcome.[I] = mu&class_.&outcome.[I];
%end;
%mend Pred_pois;

%macro Pred_zip(class,outcome);   
%local s; 
%do s=1 %to &class;
    %let class_=%scan("&class_all.", &s, ", ");
    p&class_.&outcome. = logistic(gamma&outcome._&class_.);
    Pred&class_.&outcome.[I] = (1 - p&class_.&outcome.)*mu&class_.&outcome.[I];
%end;
%mend Pred_zip;


%macro plot_prep(T,LC,result,order,equal_sigma,dist);


data parameter;set &result. ;keep  parameter estimate;run;

proc transpose data=parameter out=parameter;id parameter;run;

proc sql;
create table data_pred as
select * from base_file_srs, parameter;quit;

data pred_membership_y;set data_pred;
/*Create arrays*/
 %initiation_universal(T=&T.);
ARRAY MAX[&T.] (&max_values);
%initiation(T=&T., class=&LC.,outcome=1);
%initiation(T=&T., class=&LC.,outcome=2);

DO I=1 TO &T.;
/*Model prediction*/
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=1);
%residual( class=&LC.,outcome=1)
%float_control( float=8 ,class=&LC.,outcome=1,equal_sigma=&equal_sigma.);
%Prob_dispatch(dist=&dist., class=&LC., outcome=1, equal_sigma=&equal_sigma.);
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=2);
%residual( class=&LC.,outcome=2)
%float_control( float=8 ,class=&LC.,outcome=2,equal_sigma=&equal_sigma.);
%Prob_dispatch(dist=&dist., class=&LC., outcome=2, equal_sigma=&equal_sigma.);
end;

/*Class Memebership P(X=t)*/
 %Class_Membership( class=&LC.);

/*Likelihood logLike=sum P(X=t)*P(Y|X=t)*/
 %Posterior(class=&LC.);
run; 

data y_1;
set base_file_srs;
run;

proc iml;
Start colsum(m);
return (m[+,]);
finish;
use y_1; read all var _num_ into y;
use pred_membership_y; read all var _num_ into pred_membership_y;
sum=colsum(pred_membership_y);*print sum;
invsum=1/sum; *print invsum;
avg_y=y`*pred_membership_y;*print avg_y;
avg_y_2=avg_y#invsum; *print avg_y_2;
create avg_y1 from avg_y_2;
append from avg_y_2;
close avg_y_2;
quit;data y_2;
set base_file_srs (keep= QINP1 QINP2 QINP3 QINP4 QINP5 QINP6 QINP7 QINP8 QINP9 QINP10 QINP11 QINP12);
run;

proc iml;
Start colsum(m);
return (m[+,]);
finish;
use y_2; read all var _num_ into y;
use pred_membership_y; read all var _num_ into pred_membership_y;
sum=colsum(pred_membership_y);*print sum;
invsum=1/sum;* print invsum;
avg_y=y`*pred_membership_y;*print avg_y;
avg_y_2=avg_y#invsum; *print avg_y_2;
create avg_y2 from avg_y_2;
append from avg_y_2;
close avg_y_2;
quit;

data long;set &result. ;keep  parameter estimate;run;

proc transpose data=long out=wide;id parameter;run;
data wide; set wide;quar1=1;quar2=2;quar3=3;quar4=4;quar5=5;quar6=6;quar7=7;quar8=8;quar9=9;quar10=10;quar11=11;quar12=12;
run;

data wide; set wide;
/*Create arrays*/
ARRAY X[12] quar1-quar12;
ARRAY MAX[&T.] ((&max_values));
 %initiation_pred(T=&T., class=&LC.,outcome=1);
 %initiation_pred(T=&T., class=&LC.,outcome=2);;

DO I=1 TO 12 ;
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=1);
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=2);

%Pred_dispatch(dist=&dist., class=&LC., outcome=1, equal_sigma=&equal_sigma.);
%Pred_dispatch(dist=&dist., class=&LC., outcome=2, equal_sigma=&equal_sigma.);

end;
run;

proc transpose data=wide out=pred_temp;run;

data pred;
set pred_temp;
format outcome $char3.;
format type $char3.;
format quar 6.;
where _NAME_ contains 'Pred';
class=substr(_NAME_,5,1);
y=substr(_NAME_,6,1);
if y eq'1' then outcome = 'HH';
if y eq'2' then outcome = 'INP';
type='pred';
quar=substr(_NAME_,8,2);
drop _name_;
run;

title 'Predicated Trajectory by Latent Class';
PROC SGPANEL data=pred ;
    *styleattrs datacontrastcolors=(red green blue);
PANELBY  class/columns=4;
    series x=quar y=estimate / group=outcome;
run;
title;

data avg_y1_plot;
set avg_y1 ;
format outcome $char3.;
format type $char3.;
array col col1-col&LC.;
quar=_n_;
do i = 1 to &LC.;
Estimate=col[i];
class=byte(i+96);
outcome="HH";
type="Avg";
output;
end;

keep estimate class quar outcome type; run;

data avg_y2_plot;
set avg_y2 ;
format outcome $char3.;
format type $char3.;
array col col1-col&LC.;
quar=_n_;
do i = 1 to &LC.;
Estimate=col[i];
class=byte(i+96);
outcome="INP";
type="Avg";
output;
end;
keep estimate class quar outcome type;
run;

data avg;
format quar 6.;
set avg_y1_plot avg_y2_plot;
class=PROPCASE(class);
run;

proc sort data=avg out=test; by outcome quar class ;run;

proc sort data=pred; by outcome quar class ;run;

data data_plot;
merge avg(rename=(Estimate= Avg)) pred(rename=(Estimate= pred));
by outcome quar class ;
run;

title 'Predicted vs Averaged Days by Latent Class ';
proc sgpanel data=data_plot;
panelby outcome;
series x=quar y=pred/group=class name="pred";
scatter x=quar y=avg/group=class name="obs";
keylegend "pred" /title="Predicted";
keylegend "obs" /title="Averaged Observed";
run;
quit;
title ;

/*table the averaged posterior prob*/
title 'Averaged Posterior Class Membership by Latent Class ';
proc means data= pred_membership_y mean missing;
var ;
output out=avg_membership mean= Avg  std=SD;
run;
title ;


%mend plot_prep;


%macro nlmixed_1_cnorm(T,LC,Y,starting,output,order,equal_sigma);

proc nlmixed data= base_file_srs  itdetails  qpoints=40 noad maxiter=1000 tech=dbldog ;
 bounds 	%bounds_alpha(bounds_alpha=3,class=&LC.),
		%bounds_sigma(bounds_sigma=0,class=&LC.,outcome=&Y.,equal_sigma=&equal_sigma.);
parms 	&starting.

/*Create arrays*/
%initiation_universal(T=&T.);;
ARRAY MAX[&T.] (&max_values);

%initiation(T=&T., class=&LC.,outcome=&Y.);;

/*Fit model P(Y=y|X=t)*/
DO I=1 TO &T.;
/*Model building & Floating Control*/
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=&Y.);
%residual( class=&LC.,outcome=&Y.)
%float_control( float=8 ,class=&LC.,outcome=&Y.,equal_sigma=&equal_sigma.);
%Prob_cnorm( class=&LC.,outcome=&Y.,equal_sigma=&equal_sigma.);;
end;

/*Class Memebership P(X=t)*/
 %Class_Membership( class=&LC.);

/*Likelihood logLike=sum P(X=t)*P(Y|X=t)*/
%LogLike( class=&LC.,outcome=&Y.);;

ll_latclass=log(l_latclass);
model ll_latclass ~ general(ll_latclass);
ods output ParameterEstimates=work.&output.; 
run; 

%mend nlmixed_1_cnorm;


%macro nlmixed_MultiTraj_cnorm(T,LC,starting,output,order,equal_sigma);

proc nlmixed data= base_file_srs  itdetails  qpoints=40 noad maxiter=1000 tech=dbldog ;
 bounds  %bounds_alpha(bounds_alpha=3,class=&LC.) ,
	 %bounds_sigma(bounds_sigma=0,class=&LC.,outcome=1,equal_sigma=&equal_sigma.),
		%bounds_sigma(bounds_sigma=0,class=&LC.,outcome=2,equal_sigma=&equal_sigma.);

parms 	&starting./*use the previous result as starting point*/
		;

/*Create arrays*/
%initiation_universal(T=&T.);;
ARRAY MAX[&T.] (&max_values);
%initiation(T=&T., class=&LC.,outcome=1);;
%initiation(T=&T., class=&LC.,outcome=2);;

/*Fit model P(Y=y|X=t)*/
DO I=1 TO &T.;
/*Model building & Floating Control*/
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=1);
%residual( class=&LC.,outcome=1)
%float_control( float=8 ,class=&LC.,outcome=1,equal_sigma=&equal_sigma.);
%Prob_cnorm( class=&LC.,outcome=1,equal_sigma=&equal_sigma.);;

%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=2);
%residual( class=&LC.,outcome=2)
%float_control( float=8 ,class=&LC.,outcome=2,equal_sigma=&equal_sigma.);
%Prob_cnorm( class=&LC.,outcome=2,equal_sigma=&equal_sigma.);;
end;

/*Class Memebership P(X=t)*/
 %Class_Membership(class=&LC.);

/*Likelihood logLike=sum P(X=t)*P(Y|X=t)*/
 l_latclass = %LogLike_multi(class=&LC.);;

ll_latclass=log(l_latclass);
model ll_latclass ~ general(ll_latclass);
ods output ParameterEstimates=work.&output.; 
run; 
%mend nlmixed_MultiTraj_cnorm;



%macro nlmixed_1_pois(T,LC,Y,starting,output,order);

proc nlmixed data= base_file_srs itdetails qpoints=40 noad maxiter=1000 tech=dbldog;
    bounds %bounds_alpha(bounds_alpha=3,class=&LC.);
    parms &starting.;

    /* Create arrays */
    %initiation_universal(T=&T.);
    %initiation(T=&T., class=&LC.,outcome=&Y.);

    /* Fit trajectory model */
    DO I=1 TO &T.;
        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=&Y.);
        %Prob_pois(class=&LC.,outcome=&Y.);
    END;

    /* Class membership */
    %Class_Membership(class=&LC.);

    /* Likelihood */
    %LogLike(class=&LC.,outcome=&Y.);

    ll_latclass = log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output.; 
run; 

%mend nlmixed_1_pois;

%macro nlmixed_MultiTraj_pois(T,LC,starting,output,order);

proc nlmixed data= base_file_srs itdetails qpoints=40 noad maxiter=1000 tech=dbldog;
    bounds %bounds_alpha(bounds_alpha=3,class=&LC.);
    parms &starting.;

    /* Create arrays */
    %initiation_universal(T=&T.);
    %initiation(T=&T., class=&LC.,outcome=1);
    %initiation(T=&T., class=&LC.,outcome=2);

    DO I=1 TO &T.;
        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=1);
        %Prob_pois(class=&LC.,outcome=1);

        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=2);
        %Prob_pois(class=&LC.,outcome=2);
    END;

    /* Class membership */
    %Class_Membership(class=&LC.);

    /* Multi-trajectory likelihood */
    l_latclass = %LogLike_multi(class=&LC.);

    ll_latclass=log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output.; 
run; 

%mend nlmixed_MultiTraj_pois;

%macro nlmixed_1_zip(T,LC,Y,starting,output,order);

proc nlmixed data= base_file_srs itdetails qpoints=40 noad maxiter=1000 tech=dbldog;
    bounds %bounds_alpha(bounds_alpha=3,class=&LC.);
    parms &starting.;
    
    /* Create arrays */
    %initiation_universal(T=&T.);
    %initiation(T=&T., class=&LC.,outcome=&Y.);
    
    DO I=1 TO &T.;
        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=&Y.);
        %Prob_zip(class=&LC.,outcome=&Y.);
    END;

    %Class_Membership(class=&LC.);
    %LogLike(class=&LC.,outcome=&Y.);

    ll_latclass = log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output.; 
run; 

%mend nlmixed_1_zip;
%macro nlmixed_MultiTraj_zip(T,LC,starting,output,order);

proc nlmixed data= base_file_srs itdetails qpoints=40 noad maxiter=1000 tech=dbldog;
    bounds %bounds_alpha(bounds_alpha=3,class=&LC.);
    parms &starting.;
    
    %initiation_universal(T=&T.);
    %initiation(T=&T., class=&LC.,outcome=1);
    %initiation(T=&T., class=&LC.,outcome=2);

    DO I=1 TO &T.;
        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=1);
        %Prob_zip(class=&LC.,outcome=1);

        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=2);
        %Prob_zip(class=&LC.,outcome=2);
    END;

    %Class_Membership(class=&LC.);
    l_latclass = %LogLike_multi(class=&LC.);

    ll_latclass=log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output.; 
run; 

%mend nlmixed_MultiTraj_zip;

/* Choose probability generator based on DIST */
%macro Prob_dispatch(dist,class,outcome,equal_sigma);
    %if &dist=normal %then %do;
        %Prob_cnorm(class=&class.,outcome=&outcome.,equal_sigma=&equal_sigma.);
    %end;
    %else %if &dist=poisson %then %do;
        %Prob_pois(class=&class.,outcome=&outcome.);
    %end;
    %else %if &dist=zip %then %do;
        %Prob_zip(class=&class.,outcome=&outcome.);
    %end;
%mend Prob_dispatch;

/* Choose prediction generator based on DIST */
%macro Pred_dispatch(dist,class,outcome,equal_sigma);
    %if &dist=normal %then %do;
        %Pred_cnorm(class=&class.,outcome=&outcome.,equal_sigma=&equal_sigma.);
    %end;
    %else %if &dist=poisson %then %do;
        %Pred_pois(class=&class.,outcome=&outcome.);
    %end;
    %else %if &dist=zip %then %do;
        %Pred_zip(class=&class.,outcome=&outcome.);
    %end;
%mend Pred_dispatch;

%macro nlmixed_1(T,LC,Y,starting,output,order,equal_sigma,dist);

proc nlmixed data= base_file_srs itdetails qpoints=40 noad maxiter=1000 tech=dbldog;
    bounds %bounds_alpha(bounds_alpha=3,class=&LC.);
    parms &starting.;

    %initiation_universal(T=&T.);
    %initiation(T=&T., class=&LC.,outcome=&Y.);

    DO I=1 TO &T.;
        %model(pattern=X[I]**order*betaoutcome_classorder, order=&order.,class=&LC.,outcome=&Y.);
        %Prob_dispatch(dist=&dist., class=&LC., outcome=&Y., equal_sigma=&equal_sigma.);
    END;

    %Class_Membership(class=&LC.);
    %LogLike(class=&LC.,outcome=&Y.);

    ll_latclass=log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output.; 
run;

%mend nlmixed_1;

%macro nlmixed_MultiTraj(T,LC,starting,output,order,equal_sigma,dist);

proc nlmixed data= base_file_srs itdetails qpoints=40 noad maxiter=1000 tech=dbldog;
    bounds %bounds_alpha(bounds_alpha=3,class=&LC.);
    parms &starting.;

    %initiation_universal(T=&T.);
    %initiation(T=&T., class=&LC.,outcome=1);
    %initiation(T=&T., class=&LC.,outcome=2);

    DO I=1 TO &T.;
        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=1);
        %Prob_dispatch(dist=&dist., class=&LC., outcome=1, equal_sigma=&equal_sigma.);

        %model(pattern=X[I]**order*betaoutcome_classorder , order=&order.,class=&LC.,outcome=2);
        %Prob_dispatch(dist=&dist., class=&LC., outcome=2, equal_sigma=&equal_sigma.);
    END;

    %Class_Membership(class=&LC.);
    l_latclass = %LogLike_multi(class=&LC.);

    ll_latclass=log(l_latclass);
    model ll_latclass ~ general(ll_latclass);

    ods output ParameterEstimates=work.&output.; 
run;

%mend nlmixed_MultiTraj;

