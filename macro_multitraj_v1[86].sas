/*This is the macro used to fit the group-based trajectory model with joint trajectories and dynamic upper bounds for truncated normal data */
/* Author: Weiyi Xia */

/* Test use
%Let T=12;
%Let LC=2;
%Let order=3;
%Let Y=1;
%Let equal_sigma=T;
PROC IMPORT OUT= WORK.base_file_srs 
            DATAFILE= "D:\Research\Lin\R33\2y\data_plot.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
*/


/*==============================Step 0: Load macros================================*/

%let class_all = A, B, C, D, E, F, G, H, I, J, K;

%macro starting_value_alpha(class);   
%local s;  
/*alpha*/
%do s=2 %to &class.;
%let class_=%scan("&class_all.", &s, ", ");
alpha0_&class_.=0 
%end;
%mend starting_value_alpha;
* %put %starting_value_alpha(class=3);;

%macro starting_value_beta_sigma(class,outcome,order,equal_sigma);   
%local s;  
/*sigma*/ 
%if &equal_sigma. eq T %then %do; sigma&outcome._=30
%let class2_=%str();%end;
%else %do;
sigma&outcome._A=30
%do s=2 %to &class.;
%let class_=%scan("&class_all.", &s, ", ");
%let class2_=&class_.; 
  sigma&outcome._&class2_.=30
%end;%end;
/*beta*/
%do s=1 %to &class.;
%let class_=%scan("&class_all.", &s, ", ");
%do i=2 %to &order.;
  beta&outcome._&class_.&i.=0
%end;
%end;

%mend starting_value_beta_sigma;
*  %put %starting_value_beta_sigma(class=3,outcome=2,order=3,equal_sigma=F);;
*  %put %starting_value_beta_sigma(class=3,outcome=2,order=3,equal_sigma=T);;


%macro bounds_alpha(bounds_alpha,class);   
%local s; -&bounds_alpha.<alpha0_B<&bounds_alpha.
%do s=3 %to (&class.);
%let class_=%scan("&class_all.", &s, ", ");
,-&bounds_alpha.<alpha0_&class_.<&bounds_alpha.
%end;

%mend bounds_alpha;
*  %put  %bounds_alpha(bounds_alpha=8,class=3);

%macro bounds_sigma(bounds_sigma,class,outcome,equal_sigma);   
%local s;
%local class2_;
%if &equal_sigma. eq T %then %do; sigma&outcome._>&bounds_sigma. %end; %else %do;
sigma&outcome._A>&bounds_sigma. 
%do s=2 %to &class.;
%let class_=%scan("&class_all.", &s, ", "); %let class2_=&class_.;
,sigma&outcome._&class2_.>&bounds_sigma. 
%end;
%end;
%mend bounds_sigma;
*  %put   %bounds_sigma(bounds_sigma=0,class=3,outcome=2,equal_sigma=T);
*  %put   %bounds_sigma(bounds_sigma=0,class=4,outcome=2,equal_sigma=F);

%macro initiation_universal(T);   
%str(ARRAY) X[&T.] quar1-quar&T.%str(;)
%str(ARRAY) Y1[&T.] Q1HH Q2HH Q3HH Q4HH Q5HH Q6HH Q7HH Q8HH Q9HH Q10HH Q11HH Q12HH%str(;)
%str(ARRAY) Y2[&T.] QINP1-QINP&T.%str(;)
%mend initiation_universal;
*  %put   %initiation_universal(T=12);;

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
*  %put   %initiation(T=12, class=4,outcome=2);;

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
*  %put   %model( pattern=X[I]**order*betaoutcome_classorder ,order=3,class=1,outcome=2);
*  %put   %model( pattern=X[I]**order*betaoutcome_classorder ,order=3,class=1,outcome=2);


%macro residual(class,outcome);   
%local i; 
%local s; 
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
e&outcome._&class_.=Y&outcome.[I]-mu&class_.&outcome.[I] %str(;)
%end;
%mend residual;
*  %put   %residual( class=1,outcome=2);
*  %put   %residual( class=1,outcome=2);

%macro float_control(float,class,outcome,equal_sigma);   
%local s; 
%local class2_;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%if &equal_sigma. eq T %then %do; %let class2_=%str();%end; %else %do; %let class2_=&class_.; %end;
if e&outcome._&class_. lt -&float.*sigma&outcome._&class2_. then e&outcome._&class_.=-&float.*sigma&outcome._&class2_.%str(;)
if e&outcome._&class_. gt &float.*sigma&outcome._&class2_. then e&outcome._&class_.=&float.*sigma&outcome._&class2_.%str(;)
%end;
%mend float_control;
*  %put   %float_control( float=8 ,class=1,outcome=2,equal_sigma=F);
*  %put   %float_control( float=8 ,class=1,outcome=2,equal_sigma=T);

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
*  %put   %Prob_cnorm(class=2,outcome=2,equal_sigma=F);;
*  %put   %Prob_cnorm(class=2,outcome=2,equal_sigma=T);;

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
*  %put %Class_Membership( class=4);;

%macro LogLike(class,outcome);   
%local s; l_latclass = pie_A*exp(PRODA&outcome.)
%do s=2 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
+pie_&class_.*exp(PROD&class_.&outcome.)
%end;%str(;)
%mend LogLike;
*  %put %LogLike( class=4,outcome=2);;

%macro LogLike_multi(class);   
%local s;(pie_A*exp(PRODA1)*exp(PRODA2)
%do s=2 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
+pie_&class_.*exp(PROD&class_.1)*exp(PROD&class_.2)
%end;)%str(;)
%mend LogLike_multi;
*  %put %LogLike_multi( class=4);;

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
*  %put %Posterior(class=4);;


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
*  %put   %initiation_pred(T=12, class=4,outcome=2);;

%macro Pred_cnorm(class,outcome,equal_sigma);   
%local s; 

%local class2_;
%do s=1 %to &class;
%let class_=%scan("&class_all.", &s, ", ");
%if &equal_sigma. eq T %then %do; %let class2_=%str();%end; %else %do; %let class2_=&class_.; %end;
temp=(exp(logcdf('normal',(MAX[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.))-
            exp(logcdf('normal',(0-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)))%str(;)
if temp = 0 then temp=0.000001 %str(;)
/*Pred&class_.&outcome.[I] = (logpdf('NORMAL',e&outcome._&class_./sigma&outcome._&class2_.))-log(sigma&outcome._&class2_.)%str(;)*/
Pred&class_.&outcome.[I] =0+(exp(logcdf('NORMAL',(MAX[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.))-
exp(logcdf('NORMAL',(0-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)))*
(mu&class_.&outcome.[I]+sigma&outcome._&class2_.*(exp(logpdf('normal',-mu&class_.&outcome.[I]/sigma&outcome._&class2_.))-
            exp(logpdf('normal',(MAX[I]-mu&class_.&outcome.[I])/sigma&outcome._&class2_.)))/temp)+
MAX[I]*exp(logcdf('normal',(-MAX[i]+mu&class_.&outcome.[I])/sigma&outcome._&class2_.))

%str(;)

%end;%str(;)
%mend Pred_cnorm;
*  %put   %Pred_cnorm(class=2,outcome=2,equal_sigma=F);;
*  %put   %Pred_cnorm(class=2,outcome=2,equal_sigma=T);;



%macro plot_prep(T,LC,result,order,equal_sigma);
data parameter;set &result. ;keep  parameter estimate;run;

proc transpose data=parameter out=parameter;id parameter;run;

proc sql;
create table data_pred as
select * from base_file_srs, parameter;quit;

data pred_membership_y;set data_pred;
/*Create arrays*/
 %initiation_universal(T=&T.);
 ARRAY MAX[&T.] (92 91 91 92 92 91 91 92 91 91 91 91);
%initiation(T=&T., class=&LC.,outcome=1);
%initiation(T=&T., class=&LC.,outcome=2);

DO I=1 TO &T.;
/*Model prediction*/
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
 %Class_Membership( class=&LC.);

/*Likelihood logLike=sum P(X=t)*P(Y|X=t)*/
 %Posterior(class=&LC.);
run; 

data y_1;
set base_file_srs (keep= Q1HH Q2HH Q3HH Q4HH Q5HH Q6HH Q7HH Q8HH Q9HH Q10HH Q11HH Q12HH);
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
data wide; set wide;quar1=1;
quar2=2;
quar3=3;
quar4=4;
quar5=5;
quar6=6;
quar7=7;
quar8=8;
quar9=9;
quar10=10;
quar11=11;
quar12=12;
run;

data wide; set wide;
/*Create arrays*/
ARRAY X[12] quar1-quar12;
ARRAY MAX[&T.] (92 91 91 92 92 91 91 92 91 91 91 91);
 %initiation_pred(T=&T., class=&LC.,outcome=1);
 %initiation_pred(T=&T., class=&LC.,outcome=2);;

DO I=1 TO 12 ;
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=1);
%model( pattern=X[I]**order*betaoutcome_classorder ,order=&order.,class=&LC.,outcome=2);

%Pred_cnorm(class=&LC.,outcome=1,equal_sigma=&equal_sigma.);
%Pred_cnorm(class=&LC.,outcome=2,equal_sigma=&equal_sigma.);
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


%macro nlmixed_1(T,LC,Y,starting,output,order,equal_sigma);

proc nlmixed data= base_file_srs  itdetails  qpoints=40 noad maxiter=1000 tech=dbldog ;
 bounds 	%bounds_alpha(bounds_alpha=3,class=&LC.),
		%bounds_sigma(bounds_sigma=0,class=&LC.,outcome=&Y.,equal_sigma=&equal_sigma.);
parms 	&starting.
/*%starting_value_alpha(class=&LC.)
		%starting_value_beta_sigma(class=&LC.,outcome=&Y.,order=&order.,equal_sigma=&equal_sigma.)*/
		/*use the previous result as starting point*/;

/*Create arrays*/
%initiation_universal(T=&T.);;
ARRAY MAX[&T.] 92 91 91 92 92 91 91 92 91 91 91 91;

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

%mend nlmixed_1;


%macro nlmixed_MultiTraj(T,LC,starting,output,order,equal_sigma);

proc nlmixed data= base_file_srs  itdetails  qpoints=40 noad maxiter=1000 tech=dbldog ;
 bounds  %bounds_alpha(bounds_alpha=3,class=&LC.) ,
	 %bounds_sigma(bounds_sigma=0,class=&LC.,outcome=1,equal_sigma=&equal_sigma.),
		%bounds_sigma(bounds_sigma=0,class=&LC.,outcome=2,equal_sigma=&equal_sigma.);

parms 	&starting./*use the previous result as starting point*/
		;

/*Create arrays*/
%initiation_universal(T=&T.);;
ARRAY MAX[&T.] 92 91 91 92 92 91 91 92 91 91 91 91;
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
%mend nlmixed_MultiTraj;




/**=========================================== Ref ===================================================**/
/* 
Jones, B. L., & Nagin, D. S. (2007). Advances in Group-Based Trajectory Modeling and an SAS Procedure for Estimating Them. Sociological Methods & Research, 35(4), 542571. https://doi.org/10.1177/0049124106292364
Extension 4
*/

