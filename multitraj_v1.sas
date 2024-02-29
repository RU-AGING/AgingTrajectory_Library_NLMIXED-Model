
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

%include 'D:\Research\Lin\R33\2y\macro_multitraj_v1.sas';



/*=================================Step 1: Load hyper parameter========================================*/
%let class_all = A, B, C, D, E, F, G, H, I, J, K; /* Name of each class. (Don't change) */
%Let class=3; /*Number of latent class*/
%Let order_model=3; /*the degree of the polynomial in the model. 1 - linear, 2 - quardratic, 3 - cubic*/
%Let equal=T; /*Equal variance for each outcome across different classes, input T or F*/


/*==============================Step 2: Run the model for each outcome================================*/
/* Model 1 for HH */
%nlmixed_1(T=12,LC=&class.,	Y=1,
			/*Assign starting value or use existing parameter datafile, the content in starting= will be added after parms in nlmixed*/
			starting=%starting_value_alpha(class=&class.)
					%starting_value_beta_sigma(class=&class.,outcome=1,order=&order_model.,equal_sigma=&equal.),
			/*Define output file name*/
			output=nlm_fix_HH&class.,order=&order_model.,equal_sigma=&equal.);

/* Model 2 for INP */
%nlmixed_1(T=12,LC=&class.,	Y=2,
			starting=%starting_value_alpha(class=&class.)
					%starting_value_beta_sigma(class=&class.,outcome=2,order=&order_model.,equal_sigma=&equal.),
		output=nlm_fix_INP&class.,order=&order_model.,equal_sigma=&equal.);

/*==============================Step 3: Run the model for multitraj model================================*/
/* 3.1 Generate starting value for multitraj model*/
data work.nlm_2y_starting;
 set nlm_fix_HH&class. nlm_fix_INP&class. ;
 if parameter =: 'alpha' then delete ;
 run;

 /* 3.2 Run the multitraj model*/
 /*logLik= sum_t P(X=t)*P(Y1|X=t)*P(Y2|X=t) */
%nlmixed_MultiTraj(T=12,LC=&class.,
			/*Assign starting value or use existing parameter datafile, the content in starting= will be added after parms in nlmixed*/
			starting=%starting_value_alpha(class=&class.)/data=nlm_2y_starting,
			output=nlm_fix_HH_INP&class.,order=&order_model.,equal_sigma=&equal.);
/* 3.3 Plot */
%plot_prep(T=12,LC=&class.,result=nlm_fix_HH_INP&class.,order=&order_model.,equal_sigma=&equal.);




/**=========================================== Ref ===================================================**/
/* 
Jones, B. L., & Nagin, D. S. (2007). Advances in Group-Based Trajectory Modeling and an SAS Procedure for Estimating Them. Sociological Methods & Research, 35(4), 542–571. https://doi.org/10.1177/0049124106292364
Extension 4
*/

