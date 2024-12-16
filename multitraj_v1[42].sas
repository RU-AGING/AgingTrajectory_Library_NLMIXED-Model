/*****************************************************************************************************************************************
* Community Health and Aging Outcome (CHAO) Lab - Rutgers, The State University of New Jersey                                           *
* Title: Group-Based Trajectory Modeling using NLMIX                                                                                    *
* Purpose: This code implements group-based trajectory modeling with joint trajectories for two outcomes (T1 and T2) using the          *
*          SAS NLMIXED procedure. It includes single-outcome models, a multi-trajectory model, and a plotting macro for visualization.   *
* Data Sources: Simulation or real-world longitudinal datasets containing repeated measures for multiple time points (e.g., Q=12).      *
* Cohort: Individuals with repeated measures data for two outcomes over N time periods.                                   		*
* Outputs:                                                                                                                              *
*   1. Single-outcome trajectory model results: nlm_fix_T1&class. and nlm_fix_T2&class.                                                 *
*   2. Multitrajectory model results: nlm_fix_T1_T2&class.                                                                              *
*   3. Visualizations of predicted vs. observed trajectories and averaged posterior class memberships.                                  *
* Author: Weiyi Xia                                                                                                                     *
* Ref: Jones, B. L., & Nagin, D. S. (2007). Advances in Group-Based Trajectory Modeling and an SAS Procedure for Estimating Them.       *
*      Sociological Methods & Research, 35(4), 542–571. https://doi.org/10.1177/0049124106292364                                        *
*****************************************************************************************************************************************/;


/*============================================================Step 0: Load macros========================================================*/
*****************************************************************************************************************************************/;
* Assign libnames to libraries in SAS; 
options fullstimer;
run;
libname DUA_SPECIFIC "/sas/vrdc/data/dua/DUA_SPECIFIC/source";
*******************************************************************************************************************************************;

/*=========================================================Step 1: Load hyper parameter====================================================*/
%let class_all = A, B, C, D, E, F, G, H, I, J, K; /* Name of each class */
%Let class= ; /*Assign number of latent class*/
%Let order_model=3; /*the degree of the polynomial in the model. 1 - linear, 2 - quadratic, 3 - cubic*/
%Let equal=T; /*Equal variance for each outcome across different classes, input T or F*/


/*=========================================================Step 2: Run the model for each outcome=============================================*/
/* Model 1 for Outcome 1 */
%nlmixed_1(T=12,LC=&class.,Y=1,
			/*Assign starting value or use existing parameter datafile, the content in starting= will be added after parms in nlmixed*/
			starting=%starting_value_alpha(class=&class.)
					%starting_value_beta_sigma(class=&class.,outcome=1,order=&order_model.,equal_sigma=&equal.),
			/*Define output file name*/
			output=nlm_fix_T1&class.,order=&order_model.,equal_sigma=&equal.);

/* Model 2 for Outcome 2 */
%nlmixed_1(T=12,LC=&class.,	Y=2,
			starting=%starting_value_alpha(class=&class.)
					%starting_value_beta_sigma(class=&class.,outcome=2,order=&order_model.,equal_sigma=&equal.),
		output=nlm_fix_T2&class.,order=&order_model.,equal_sigma=&equal.);

/*=========================================================Step 3: Run the model for multitraj model=============================================*/
*Combines results from the two individual models, excluding alpha parameters, to initialize the multitrajectory model;
/* 3.1 Generate starting value for multitraj model*/
data work.nlm_2y_starting;
 set nlm_fix_T1&class. nlm_fix_T2&class. ;
 if parameter =: 'alpha' then delete ;
 run;

 /* 3.2 Run the multitraj model*/
 * Fits a multitrajectory model using starting values derived from earlier models;
%nlmixed_MultiTraj(T=12,LC=&class.,
			/*Assign starting value or use existing parameter datafile, the content in starting= will be added after parms in nlmixed*/
			starting=%starting_value_alpha(class=&class.)/data=nlm_2y_starting,
			output=nlm_fix_T1_T2&class.,order=&order_model.,equal_sigma=&equal.);
/* 3.3 Plot */
*Prepares data for visualization of the multitrajectory results;
%plot_prep(T=12,LC=&class.,result=nlm_fix_T1_T2&class.,order=&order_model.,equal_sigma=&equal.);



