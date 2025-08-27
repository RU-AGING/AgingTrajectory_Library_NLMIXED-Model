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
/*============================================================Step 0: Load macros========================================================*/
*****************************************************************************************************************************************/;
* Assign libnames to libraries in SAS; 
options fullstimer;
run;
libname DUA_SPECIFIC "/sas/vrdc/data/dua/DUA_SPECIFIC/source";
*******************************************************************************************************************************************;

/*=========================================================Step 1: Load hyper parameter====================================================*/
%Let class=3 ; /*Assign number of latent class*/
%Let order_model=2; /*the degree of the polynomial in the model. 1 - linear, 2 - quadratic, 3 - cubic*/
%Let equal=T; /*Equal variance for each outcome across different classes, input T or F*/
%let dist=zip;            /* distribution: normal | poisson | zip */

/*=========================================================Step 2: Run Models=============================================*/

/* Fit models */
%nlmixed_1(T=12,LC=&class.,Y=1,
           starting=%starting_value_alpha(class=&class.)
                    %starting_value_beta_sigma(class=&class.,outcome=1,order=&order_model.,equal_sigma=&equal.),
           output=nlm_fix_T1&class.,order=&order_model.,equal_sigma=&equal.,dist=&dist.);

%nlmixed_1(T=12,LC=&class.,Y=2,
           starting=%starting_value_alpha(class=&class.)
                    %starting_value_beta_sigma(class=&class.,outcome=2,order=&order_model.,equal_sigma=&equal.),
           output=nlm_fix_T2&class.,order=&order_model.,equal_sigma=&equal.,dist=&dist.);

%nlmixed_MultiTraj(T=12,LC=&class.,
           starting=%starting_value_alpha(class=&class.)/data=nlm_2y_starting,
           output=nlm_fix_T1_T2&class.,order=&order_model.,equal_sigma=&equal.,dist=&dist.);

%plot_prep(T=12,LC=&class.,result=nlm_fix_T1_T2&class.,order=&order_model.,equal_sigma=&equal.,dist=&dist.);





