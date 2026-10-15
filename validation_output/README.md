# Validation output

Output files from the replication scripts, kept so the published tables and
figures can be checked without rerunning SAS or R.

All results come from simulated data with fixed seeds, so rerunning the scripts
reproduces these files exactly.

## ordinal/

From `traj2_replication_ordinal.sas` STEP 4 and `traj2_lcmm_crosscheck.R`.
Seed `SEED_CC=20260601`. R 4.6.1 with lcmm 2.2.2.

| File | What it holds |
|---|---|
| `ordinal_data_for_R.csv` | Long-format simulated data, the input to `traj2_lcmm_crosscheck.R` |
| `traj2_EY_by_class.csv` | Traj2 fitted E[Y_t given class] and mixing proportions, the Traj2 column of MethodsX Table 15 |
| `table15_EY_comparison.csv` | Traj2 and lcmm fitted E[Y_t given class] side by side, 36 class-by-time points, with the absolute difference |
| `fig5_ordinal_traj2_vs_lcmm.pdf` | MethodsX Figure 5 |

## thresholds/

From `traj2_threshold_check.sas` and `traj2_threshold_lcmm.R`.
Seeds `SEED_CC=20260930` and `SEED_MC=20261000`. R 4.6.1 with lcmm 2.2.2.

| File | What it holds |
|---|---|
| `cst_data_for_R.csv` | Long-format simulated data, the input to `traj2_threshold_lcmm.R` |
| `cst_recovery_metrics.csv` | One row per class per replication, the per-replication metrics behind MethodsX Table 16 |
| `traj2_threshold_recovery.pdf` | Threshold recovery summary, MethodsX Table 16 |
| `cst_traj2_EY.csv` | Traj2 fitted E[Y_t given class] on the single data set |
| `cst_true_EY.csv` | True E[Y_t given class] under the generating parameters |
| `cst_traj2_fit.csv` | Traj2 -2 log L for the single data set |
| `cst_EY_comparison.csv` | Truth, Traj2 and lcmm side by side, the data behind MethodsX Figure 6 |
| `cst_table.csv` | Class proportions and RMSE for both programs, MethodsX Table 17 and CMPB Table 4 |
| `fig_cst_traj2_lcmm_truth.pdf` | MethodsX Figure 6 |
