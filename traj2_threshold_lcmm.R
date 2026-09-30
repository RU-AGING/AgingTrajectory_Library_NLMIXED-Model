#Copyright (c) 2026 Community Health and Aging Outcomes (CHAO) Lab, Rutgers University.
#Released under the MIT License. Full text in LICENSE at
#https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model;
##########################################################################################################################
# PROJECT NAME: Traj2 Class-Specific Threshold Check, R Side
# LAST UPDATED DATE: 30 SEP 2026 (pi fix, wider gridsearch, duplicate-class check, Figure 6 layout)
# DATA SOURCES: four CSV files written by traj2_threshold_check.sas STEP 2 and STEP 3. No other input
# PURPOSE: Companion to traj2_threshold_check.sas. Reproduces Table 16 and Figure 6 of Traj2:
# Censored-normal and ordinal-probit group-based trajectory modeling in SAS for restricted data
# environments, Zafar, Xia, Lin, Jarrin, MethodsX, Section 3.3. Fits the same data in lcmm, whose
# ordinal link uses one threshold set for all classes, and compares Traj2 and lcmm against the TRUTH:
# (a)class proportions from each program next to the true 0.45 / 0.35 / 0.20
# (b)worst-class and mean RMSE of fitted E[Y_t | class] against the true E[Y_t | class]
# (c)likelihood ratio test, lcmm nested in Traj2, on (C-1)(M-2)=4 degrees of freedom
# (d)overlay figure: truth, Traj2 and lcmm by class
# AUTHOR: Anum Zafar, Weiyi Xia, Haiqun Lin, Olga F. Jarrin
##########################################################################################################################
# Execution Environment: R 4.x with lcmm version 1.9 or later. No other packages required        #
# Run traj2_threshold_check.sas FIRST with RUN_SINGLE=1                                          #
# outpath below must match OUTPATH in the SAS script                                             #
##########################################################################################################################
### CODE OVERVIEW ##
# STEP 0:Paths and run controls. Edit only this block
# STEP 1:Read the SAS output
# STEP 2:Fit the same model in lcmm
# STEP 3:lcmm class proportions and log likelihood
# STEP 4:lcmm fitted expected outcome by class
# STEP 5:Match lcmm classes to the true classes, then compare both programs to the truth
# STEP 6:Overlay figure
##########################################################################################################################

##########################################################################################################################
# STEP 0: PATHS AND RUN CONTROLS
##########################################################################################################################
library(lcmm)

outpath <- "T:/datasets_t/replication_1"   # must match OUTPATH in the SAS script

KTRUE <- 3    # latent classes
TT    <- 12   # time points
MCAT  <- 4    # ordinal categories, coded 0 1 2 3
PI_TRUE <- c(0.45, 0.35, 0.20)   # true class proportions, true class order

##########################################################################################################################
# STEP 1: READ THE SAS OUTPUT
##########################################################################################################################
# cst_data_for_R.csv  long format, one row per subject and time point
# cst_traj2_EY.csv    Traj2 fitted E[Y_t | class], already matched to the true class by SAS
# cst_traj2_fit.csv   Traj2 -2 log L
# cst_true_EY.csv     true E[Y_t | class]
##########################################################################################################################
dat <- read.csv(file.path(outpath, "cst_data_for_R.csv"))
names(dat) <- tolower(names(dat))
dat$id <- dat$bene_id

traj2 <- read.csv(file.path(outpath, "cst_traj2_EY.csv"))
names(traj2) <- tolower(names(traj2))

truth <- read.csv(file.path(outpath, "cst_true_EY.csv"))
names(truth) <- tolower(names(truth))

m2ll_traj2 <- read.csv(file.path(outpath, "cst_traj2_fit.csv"))[1, 1]

##########################################################################################################################
# STEP 2: FIT THE SAME MODEL IN LCMM
##########################################################################################################################
# Same call as traj2_lcmm_crosscheck.R: quadratic in time, class-specific coefficients through
# mixture(), thresholds link, no random effect, one-class fit seeding a gridsearch.
##########################################################################################################################
m1 <- lcmm(yval ~ t + I(t^2),
           subject = "id",
           ng = 1,
           link = "thresholds",
           data = dat)

# rep=30 settled on a local maximum at -2 log L 71,344 in which two classes had identical
# coefficients. 100 starts of 50 iterations each are used so the reported lcmm fit is its best.
set.seed(20260930)
m3 <- gridsearch(
  rep = 100, maxiter = 50, minit = m1,
  lcmm(yval ~ t + I(t^2),
       mixture = ~ t + I(t^2),
       subject = "id",
       ng = KTRUE,
       link = "thresholds",
       data = dat)
)

summary(m3)

##########################################################################################################################
# STEP 3: LCMM CLASS PROPORTIONS AND LOG LIKELIHOOD
##########################################################################################################################
# Estimated mixing proportions from the class-membership intercepts, so the quantity matches the
# Traj2 column (softmax of alpha0_*). lcmm uses the last class as the reference.
##########################################################################################################################
# the class-membership intercepts are the first KTRUE-1 entries of m3$best. Selecting them by
# name also picks up the longitudinal "intercept class2/3" terms, so select by position
b_int <- m3$best[1:(KTRUE - 1)]
pi_lcmm <- exp(c(b_int, 0)) / sum(exp(c(b_int, 0)))

m2ll_lcmm <- -2 * m3$loglik

# duplicate-class check: in the summary(m3) output above, no two classes should have the same
# intercept, t and t^2 coefficients. If two do, lcmm stopped at a degenerate local maximum; rerun
# with more starts (larger rep in STEP 2) before using the numbers

##########################################################################################################################
# STEP 4: LCMM FITTED EXPECTED OUTCOME BY CLASS
##########################################################################################################################
py <- predictY(m3, data.frame(t = 1:TT), draws = FALSE)

ey_lcmm <- data.frame(
  t       = rep(1:TT, times = KTRUE),
  class   = rep(1:KTRUE, each = TT),
  EY_lcmm = as.vector(as.matrix(py$pred)[, 1:KTRUE])
)

##########################################################################################################################
# STEP 5: MATCH LCMM CLASSES TO THE TRUE CLASSES, THEN COMPARE BOTH PROGRAMS TO THE TRUTH
##########################################################################################################################
# Traj2 classes arrive already matched to the true classes. lcmm classes are matched here by the
# same rule, the smallest total squared distance between fitted and true E[Y_t | class].
#
# SAMPLE DATA OUTPUT (cst_EY_comparison.csv):
# trueclass   t   EY_true   EY_traj2   EY_lcmm
# 1           1    0.8912     0.8870    1.0431
##########################################################################################################################
perms <- as.matrix(expand.grid(rep(list(1:KTRUE), KTRUE)))
perms <- perms[apply(perms, 1, function(r) length(unique(r)) == KTRUE), , drop = FALSE]

sse_of <- function(map) {
  tot <- 0
  for (k in 1:KTRUE) {
    a <- truth$ey_true[truth$trueclass == k]
    b <- ey_lcmm$EY_lcmm[ey_lcmm$class == map[k]]
    tot <- tot + sum((a - b)^2)
  }
  tot
}
best <- perms[which.min(apply(perms, 1, sse_of)), ]

cmp <- do.call(rbind, lapply(1:KTRUE, function(k) {
  data.frame(
    trueclass = k,
    t         = 1:TT,
    EY_true   = truth$ey_true[truth$trueclass == k],
    EY_traj2  = traj2$ey_traj2[traj2$trueclass == k],
    EY_lcmm   = ey_lcmm$EY_lcmm[ey_lcmm$class == best[k]]
  )
}))

rmse_by <- function(col) sapply(1:KTRUE, function(k) {
  s <- cmp[cmp$trueclass == k, ]
  sqrt(mean((s[[col]] - s$EY_true)^2))
})
rmse_traj2 <- rmse_by("EY_traj2")
rmse_lcmm  <- rmse_by("EY_lcmm")

pi_traj2 <- sapply(1:KTRUE, function(k) traj2$pi_traj2[traj2$trueclass == k][1])

lrt <- m2ll_lcmm - m2ll_traj2
df_lrt <- (KTRUE - 1) * (MCAT - 2)
p_lrt <- pchisq(lrt, df = df_lrt, lower.tail = FALSE)

cst_table <- data.frame(
  trueclass  = 1:KTRUE,
  pi_true    = PI_TRUE,
  pi_traj2   = round(pi_traj2, 3),
  pi_lcmm    = round(pi_lcmm[best], 3),
  rmse_traj2 = round(rmse_traj2, 4),
  rmse_lcmm  = round(rmse_lcmm, 4)
)

cat("\n--- Table 16: class-specific thresholds, Traj2 versus lcmm against the truth ---\n")
print(cst_table, row.names = FALSE)
cat("\n-2 log L  Traj2:", round(m2ll_traj2, 1), "  lcmm:", round(m2ll_lcmm, 1), "\n")
cat("LRT (lcmm nested in Traj2):", round(lrt, 1), "on", df_lrt, "df, p =", signif(p_lrt, 3), "\n")
cat("worst-class RMSE of E[Y|class] vs truth  Traj2:", round(max(rmse_traj2), 4),
    "  lcmm:", round(max(rmse_lcmm), 4), "\n")

write.csv(cmp, file.path(outpath, "cst_EY_comparison.csv"), row.names = FALSE)
write.csv(cst_table, file.path(outpath, "cst_table.csv"), row.names = FALSE)

##########################################################################################################################
# STEP 6: OVERLAY FIGURE
##########################################################################################################################
# One panel per true class. Truth grey, Traj2 solid, lcmm dashed.
##########################################################################################################################
pdf(file.path(outpath, "fig_cst_traj2_lcmm_truth.pdf"), width = 7.5, height = 3)
par(mfrow = c(1, KTRUE), mar = c(4, 4, 2, 0.5), cex = 0.9)
cols <- c("darkgreen", "orange3", "purple3")
for (k in 1:KTRUE) {
  s <- cmp[cmp$trueclass == k, ]
  plot(s$t, s$EY_true, type = "l", lwd = 5, col = "grey75", ylim = c(0, 3),
       xlab = "Quarter (t)", ylab = if (k == 1) "E[Y | class]" else "",
       main = paste("Class", k), xaxt = "n")
  axis(1, at = c(1, 4, 8, 12))
  lines(s$t, s$EY_traj2, col = cols[k], lwd = 2, lty = 1)
  lines(s$t, s$EY_lcmm,  col = cols[k], lwd = 2, lty = 2)
  if (k == 1) legend("topleft", legend = c("Truth", "Traj2", "lcmm"),
                     col = c("grey75", "black", "black"), lwd = c(5, 2, 2),
                     lty = c(1, 1, 2), bty = "n", cex = 0.9)
}
dev.off()

##########################################################################################################################
# END
##########################################################################################################################
# OUTPUT FILES, all written to outpath
# cst_EY_comparison.csv            one row per true class and time. EY_true, EY_traj2, EY_lcmm
# cst_table.csv                    class proportions and E[Y|class] RMSE for both programs
# fig_cst_traj2_lcmm_truth.pdf     truth, Traj2 and lcmm by class
#
# CONSOLE OUTPUT TO CARRY INTO THE PAPER
#   the cst_table rows, both -2 log L values and the likelihood ratio test
##########################################################################################################################
