

##########################################################################################################################
# PROJECT NAME: Traj2 Ordinal Cross-Check, R Side
# LAST UPDATED DATE: 24 JUL 2026
# DATA SOURCES: two CSV files written by traj2_replication_ordinal.sas STEP 4. No other input
# PURPOSE: Companion to traj2_replication_ordinal.sas. Produces the lcmm side of the
# cross-implementation comparison in Traj2: A Native Macro Library for Single and Multi-Outcome
# Group-Based Trajectory Modeling in SAS, Zafar, Xia, Lin, Jarrin, Journal of Statistical Software.
# (a)fits the same three-class ordinal model in lcmm
# (b)reports lcmm class proportions and maximised log likelihood
# (c)matches classes across the two programs and reports fitted-trajectory agreement
# (d)draws the Traj2 versus lcmm overlay figure
# COMPARISON BASIS: only quantities that do not depend on parameterisation are compared, namely
# E[Y | class, t], the class proportions, and the maximised log likelihood. Raw thresholds and
# regression coefficients are NOT compared, because the two programs parameterise them differently.
# AUTHOR: Anum Zafar, Weiyi Xia, Haiqun Lin, Olga F. Jarrin
##########################################################################################################################
# Execution Environment: R 4.x with lcmm version 1.9 or later. No other packages required        #
# Run traj2_replication_ordinal.sas FIRST. This script reads the two CSV files it writes         #
# outpath below must match OUTPATH in the SAS script                                             #
##########################################################################################################################
### CODE OVERVIEW ##
# STEP 0:Paths and run controls. Edit only this block
# STEP 1:Read the SAS output
# STEP 2:Fit the same model in lcmm
# STEP 3:lcmm class proportions and log likelihood
# STEP 4:lcmm fitted expected outcome by class
# STEP 5:Match classes across programs, then compare
# STEP 6:Overlay figure
##########################################################################################################################

##########################################################################################################################
# STEP 0: PATHS AND RUN CONTROLS
##########################################################################################################################
# NOTE ON THE PUBLISHED TABLE 15: the paper reports lcmm proportions of 50.4, 33.0 and 16.6 percent
# against Traj2 values of 49.7, 33.5 and 16.8, with -2 log L of 68,520 and 68,517. The current SAS
# cross-check run does NOT return those Traj2 values, so the published pair came from an earlier
# configuration. Regenerate BOTH columns of Table 15 from one run rather than trying to reproduce
# the published figures, and recompute the agreement summary from the new comparison.
##########################################################################################################################

library(lcmm)

outpath <- "T:/datasets_t/replication"   # must match OUTPATH in the SAS script

f_data <- file.path(outpath, "ordinal_data_for_R.csv")
f_traj2 <- file.path(outpath, "traj2_EY_by_class.csv")

KTRUE <- 3    # latent classes
TT    <- 12   # time points
MCAT  <- 4    # ordinal categories, coded 0 1 2 3

##########################################################################################################################
# STEP 1: READ THE SAS OUTPUT
##########################################################################################################################
# ordinal_data_for_R.csv  long format, one row per subject and time point
# traj2_EY_by_class.csv   Traj2 fitted E[Y_t | class] on the 0 to 3 scale
##########################################################################################################################

dat <- read.csv(f_data)
names(dat) <- tolower(names(dat))
dat$id <- dat$bene_id

traj2 <- read.csv(f_traj2)
names(traj2) <- tolower(names(traj2))

##########################################################################################################################
# STEP 2: FIT THE SAME MODEL IN LCMM
##########################################################################################################################
# Ordinal outcome via link = thresholds. The fixed part is a quadratic in time; mixture() puts that
# quadratic under class-specific coefficients, which is the group-based trajectory structure. No
# random effect, matching both the data generating process and the Traj2 fit.
#
# lcmm is sensitive to starting values with more than one class, so the one-class fit seeds the
# multi-class fit through gridsearch, as the package documentation recommends.
##########################################################################################################################

m1 <- lcmm(yval ~ t + I(t^2),
           subject = "id",
           ng = 1,
           link = "thresholds",
           data = dat)

set.seed(20260601)
m3 <- gridsearch(
  rep = 30, maxiter = 15, minit = m1,
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
# CAUTION: pi_lcmm below is the proportion of subjects ASSIGNED to each class by maximum posterior
# probability, not the estimated mixing parameter. The Traj2 column of Table 15 reports the model
# mixing proportions from the softmax of alpha0_*. The two agree closely when classes are well
# separated but they are not the same quantity, so state in the table note which one is reported,
# or take lcmm's estimated mixing proportions from summary(m3) instead.
##########################################################################################################################

pi_lcmm <- as.numeric(table(m3$pprob$class)) / nrow(m3$pprob)

loglik_lcmm  <- m3$loglik
m2loglik_lcmm <- -2 * loglik_lcmm

cat("\nlcmm class proportions:", round(100 * pi_lcmm, 1), "percent\n")
cat("lcmm -2 log L:", round(m2loglik_lcmm, 0), "\n")

##########################################################################################################################
# STEP 4: LCMM FITTED EXPECTED OUTCOME BY CLASS
##########################################################################################################################
# For a thresholds link, predictY on the natural outcome scale returns the expected value of the
# ordinal outcome. Evaluated on a regular grid t = 1 to TT. py$pred carries one column per class,
# in lcmm's own class order, which STEP 5 then matches to the Traj2 order.
##########################################################################################################################

newd <- data.frame(t = 1:TT)
py   <- predictY(m3, newd, draws = FALSE)

ey_lcmm <- data.frame(
  t     = rep(1:TT, times = KTRUE),
  class = rep(1:KTRUE, each = TT),
  EY_lcmm = as.vector(as.matrix(py$pred)[, 1:KTRUE])
)

##########################################################################################################################
# STEP 5: MATCH CLASSES ACROSS PROGRAMS, THEN COMPARE
##########################################################################################################################
# Classes carry no inherent order in a finite mixture, so they are matched by minimising the total
# squared distance between the two programs' fitted trajectories, which is what the paper describes.
#
# SAMPLE DATA OUTPUT (table15_EY_comparison.csv):
# class   t   EY_traj2   EY_lcmm   absdiff
# 1       1     1.8721    1.8698    0.0023
##########################################################################################################################

t2 <- traj2[, c("classnum", "t", "ey_traj2")]
names(t2) <- c("class_t2", "t", "EY_traj2")

perms <- as.matrix(expand.grid(rep(list(1:KTRUE), KTRUE)))
perms <- perms[apply(perms, 1, function(r) length(unique(r)) == KTRUE), , drop = FALSE]

sse_of <- function(map) {
  tot <- 0
  for (k in 1:KTRUE) {
    a <- t2$EY_traj2[t2$class_t2 == k]
    b <- ey_lcmm$EY_lcmm[ey_lcmm$class == map[k]]
    tot <- tot + sum((a - b)^2)
  }
  tot
}
best <- perms[which.min(apply(perms, 1, sse_of)), ]

cmp <- do.call(rbind, lapply(1:KTRUE, function(k) {
  data.frame(
    class    = k,
    t        = 1:TT,
    EY_traj2 = t2$EY_traj2[t2$class_t2 == k],
    EY_lcmm  = ey_lcmm$EY_lcmm[ey_lcmm$class == best[k]]
  )
}))
cmp$absdiff <- abs(cmp$EY_traj2 - cmp$EY_lcmm)

cat("\n--- Table 15: fitted E[Y_t | class] agreement ---\n")
cat("max absolute difference:", round(max(cmp$absdiff), 4), "\n")
cat("mean absolute difference:", round(mean(cmp$absdiff), 4), "\n")
cat("across", KTRUE * TT, "class-by-time points\n")

cat("\nlcmm proportions in Traj2 class order:",
    round(100 * pi_lcmm[best], 1), "percent\n")

write.csv(cmp, file.path(outpath, "table15_EY_comparison.csv"), row.names = FALSE)

##########################################################################################################################
# STEP 6: OVERLAY FIGURE
##########################################################################################################################
# Traj2 as solid lines with markers, lcmm as dashed, one colour per matched class pair.
##########################################################################################################################

pdf(file.path(outpath, "fig2_ordinal_traj2_vs_lcmm.pdf"), width = 7, height = 5)
plot(range(cmp$t), range(c(cmp$EY_traj2, cmp$EY_lcmm)), type = "n",
     xlab = "Time", ylab = "E[Y | class] (0-3 scale)",
     main = "Ordinal GBTM: Traj2 (solid) vs lcmm (dashed)")
cols <- c("darkgreen", "orange3", "purple3")
for (k in 1:KTRUE) {
  s <- cmp[cmp$class == k, ]
  lines(s$t, s$EY_traj2, col = cols[k], lwd = 2, lty = 1)
  lines(s$t, s$EY_lcmm,  col = cols[k], lwd = 2, lty = 2)
  points(s$t, s$EY_traj2, col = cols[k], pch = 16, cex = 0.6)
}
legend("bottomright", legend = paste("Class", 1:KTRUE),
       col = cols, lwd = 2, bty = "n")
legend("topleft", legend = c("Traj2 (solid)", "lcmm (dashed)"),
       lty = c(1, 2), bty = "n")
dev.off()

cat("\nWrote table15_EY_comparison.csv and fig2_ordinal_traj2_vs_lcmm.pdf to",
    outpath, "\n")

##########################################################################################################################
# END
##########################################################################################################################
# OUTPUT FILES, both written to outpath
# table15_EY_comparison.csv        one row per class and time point. EY_traj2, EY_lcmm, absdiff
# fig2_ordinal_traj2_vs_lcmm.pdf   the Traj2 versus lcmm overlay
#
# CONSOLE OUTPUT TO CARRY INTO TABLE 15
#   lcmm class proportions, reported in Traj2 class order
#   lcmm -2 log L
#   maximum and mean absolute difference in fitted E[Y | class]
##########################################################################################################################