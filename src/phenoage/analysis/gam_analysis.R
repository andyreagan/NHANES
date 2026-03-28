#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════
# GAM survival analysis for PhenoAge biomarkers
#
# Uses mgcv::gam() with the cox.ph family to fit an additive Cox PH
# model with smooth (spline) terms on each biomarker.  This captures
# non-linear (U-shaped, J-shaped) relationships that the linear
# Gompertz model cannot.
#
# Why GAMs?
#   - Each smooth term is biologically interpretable: you can plot
#     the partial effect of albumin and SEE the U-shape.
#   - Penalized splines prevent overfitting (automatic smoothness
#     selection via REML or GCV).
#   - bam() scales to large datasets; gam() is fine for n ≈ 10k.
#
# Usage:
#   Rscript src/phenoage/gam_analysis.R
#   # or via Make:
#   make gam-analysis
# ══════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(mgcv)
  library(survival)
})

cat("========================================================================\n")
cat("  GAM SURVIVAL ANALYSIS (mgcv::gam with cox.ph family)\n")
cat("========================================================================\n\n")

# ──────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ──────────────────────────────────────────────────────────────────────

train <- read.csv("research_agents/phenoage_train_2015.csv")
cat("Training data: n =", nrow(train), ", deaths =", sum(train$event), "\n")

# Features
features <- c("albumin_gL", "creat_umol", "glucose_mmol", "lncrp",
               "lymph", "mcv", "rdw", "alp", "wbc", "age")

# Add sex
train$sex <- ifelse(train$HSSEX == 2, 1, 0)

# Ensure no NAs in features + outcome
complete <- complete.cases(train[, c(features, "sex", "time", "event")])
train <- train[complete, ]
cat("After completeness filter: n =", nrow(train), ", deaths =", sum(train$event), "\n\n")

# ──────────────────────────────────────────────────────────────────────
# 2. FIT COX PH GAM (all smooths)
# ──────────────────────────────────────────────────────────────────────

# mgcv::cox.ph expects:
#   response = time (event time)
#   weights  = event indicator (1 = died, 0 = censored)

cat("Fitting GAM with smooth terms on all 10 PhenoAge biomarkers + sex...\n")
cat("  (Using gam() with cox.ph family, REML smoothness selection)\n\n")

# Convert time from months to years for nicer interpretation
train$time_yr <- train$time / 12

# k = basis dimension. Default k=10 is usually fine for these sample sizes.
# We use k=10 for biomarkers and k=20 for age (which has a wider range
# and a more complex relationship with mortality).

t0 <- proc.time()
gam_full <- gam(
  time_yr ~
    s(albumin_gL, k = 10) +
    s(creat_umol, k = 10) +
    s(glucose_mmol, k = 10) +
    s(lncrp, k = 10) +
    s(lymph, k = 10) +
    s(mcv, k = 10) +
    s(rdw, k = 10) +
    s(alp, k = 10) +
    s(wbc, k = 10) +
    s(age, k = 20) +
    sex,
  family = cox.ph(),
  data = train,
  weights = event,
  method = "REML"
)
t1 <- proc.time()
cat("  GAM fit complete in", round((t1 - t0)[3], 1), "seconds\n\n")

# ──────────────────────────────────────────────────────────────────────
# 3. MODEL SUMMARY
# ──────────────────────────────────────────────────────────────────────

cat("========================================================================\n")
cat("  GAM SUMMARY\n")
cat("========================================================================\n\n")
print(summary(gam_full))

# ──────────────────────────────────────────────────────────────────────
# 4. EFFECTIVE DEGREES OF FREEDOM (non-linearity check)
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  SMOOTH TERM EFFECTIVE DEGREES OF FREEDOM\n")
cat("  (edf ≈ 1 means linear; edf > 1 means non-linear)\n")
cat("========================================================================\n\n")

edf_table <- data.frame(
  term = rownames(summary(gam_full)$s.table),
  edf = summary(gam_full)$s.table[, "edf"],
  ref_df = summary(gam_full)$s.table[, "Ref.df"],
  chi_sq = summary(gam_full)$s.table[, "Chi.sq"],
  p_value = summary(gam_full)$s.table[, "p-value"]
)
# Sort by edf descending
edf_table <- edf_table[order(-edf_table$edf), ]
for (i in seq_len(nrow(edf_table))) {
  r <- edf_table[i, ]
  lin_flag <- ifelse(r$edf < 1.5, " (≈ linear)", ifelse(r$edf > 3, " (highly non-linear)", ""))
  cat(sprintf("  %-20s  edf = %5.2f  χ² = %8.1f  p = %.1e%s\n",
              r$term, r$edf, r$chi_sq, r$p_value, lin_flag))
}

# ──────────────────────────────────────────────────────────────────────
# 5. COMPARISON: GAM vs LINEAR COX PH
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  LINEAR COX PH COMPARISON\n")
cat("========================================================================\n\n")

# Fit standard linear Cox PH for comparison
cox_linear <- coxph(
  Surv(time, event) ~
    albumin_gL + creat_umol + glucose_mmol + lncrp +
    lymph + mcv + rdw + alp + wbc + age + sex,
  data = train
)
cat("Linear Cox PH:\n")
print(summary(cox_linear))

# AIC comparison (note: AIC for gam cox.ph uses deviance-based approximation)
cat("\n  Model comparison:\n")
cat(sprintf("    Linear Cox PH:  AIC = %.1f\n", AIC(cox_linear)))
cat(sprintf("    GAM Cox PH:     AIC = %.1f\n", AIC(gam_full)))
cat(sprintf("    ΔAIC = %.1f (negative = GAM better)\n",
            AIC(gam_full) - AIC(cox_linear)))

# ──────────────────────────────────────────────────────────────────────
# 6. CONCORDANCE INDEX (in-sample and cross-validated)
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  CONCORDANCE INDEX\n")
cat("========================================================================\n\n")

# Linear Cox C-index (built-in)
lin_conc <- cox_linear$concordance
cat(sprintf("  Linear Cox PH C-index:  %.4f (se=%.4f)\n",
            lin_conc["concordance"], lin_conc["std(c)"]))

# GAM: compute C-index from linear predictor
gam_lp <- predict(gam_full, type = "link")
# mgcv cox.ph: higher linear predictor = higher risk
# concordance() expects higher risk = higher predictor, which is the default
surv_obj <- Surv(train$time, train$event)
gam_conc <- concordance(surv_obj ~ gam_lp)
# concordance() returns c-index where 1 = perfect concordance
# but we need to check the sign convention
cat(sprintf("  GAM Cox PH C-index:     %.4f (se=%.4f)\n",
            gam_conc$concordance, sqrt(gam_conc$var)))
# If C < 0.5, the sign convention is flipped; use 1 - C
if (gam_conc$concordance < 0.5) {
  cat("  (Note: concordance < 0.5 indicates reversed sign; using 1-C)\n")
  gam_cindex <- 1 - gam_conc$concordance
} else {
  gam_cindex <- gam_conc$concordance
}

# ── 5-fold cross-validated C-index ──
cat("\n  5-fold cross-validated C-index:\n")
set.seed(42)
n <- nrow(train)
folds <- sample(rep(1:5, length.out = n))

cv_cindex_gam <- numeric(5)
cv_cindex_lin <- numeric(5)

for (k in 1:5) {
  trn <- train[folds != k, ]
  tst <- train[folds == k, ]

  # GAM
  gam_k <- tryCatch(
    gam(
      time_yr ~
        s(albumin_gL, k = 10) + s(creat_umol, k = 10) +
        s(glucose_mmol, k = 10) + s(lncrp, k = 10) +
        s(lymph, k = 10) + s(mcv, k = 10) +
        s(rdw, k = 10) + s(alp, k = 10) +
        s(wbc, k = 10) + s(age, k = 20) + sex,
      family = cox.ph(),
      data = trn,
      weights = event,
      method = "REML"
    ),
    error = function(e) NULL
  )

  if (!is.null(gam_k)) {
    pred_gam <- predict(gam_k, newdata = tst, type = "link")
    conc_gam <- concordance(Surv(tst$time, tst$event) ~ pred_gam)
    c_gam <- conc_gam$concordance
    if (c_gam < 0.5) c_gam <- 1 - c_gam
    cv_cindex_gam[k] <- c_gam
  } else {
    cv_cindex_gam[k] <- NA
  }

  # Linear Cox
  cox_k <- coxph(
    Surv(time, event) ~
      albumin_gL + creat_umol + glucose_mmol + lncrp +
      lymph + mcv + rdw + alp + wbc + age + sex,
    data = trn
  )
  pred_lin <- predict(cox_k, newdata = tst, type = "lp")
  conc_lin <- concordance(Surv(tst$time, tst$event) ~ pred_lin)
  c_lin <- conc_lin$concordance
  if (c_lin < 0.5) c_lin <- 1 - c_lin
  cv_cindex_lin[k] <- c_lin

  cat(sprintf("    Fold %d: GAM=%.4f  Linear=%.4f\n",
              k, cv_cindex_gam[k], cv_cindex_lin[k]))
}

cat(sprintf("\n    CV Mean:  GAM=%.4f (sd=%.4f)  Linear=%.4f (sd=%.4f)\n",
            mean(cv_cindex_gam, na.rm = TRUE), sd(cv_cindex_gam, na.rm = TRUE),
            mean(cv_cindex_lin, na.rm = TRUE), sd(cv_cindex_lin, na.rm = TRUE)))

# ──────────────────────────────────────────────────────────────────────
# 7. PARTIAL EFFECT PLOTS (smooth curves)
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  GENERATING PARTIAL EFFECT PLOTS\n")
cat("========================================================================\n\n")

outdir <- "output/phenoage_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Save smooth curves as CSV for plotting in Python/Altair
# (or use the PDF plot for quick visualization)

# PDF of all smooth terms
pdf_path <- file.path(outdir, "gam_smooth_effects.pdf")
pdf(pdf_path, width = 12, height = 10)
par(mfrow = c(3, 4))
plot(gam_full, shade = TRUE, shade.col = "lightblue",
     seWithMean = TRUE, scale = 0, pages = 0)
dev.off()
cat("  Saved smooth effects plot:", pdf_path, "\n")

# Export smooth curve data as CSV for each term
smooth_data <- list()
for (feat in c(features, "age")) {
  if (feat == "age" && "age" %in% features) next  # avoid duplicate

  # Create prediction grid
  feat_vals <- seq(
    quantile(train[[feat]], 0.01, na.rm = TRUE),
    quantile(train[[feat]], 0.99, na.rm = TRUE),
    length.out = 200
  )

  # Newdata: set all other variables to their median, vary this one
  newdata <- data.frame(matrix(0, nrow = 200, ncol = length(features) + 1))
  colnames(newdata) <- c(features, "sex")
  for (f in c(features, "sex")) {
    newdata[[f]] <- median(train[[f]], na.rm = TRUE)
  }
  newdata[[feat]] <- feat_vals
  # time_yr needed for prediction (won't affect link prediction)
  newdata$time_yr <- median(train$time_yr)

  # Predict on the link (log-hazard ratio) scale
  pred <- predict(gam_full, newdata = newdata, type = "terms",
                  se.fit = TRUE)

  # Find the column corresponding to this feature's smooth
  smooth_col <- grep(feat, colnames(pred$fit), value = TRUE)
  if (length(smooth_col) == 0) next

  smooth_data[[feat]] <- data.frame(
    feature = feat,
    value = feat_vals,
    effect = pred$fit[, smooth_col[1]],
    se = pred$se.fit[, smooth_col[1]]
  )
}

all_smooth <- do.call(rbind, smooth_data)
csv_path <- file.path(outdir, "gam_smooth_curves.csv")
write.csv(all_smooth, csv_path, row.names = FALSE)
cat("  Saved smooth curve data:", csv_path, "\n")

# ──────────────────────────────────────────────────────────────────────
# 8. CONCURVITY CHECK
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  CONCURVITY (analogue of multicollinearity for GAMs)\n")
cat("========================================================================\n\n")

conc_est <- concurvity(gam_full, full = TRUE)
cat("  Estimated concurvity (worst case):\n")
print(round(conc_est, 3))

conc_obs <- concurvity(gam_full, full = FALSE)
cat("\n  Observed pairwise concurvity:\n")
print(round(conc_obs$estimate, 3))

# ──────────────────────────────────────────────────────────────────────
# 9. PENALIZED COX PH (survival::coxph with pspline) FOR COMPARISON
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  PENALIZED COX PH WITH PSPLINES (survival::coxph + pspline)\n")
cat("========================================================================\n\n")

# Note: bam() does not support cox.ph family.
# Instead, compare with coxph + pspline() which also fits penalized splines
# but within the survival package framework.

t0 <- proc.time()
cox_pspline <- coxph(
  Surv(time, event) ~
    pspline(albumin_gL, df = 4) +
    pspline(creat_umol, df = 4) +
    pspline(glucose_mmol, df = 4) +
    pspline(lncrp, df = 4) +
    pspline(lymph, df = 4) +
    pspline(mcv, df = 4) +
    pspline(rdw, df = 4) +
    pspline(alp, df = 4) +
    pspline(wbc, df = 4) +
    pspline(age, df = 6) +
    sex,
  data = train
)
t1 <- proc.time()
cat("  Penalized Cox PH fit in", round((t1 - t0)[3], 1), "seconds\n\n")

pspline_conc <- cox_pspline$concordance
cat(sprintf("  Penalized Cox PH C-index: %.4f\n", pspline_conc["concordance"]))
cat(sprintf("  Penalized Cox PH AIC:     %.1f\n", AIC(cox_pspline)))

# ──────────────────────────────────────────────────────────────────────
# 10. INTERACTION MODELS (age × biomarker)
# ──────────────────────────────────────────────────────────────────────

cat("\n========================================================================\n")
cat("  TENSOR PRODUCT INTERACTIONS (age × key biomarkers)\n")
cat("========================================================================\n\n")

# For the most interesting biomarkers, fit age interactions
# using tensor product smooths ti()
cat("  Fitting model with age×albumin and age×glucose interactions...\n")

t0 <- proc.time()
gam_interact <- gam(
  time_yr ~
    s(albumin_gL, k = 10) +
    s(creat_umol, k = 10) +
    s(glucose_mmol, k = 10) +
    s(lncrp, k = 10) +
    s(lymph, k = 10) +
    s(mcv, k = 10) +
    s(rdw, k = 10) +
    s(alp, k = 10) +
    s(wbc, k = 10) +
    s(age, k = 20) +
    ti(age, albumin_gL, k = c(5, 5)) +
    ti(age, glucose_mmol, k = c(5, 5)) +
    ti(age, lncrp, k = c(5, 5)) +
    sex,
  family = cox.ph(),
  data = train,
  weights = event,
  method = "REML"
)
t1 <- proc.time()
cat("  Interaction model fit in", round((t1 - t0)[3], 1), "seconds\n\n")

cat("  Interaction model summary:\n")
int_summary <- summary(gam_interact)
int_edf <- int_summary$s.table
cat(sprintf("  %-35s  %8s  %10s  %10s\n", "Term", "edf", "Chi.sq", "p-value"))
cat(sprintf("  %-35s  %8s  %10s  %10s\n", "----", "---", "------", "-------"))
for (i in seq_len(nrow(int_edf))) {
  cat(sprintf("  %-35s  %8.2f  %10.1f  %10.1e\n",
              rownames(int_edf)[i], int_edf[i, "edf"],
              int_edf[i, "Chi.sq"], int_edf[i, "p-value"]))
}

# Compare AIC
cat(sprintf("\n  AIC comparison:\n"))
cat(sprintf("    Linear Cox PH:      %.1f\n", AIC(cox_linear)))
cat(sprintf("    GAM (main effects): %.1f\n", AIC(gam_full)))
cat(sprintf("    GAM (interactions): %.1f\n", AIC(gam_interact)))

# Interaction model C-index
int_lp <- predict(gam_interact, type = "link")
int_conc <- concordance(surv_obj ~ int_lp)
int_cindex <- int_conc$concordance
if (int_cindex < 0.5) int_cindex <- 1 - int_cindex
cat(sprintf("    Interaction GAM C-index: %.4f\n", int_cindex))

# ──────────────────────────────────────────────────────────────────────
# 11. GENERATE ALTAIR-COMPATIBLE JSON FOR SMOOTH PLOTS
# ──────────────────────────────────────────────────────────────────────

# Export predictions from GAM for Altair visualization in Python
cat("\n========================================================================\n")
cat("  EXPORTING PREDICTIONS FOR VISUALIZATION\n")
cat("========================================================================\n\n")

# For each biomarker: sweep its value, compute GAM predicted log-hazard
# and linear Cox predicted log-hazard, for comparison
sweep_results <- list()

for (feat in features) {
  feat_vals <- seq(
    quantile(train[[feat]], 0.01, na.rm = TRUE),
    quantile(train[[feat]], 0.99, na.rm = TRUE),
    length.out = 200
  )

  # Create newdata at medians
  newdata <- as.data.frame(lapply(train[, c(features, "sex", "time_yr")],
                                   median, na.rm = TRUE))
  newdata <- newdata[rep(1, 200), ]
  newdata[[feat]] <- feat_vals

  # GAM predictions (log hazard ratio relative to baseline)
  gam_pred <- predict(gam_full, newdata = newdata, type = "link", se.fit = TRUE)

  # Linear Cox predictions
  lin_pred <- predict(cox_linear, newdata = newdata, type = "lp")

  sweep_results[[feat]] <- data.frame(
    feature = feat,
    value = feat_vals,
    gam_log_hr = gam_pred$fit,
    gam_se = gam_pred$se.fit,
    linear_log_hr = lin_pred
  )
}

sweep_df <- do.call(rbind, sweep_results)
sweep_csv <- file.path(outdir, "gam_vs_linear_sweeps.csv")
write.csv(sweep_df, sweep_csv, row.names = FALSE)
cat("  Saved sweep comparison data:", sweep_csv, "\n")

# ──────────────────────────────────────────────────────────────────────
# 12. FINAL SUMMARY TABLE
# ──────────────────────────────────────────────────────────────────────

cat("\n")
cat("========================================================================\n")
cat("  FINAL SUMMARY\n")
cat("========================================================================\n\n")

cat(sprintf("  %-35s  %8s  %8s\n", "Model", "C-index", "AIC"))
cat(sprintf("  %-35s  %8s  %8s\n", "-----", "-------", "---"))
cat(sprintf("  %-35s  %8.4f  %8.1f\n", "Linear Cox PH",
            lin_conc["concordance"], AIC(cox_linear)))
cat(sprintf("  %-35s  %8.4f  %8.1f\n", "GAM Cox PH (main effects)",
            gam_cindex, AIC(gam_full)))
cat(sprintf("  %-35s  %8.4f  %8.1f\n", "GAM Cox PH (interactions)",
            int_cindex, AIC(gam_interact)))
cat(sprintf("  %-35s  %8.4f  %8.1f\n", "Penalized Cox PH (pspline)",
            pspline_conc["concordance"], AIC(cox_pspline)))

cat(sprintf("\n  CV 5-fold C-index:\n"))
cat(sprintf("    GAM:    %.4f ± %.4f\n",
            mean(cv_cindex_gam, na.rm = TRUE), sd(cv_cindex_gam, na.rm = TRUE)))
cat(sprintf("    Linear: %.4f ± %.4f\n",
            mean(cv_cindex_lin, na.rm = TRUE), sd(cv_cindex_lin, na.rm = TRUE)))

cat("\n  Key findings — non-linear biomarkers (edf > 2):\n")
nonlinear <- edf_table[edf_table$edf > 2, ]
if (nrow(nonlinear) > 0) {
  for (i in seq_len(nrow(nonlinear))) {
    cat(sprintf("    %s (edf=%.1f)\n", nonlinear$term[i], nonlinear$edf[i]))
  }
} else {
  cat("    (none — all relationships approximately linear)\n")
}

cat("\n  Interpretation:\n")
cat("    - edf ≈ 1: relationship is linear (GAM agrees with PhenoAge)\n")
cat("    - edf > 2: relationship is non-linear (GAM captures curvature)\n")
cat("    - Smooth plots show the actual shape of each relationship\n")
cat("    - U-shaped curves mean BOTH very low AND very high values\n")
cat("      are associated with increased mortality risk\n")

cat("\n  Output files:\n")
cat(sprintf("    %s\n", pdf_path))
cat(sprintf("    %s\n", csv_path))
cat(sprintf("    %s\n", sweep_csv))

cat("\n========================================================================\n")
cat("  DONE\n")
cat("========================================================================\n")
