#!/usr/bin/env Rscript
# Compare R ranger and randomForestSRC with Python scikit-survival RSF.
#
# Usage: Rscript src/phenoage/rsf_r_comparison.R

library(survival)
library(ranger)
library(randomForestSRC)

cat("========================================================================\n")
cat("  RSF MODEL COMPARISON IN R (ranger + randomForestSRC)\n")
cat("========================================================================\n\n")

# ── Load training data ──
train <- read.csv("research_agents/phenoage_train_2015.csv")
cat("Training data:", nrow(train), "rows,", sum(train$event), "deaths\n")

# Features
features <- c("albumin_gL", "creat_umol", "glucose_mmol", "lncrp",
               "lymph", "mcv", "rdw", "alp", "wbc", "age")

# Add sex (HSSEX: 1=Male, 2=Female -> binary)
if ("HSSEX" %in% colnames(train)) {
  train$sex <- ifelse(train$HSSEX == 2, 1, 0)
  features <- c(features, "sex")
}

cat("Features:", paste(features, collapse=", "), "\n\n")

# ── Prepare survival formula ──
surv_formula <- as.formula(paste(
  "Surv(time, event) ~", paste(features, collapse=" + ")
))

# ── 1. ranger ──
cat("Fitting ranger RSF (500 trees, min.node.size=15)...\n")
t0 <- proc.time()
rf_ranger <- ranger(
  surv_formula,
  data=train,
  num.trees=500,
  min.node.size=15,
  importance="permutation",
  seed=42,
  num.threads=4
)
t1 <- proc.time()
cat("  ranger done in", round((t1 - t0)[3], 1), "seconds\n")

# C-index on training data
cat("  ranger C-index (OOB):", round(1 - rf_ranger$prediction.error, 4), "\n")

# Variable importance
cat("\n  ranger variable importance (permutation):\n")
vimp <- sort(rf_ranger$variable.importance, decreasing=TRUE)
for (v in names(vimp)) {
  cat(sprintf("    %-15s: %.4f\n", v, vimp[v]))
}

# ── 2. randomForestSRC ──
cat("\nFitting randomForestSRC (500 trees, nodesize=15)...\n")
t0 <- proc.time()
rf_rfsrc <- rfsrc(
  surv_formula,
  data=train,
  ntree=500,
  nodesize=15,
  importance=TRUE,
  seed=42
)
t1 <- proc.time()
cat("  rfsrc done in", round((t1 - t0)[3], 1), "seconds\n")

# C-index
rfsrc_err <- rf_rfsrc$err.rate[length(rf_rfsrc$err.rate)]
cat("  rfsrc C-index (OOB):", round(1 - rfsrc_err, 4), "\n")

# Variable importance
cat("\n  rfsrc variable importance:\n")
vimp2 <- sort(rf_rfsrc$importance, decreasing=TRUE)
for (v in names(vimp2)) {
  cat(sprintf("    %-15s: %.4f\n", v, vimp2[v]))
}

# ── Summary ──
cat("\n========================================================================\n")
cat("  SUMMARY\n")
cat("========================================================================\n")
cat(sprintf("  %-25s C-index\n", "Model"))
cat(sprintf("  %-25s -------\n", ""))
cat(sprintf("  %-25s %.4f (OOB)\n", "ranger", 1 - rf_ranger$prediction.error))
cat(sprintf("  %-25s %.4f (OOB)\n", "randomForestSRC", 1 - rfsrc_err))
cat(sprintf("  %-25s %.4f (test, from Python)\n", "Python sksurv RSF", 0.8512))
cat(sprintf("  %-25s %.4f (test, from Python)\n", "Linear PhenoAge", 0.8551))
