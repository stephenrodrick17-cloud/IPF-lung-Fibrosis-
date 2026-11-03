#!/usr/bin/env Rscript
## ==============================================================================
## 10_prespecified_classifier.R
## PART 3 & PART 4: Multi-Gene Classifier & Feature Importance (15 Pre-Specified Genes)
## Seed: 20261102 | R version: 4.6.2 (Windows)
## ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0: Pre-Specified Rules & Parameters
# ------------------------------------------------------------------------------
# Models:
# 1. Model 1: L2-regularized logistic regression (Ridge, glmnet alpha=0, family="binomial")
# 2. Model 2: Unweighted directional score (sum of z-scores * sign(logFC))
# Validation:
# - Repeated Stratified 5-Fold Cross-Validation x 20 repeats (100 folds) on GSE32539 (N=168)
# - External Evaluation: GSE110147 (N=33), GSE24206 (N=15), GSE135065 (N=18, independent)
# Baseline Comparisons (1,000 random sets each):
# - 1,000 random 15-gene sets from 8,000 network genes
# - 1,000 random 15-gene sets from 639 matrisome genes
# - 1,000 random 15-gene sets from 45 shared DEGs
# - 200 label permutations
# Leave-One-Gene-Out Analysis: Delta AUC per gene
# Feature Importance: Random Forest (500 trees) & Linear SVM weights (no gene pruning)
# ------------------------------------------------------------------------------

set.seed(20261102)
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(randomForest)
  library(e1071)
  library(ggplot2)
  library(patchwork)
})

dirs <- c("logs", "results/tables/prespecified", "results/figures/prespecified")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- "logs/pre_specified_test.log"
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("\n=== Initializing PART 3 & 4: Multi-Gene Classifier Testing ===")

save_plot <- function(p, file_base, width = 9, height = 6) {
  ggsave(paste0(file_base, ".png"), plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

prespecified_15 <- c(
  "ASPN", "COL14A1", "COL15A1", "COL1A1", "COL3A1", "COL5A2", "COMP",
  "IL6", "MMP7", "MMP8", "MMP10", "OGN", "SFRP4", "TGFB3", "TIMP3"
)

# ------------------------------------------------------------------------------
# STEP 1: Load and Prepare Datasets
# ------------------------------------------------------------------------------
datasets <- c("GSE32539", "GSE110147", "GSE24206", "GSE135065")
data_list <- list()

for (gse in datasets) {
  e <- read.csv(paste0("data/processed/", gse, "_expr_gene_level.csv"), check.names = FALSE)
  rownames(e) <- e$gene
  e$gene <- NULL
  
  meta <- read.csv(paste0("data/processed/", gse, "_metadata.csv"))
  if ("group" %in% colnames(meta)) {
    grp <- meta$group
  } else if ("disease_state" %in% colnames(meta)) {
    grp <- meta$disease_state
  } else {
    grp <- meta$label
  }
  names(grp) <- meta$sample_id
  
  if (gse == "GSE32539") {
    valid_samples <- intersect(colnames(e), setdiff(meta$sample_id, "GSM806411"))
  } else {
    valid_samples <- intersect(colnames(e), meta$sample_id)
  }
  
  grp_sub <- grp[valid_samples]
  keep <- grp_sub %in% c("IPF", "Control")
  valid_samples <- valid_samples[keep]
  grp_sub <- factor(grp_sub[keep], levels = c("Control", "IPF"))
  
  # Standardize (z-score) within dataset across valid samples
  e_sub <- as.matrix(e[prespecified_15, valid_samples])
  z_sub <- t(scale(t(e_sub)))
  
  data_list[[gse]] <- list(
    expr_raw = e_sub,
    expr_z = z_sub,
    y = ifelse(grp_sub == "IPF", 1, 0),
    group = grp_sub,
    meta = meta[meta$sample_id %in% valid_samples, ]
  )
}

log_msg("Data loaded and standardized across all 4 cohorts.")

# Training dataset: GSE32539
X_train_z <- t(data_list[["GSE32539"]]$expr_z) # samples in rows
y_train <- data_list[["GSE32539"]]$y

# Weights for Model 2: sign of logFC in GSE32539
lim32 <- read.csv("results/tables/GSE32539_limma_full.csv")
rownames(lim32) <- lim32$gene
weights_m2 <- setNames(sign(lim32[prespecified_15, "logFC"]), prespecified_15)
log_msg("Model 2 weights (sign of logFC): ", paste(paste(names(weights_m2), weights_m2, sep = "="), collapse = ", "))

# Function to compute classification metrics
get_metrics <- function(y_true, y_prob, y_score = NULL) {
  # If y_score is supplied (e.g. unweighted score), evaluate ROC on score
  pred_val <- if (!is.null(y_score)) y_score else y_prob
  
  roc_res <- suppressMessages(pROC::roc(response = y_true, predictor = pred_val, levels = c(0, 1), direction = "<"))
  auc_val <- as.numeric(pROC::auc(roc_res))
  ci_val <- suppressMessages(pROC::ci.auc(roc_res, method = "delong"))
  
  # Optimal threshold by Youden index
  best_coords <- suppressMessages(pROC::coords(roc_res, "best", ret = c("threshold", "sensitivity", "specificity", "accuracy"), best.method = "youden"))
  if (is.matrix(best_coords) || is.data.frame(best_coords)) best_coords <- best_coords[1, ]
  
  sens <- as.numeric(best_coords["sensitivity"])
  spec <- as.numeric(best_coords["specificity"])
  bal_acc <- (sens + spec) / 2
  
  return(list(
    roc = roc_res,
    AUC = auc_val,
    AUC_CI_lower = as.numeric(ci_val[1]),
    AUC_CI_upper = as.numeric(ci_val[3]),
    sensitivity = sens,
    specificity = spec,
    balanced_accuracy = bal_acc
  ))
}

# ------------------------------------------------------------------------------
# STEP 2: Repeated Stratified 5-Fold Cross-Validation on GSE32539
# ------------------------------------------------------------------------------
log_msg("\n--- Running Repeated Stratified 5-Fold CV (20 repeats = 100 folds) on GSE32539 ---")

n_repeats <- 20
k_folds <- 5
n_samples <- nrow(X_train_z)

# Store out-of-fold predictions
oof_prob_m1 <- matrix(NA, nrow = n_samples, ncol = n_repeats)
oof_score_m2 <- matrix(NA, nrow = n_samples, ncol = n_repeats)

ipf_indices <- which(y_train == 1)
ctrl_indices <- which(y_train == 0)

for (r in 1:n_repeats) {
  set.seed(20261102 + r)
  folds_ipf <- sample(rep(1:k_folds, length.out = length(ipf_indices)))
  folds_ctrl <- sample(rep(1:k_folds, length.out = length(ctrl_indices)))
  
  fold_assign <- numeric(n_samples)
  fold_assign[ipf_indices] <- folds_ipf
  fold_assign[ctrl_indices] <- folds_ctrl
  
  for (k in 1:k_folds) {
    val_idx <- which(fold_assign == k)
    train_idx <- which(fold_assign != k)
    
    # Standardize fold-specifically to prevent data leakage
    raw_train <- t(data_list[["GSE32539"]]$expr_raw[, train_idx])
    raw_val <- t(data_list[["GSE32539"]]$expr_raw[, val_idx])
    
    train_mean <- colMeans(raw_train)
    train_sd <- apply(raw_train, 2, sd)
    
    z_tr <- scale(raw_train, center = train_mean, scale = train_sd)
    z_vl <- scale(raw_val, center = train_mean, scale = train_sd)
    
    y_tr <- y_train[train_idx]
    
    # Model 1: Ridge Logistic Regression
    cv_ridge <- cv.glmnet(z_tr, y_tr, alpha = 0, family = "binomial", type.measure = "auc", nfolds = 5)
    prob_vl <- as.numeric(predict(cv_ridge, newx = z_vl, type = "response", s = "lambda.min"))
    oof_prob_m1[val_idx, r] <- prob_vl
    
    # Model 2: Unweighted Score using training fold signs
    # Estimate train sign(diff)
    diff_train <- colMeans(raw_train[y_tr == 1, , drop = FALSE]) - colMeans(raw_train[y_tr == 0, , drop = FALSE])
    w_fold <- sign(diff_train)
    score_vl <- as.numeric(z_vl %*% w_fold)
    oof_score_m2[val_idx, r] <- score_vl
  }
}

mean_oof_prob_m1 <- rowMeans(oof_prob_m1)
mean_oof_score_m2 <- rowMeans(oof_score_m2)

cv_metrics_m1 <- get_metrics(y_train, mean_oof_prob_m1)
cv_metrics_m2 <- get_metrics(y_train, y_score = mean_oof_score_m2)

log_msg(sprintf("GSE32539 20x5-fold CV - Model 1 (Ridge): AUC = %.4f (95%% CI: %.4f - %.4f), BalAcc = %.4f",
                cv_metrics_m1$AUC, cv_metrics_m1$AUC_CI_lower, cv_metrics_m1$AUC_CI_upper, cv_metrics_m1$balanced_accuracy))
log_msg(sprintf("GSE32539 20x5-fold CV - Model 2 (Unweighted): AUC = %.4f (95%% CI: %.4f - %.4f), BalAcc = %.4f",
                cv_metrics_m2$AUC, cv_metrics_m2$AUC_CI_lower, cv_metrics_m2$AUC_CI_upper, cv_metrics_m2$balanced_accuracy))

# ------------------------------------------------------------------------------
# STEP 3: External Evaluation across Datasets
# ------------------------------------------------------------------------------
log_msg("\n--- Training Full Models on GSE32539 and Evaluating Externally ---")

# Full Model 1 fit on all GSE32539
full_ridge <- cv.glmnet(X_train_z, y_train, alpha = 0, family = "binomial", type.measure = "auc", nfolds = 5)
ridge_coefs <- coef(full_ridge, s = "lambda.min")

# Full Model 2 score function
compute_m2_score <- function(z_matrix) {
  as.numeric(z_matrix %*% weights_m2)
}

perf_results <- list()

# Record GSE32539 CV
perf_results[["GSE32539_CV_Model1"]] <- data.frame(
  dataset = "GSE32539 (CV)",
  cohort_type = "Internal_CV",
  model = "Model 1: Ridge Logistic",
  n_samples = n_samples,
  n_IPF = sum(y_train == 1),
  n_Control = sum(y_train == 0),
  AUC = round(cv_metrics_m1$AUC, 4),
  AUC_CI_lower = round(cv_metrics_m1$AUC_CI_lower, 4),
  AUC_CI_upper = round(cv_metrics_m1$AUC_CI_upper, 4),
  Balanced_Accuracy = round(cv_metrics_m1$balanced_accuracy, 4),
  Sensitivity = round(cv_metrics_m1$sensitivity, 4),
  Specificity = round(cv_metrics_m1$specificity, 4)
)

perf_results[["GSE32539_CV_Model2"]] <- data.frame(
  dataset = "GSE32539 (CV)",
  cohort_type = "Internal_CV",
  model = "Model 2: Unweighted Directional Score",
  n_samples = n_samples,
  n_IPF = sum(y_train == 1),
  n_Control = sum(y_train == 0),
  AUC = round(cv_metrics_m2$AUC, 4),
  AUC_CI_lower = round(cv_metrics_m2$AUC_CI_lower, 4),
  AUC_CI_upper = round(cv_metrics_m2$AUC_CI_upper, 4),
  Balanced_Accuracy = round(cv_metrics_m2$balanced_accuracy, 4),
  Sensitivity = round(cv_metrics_m2$sensitivity, 4),
  Specificity = round(cv_metrics_m2$specificity, 4)
)

ext_cohorts <- c("GSE110147", "GSE24206", "GSE135065")
ext_roc_list <- list()

for (gse in ext_cohorts) {
  d <- data_list[[gse]]
  X_ext_z <- t(d$expr_z)
  y_ext <- d$y
  
  # Model 1 external predictions
  prob_m1 <- as.numeric(predict(full_ridge, newx = X_ext_z, type = "response", s = "lambda.min"))
  met_m1 <- get_metrics(y_ext, prob_m1)
  
  # Model 2 external predictions
  score_m2 <- compute_m2_score(X_ext_z)
  met_m2 <- get_metrics(y_ext, y_score = score_m2)
  
  ext_roc_list[[paste(gse, "Model1", sep = "_")]] <- met_m1$roc
  ext_roc_list[[paste(gse, "Model2", sep = "_")]] <- met_m2$roc
  
  is_indep <- (gse == "GSE135065")
  cohort_desc <- ifelse(is_indep, "Independent_Validation", "External_Discovery_Cohort")
  
  perf_results[[paste(gse, "Model1", sep = "_")]] <- data.frame(
    dataset = gse,
    cohort_type = cohort_desc,
    model = "Model 1: Ridge Logistic",
    n_samples = length(y_ext),
    n_IPF = sum(y_ext == 1),
    n_Control = sum(y_ext == 0),
    AUC = round(met_m1$AUC, 4),
    AUC_CI_lower = round(met_m1$AUC_CI_lower, 4),
    AUC_CI_upper = round(met_m1$AUC_CI_upper, 4),
    Balanced_Accuracy = round(met_m1$balanced_accuracy, 4),
    Sensitivity = round(met_m1$sensitivity, 4),
    Specificity = round(met_m1$specificity, 4)
  )
  
  perf_results[[paste(gse, "Model2", sep = "_")]] <- data.frame(
    dataset = gse,
    cohort_type = cohort_desc,
    model = "Model 2: Unweighted Directional Score",
    n_samples = length(y_ext),
    n_IPF = sum(y_ext == 1),
    n_Control = sum(y_ext == 0),
    AUC = round(met_m2$AUC, 4),
    AUC_CI_lower = round(met_m2$AUC_CI_lower, 4),
    AUC_CI_upper = round(met_m2$AUC_CI_upper, 4),
    Balanced_Accuracy = round(met_m2$balanced_accuracy, 4),
    Sensitivity = round(met_m2$sensitivity, 4),
    Specificity = round(met_m2$specificity, 4)
  )
  
  log_msg(sprintf("%s (%s) - Model 1 AUC: %.4f (95%% CI: %.4f - %.4f), BalAcc: %.4f",
                  gse, cohort_desc, met_m1$AUC, met_m1$AUC_CI_lower, met_m1$AUC_CI_upper, met_m1$balanced_accuracy))
  log_msg(sprintf("%s (%s) - Model 2 AUC: %.4f (95%% CI: %.4f - %.4f), BalAcc: %.4f",
                  gse, cohort_desc, met_m2$AUC, met_m2$AUC_CI_lower, met_m2$AUC_CI_upper, met_m2$balanced_accuracy))
}

# Sub-analysis: GSE135065 Untreated Baseline (treatment: none, N=6: 3 IPF vs 3 Control)
base_meta <- data_list[["GSE135065"]]$meta
base_idx <- which(base_meta$treatment == "treatment: none")
if (length(base_idx) == 6) {
  X_base_z <- t(data_list[["GSE135065"]]$expr_z[, base_idx])
  y_base <- data_list[["GSE135065"]]$y[base_idx]
  
  prob_base_m1 <- as.numeric(predict(full_ridge, newx = X_base_z, type = "response", s = "lambda.min"))
  met_base_m1 <- get_metrics(y_base, prob_base_m1)
  
  score_base_m2 <- compute_m2_score(X_base_z)
  met_base_m2 <- get_metrics(y_base, y_score = score_base_m2)
  
  perf_results[["GSE135065_Baseline_Model1"]] <- data.frame(
    dataset = "GSE135065 (Baseline untreated)",
    cohort_type = "Independent_Subanalysis",
    model = "Model 1: Ridge Logistic",
    n_samples = 6,
    n_IPF = 3,
    n_Control = 3,
    AUC = round(met_base_m1$AUC, 4),
    AUC_CI_lower = round(met_base_m1$AUC_CI_lower, 4),
    AUC_CI_upper = round(met_base_m1$AUC_CI_upper, 4),
    Balanced_Accuracy = round(met_base_m1$balanced_accuracy, 4),
    Sensitivity = round(met_base_m1$sensitivity, 4),
    Specificity = round(met_base_m1$specificity, 4)
  )
  
  perf_results[["GSE135065_Baseline_Model2"]] <- data.frame(
    dataset = "GSE135065 (Baseline untreated)",
    cohort_type = "Independent_Subanalysis",
    model = "Model 2: Unweighted Directional Score",
    n_samples = 6,
    n_IPF = 3,
    n_Control = 3,
    AUC = round(met_base_m2$AUC, 4),
    AUC_CI_lower = round(met_base_m2$AUC_CI_lower, 4),
    AUC_CI_upper = round(met_base_m2$AUC_CI_upper, 4),
    Balanced_Accuracy = round(met_base_m2$balanced_accuracy, 4),
    Sensitivity = round(met_base_m2$sensitivity, 4),
    Specificity = round(met_base_m2$specificity, 4)
  )
  
  log_msg(sprintf("GSE135065 Untreated Baseline (N=6) - Model 1 AUC: %.4f, Model 2 AUC: %.4f",
                  met_base_m1$AUC, met_base_m2$AUC))
}

perf_df <- do.call(rbind, perf_results)
write.csv(perf_df, "results/tables/prespecified/classifier_performance.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/classifier_performance.csv")

# ------------------------------------------------------------------------------
# STEP 4: Random-Set Comparisons (1,000 sets each) & Label Permutations
# ------------------------------------------------------------------------------
log_msg("\n--- Running Random-Set Comparisons (1,000 sets each) and 200 Permutations ---")

# Universe sets
net_genes <- read.csv("results/tables/wgcna/module_assignment.csv")$gene
mat_df <- read.csv("results/tables/GSE32539_matrisome_all_genes.csv", check.names = FALSE)
mat_genes <- intersect(mat_df[["Annotated Gene"]], rownames(lim32))
deg45_genes <- intersect(read.csv("results/tables/common_matrisome_DEGs_preview.csv")$gene, rownames(lim32))

# Load full GSE32539 matrix for random sampling
e32_full <- read.csv("data/processed/GSE32539_expr_gene_level.csv", check.names = FALSE)
rownames(e32_full) <- e32_full$gene
e32_full$gene <- NULL
keep32 <- setdiff(colnames(data_list[["GSE32539"]]$expr_raw), "GSM806411")
e32_mat <- as.matrix(e32_full[, colnames(data_list[["GSE32539"]]$expr_raw)])

# Function to evaluate directional score AUC for any 15-gene set on GSE32539
eval_set_auc <- function(gene_set) {
  sub_e <- e32_mat[gene_set, , drop = FALSE]
  z_e <- t(scale(t(sub_e)))
  w <- sign(lim32[gene_set, "logFC"])
  score <- as.numeric(t(z_e) %*% w)
  roc_s <- suppressMessages(pROC::roc(response = y_train, predictor = score, levels = c(0, 1), direction = "<"))
  return(as.numeric(pROC::auc(roc_s)))
}

N_PERM <- 1000
set.seed(20261102)

log_msg("Sampling 1,000 random sets of 15 genes from Network (8,000 genes)...")
auc_null_network <- numeric(N_PERM)
for (i in 1:N_PERM) {
  g_samp <- sample(net_genes, 15)
  auc_null_network[i] <- eval_set_auc(g_samp)
}

log_msg("Sampling 1,000 random sets of 15 genes from Matrisome (639 genes)...")
auc_null_matrisome <- numeric(N_PERM)
for (i in 1:N_PERM) {
  g_samp <- sample(mat_genes, 15)
  auc_null_matrisome[i] <- eval_set_auc(g_samp)
}

log_msg("Sampling 1,000 random sets of 15 genes from Shared 45 DEGs...")
auc_null_deg45 <- numeric(N_PERM)
for (i in 1:N_PERM) {
  g_samp <- sample(deg45_genes, 15)
  auc_null_deg45[i] <- eval_set_auc(g_samp)
}

# Label Permutation test (200 permutations)
log_msg("Running 200 label permutations on GSE32539 with the 15-gene signature...")
auc_label_perm <- numeric(200)
for (p in 1:200) {
  y_perm <- sample(y_train)
  score_15 <- as.numeric(X_train_z %*% weights_m2)
  roc_p <- suppressMessages(pROC::roc(response = y_perm, predictor = score_15, levels = c(0, 1), direction = "<"))
  auc_label_perm[p] <- as.numeric(pROC::auc(roc_p))
}

# Observed 15-gene score AUC on GSE32539
obs_auc_15 <- eval_set_auc(prespecified_15)

# Empirical p-values
p_emp_net <- (1 + sum(auc_null_network >= obs_auc_15)) / (1 + N_PERM)
p_emp_mat <- (1 + sum(auc_null_matrisome >= obs_auc_15)) / (1 + N_PERM)
p_emp_deg45 <- (1 + sum(auc_null_deg45 >= obs_auc_15)) / (1 + N_PERM)
p_emp_perm <- (1 + sum(auc_label_perm >= obs_auc_15)) / (1 + 200)

log_msg(sprintf("Observed 15-gene Score AUC: %.4f", obs_auc_15))
log_msg(sprintf("Empirical p vs 1,000 Network random sets: p = %.4f (mean null AUC: %.4f)", p_emp_net, mean(auc_null_network)))
log_msg(sprintf("Empirical p vs 1,000 Matrisome random sets: p = %.4f (mean null AUC: %.4f)", p_emp_mat, mean(auc_null_matrisome)))
log_msg(sprintf("Empirical p vs 1,000 45-DEG random sets: p = %.4f (mean null AUC: %.4f)", p_emp_deg45, mean(auc_null_deg45)))
log_msg(sprintf("Empirical p vs 200 Label permutations: p = %.4f (mean null AUC: %.4f)", p_emp_perm, mean(auc_label_perm)))

random_summary <- data.frame(
  comparison_type = c("Network_Random_15", "Matrisome_Random_15", "Shared_45_DEG_Random_15", "Label_Permutations"),
  n_iterations = c(N_PERM, N_PERM, N_PERM, 200),
  observed_15gene_AUC = round(obs_auc_15, 4),
  mean_null_AUC = round(c(mean(auc_null_network), mean(auc_null_matrisome), mean(auc_null_deg45), mean(auc_label_perm)), 4),
  sd_null_AUC = round(c(sd(auc_null_network), sd(auc_null_matrisome), sd(auc_null_deg45), sd(auc_label_perm)), 4),
  q95_null_AUC = round(c(quantile(auc_null_network, 0.95), quantile(auc_null_matrisome, 0.95), quantile(auc_null_deg45, 0.95), quantile(auc_label_perm, 0.95)), 4),
  empirical_p_value = c(p_emp_net, p_emp_mat, p_emp_deg45, p_emp_perm)
)

write.csv(random_summary, "results/tables/prespecified/random_set_comparison.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/random_set_comparison.csv")

# ------------------------------------------------------------------------------
# STEP 5: Leave-One-Gene-Out Sensitivity Analysis
# ------------------------------------------------------------------------------
log_msg("\n--- Running Leave-One-Gene-Out Analysis ---")

loo_list <- list()

# Baseline AUCs with all 15 genes
auc_base_32 <- cv_metrics_m2$AUC
auc_base_110 <- get_metrics(data_list[["GSE110147"]]$y, y_score = compute_m2_score(t(data_list[["GSE110147"]]$expr_z)))$AUC
auc_base_24 <- get_metrics(data_list[["GSE24206"]]$y, y_score = compute_m2_score(t(data_list[["GSE24206"]]$expr_z)))$AUC
auc_base_135 <- get_metrics(data_list[["GSE135065"]]$y, y_score = compute_m2_score(t(data_list[["GSE135065"]]$expr_z)))$AUC

for (omitted in prespecified_15) {
  remain_genes <- setdiff(prespecified_15, omitted)
  w_remain <- weights_m2[remain_genes]
  
  # Helper to compute score with remaining genes
  eval_remain <- function(gse) {
    z_mat <- t(data_list[[gse]]$expr_z[remain_genes, , drop = FALSE])
    sc <- as.numeric(z_mat %*% w_remain)
    roc_r <- suppressMessages(pROC::roc(response = data_list[[gse]]$y, predictor = sc, levels = c(0, 1), direction = "<"))
    return(as.numeric(pROC::auc(roc_r)))
  }
  
  auc_loo_32 <- eval_remain("GSE32539")
  auc_loo_110 <- eval_remain("GSE110147")
  auc_loo_24 <- eval_remain("GSE24206")
  auc_loo_135 <- eval_remain("GSE135065")
  
  loo_list[[omitted]] <- data.frame(
    omitted_gene = omitted,
    AUC_GSE32539_all15 = round(auc_base_32, 4),
    AUC_GSE32539_minus_gene = round(auc_loo_32, 4),
    delta_AUC_GSE32539 = round(auc_loo_32 - auc_base_32, 4),
    AUC_GSE110147_minus_gene = round(auc_loo_110, 4),
    delta_AUC_GSE110147 = round(auc_loo_110 - auc_base_110, 4),
    AUC_GSE24206_minus_gene = round(auc_loo_24, 4),
    delta_AUC_GSE24206 = round(auc_loo_24 - auc_base_24, 4),
    AUC_GSE135065_minus_gene = round(auc_loo_135, 4),
    delta_AUC_GSE135065 = round(auc_loo_135 - auc_base_135, 4)
  )
}

loo_df <- do.call(rbind, loo_list)
write.csv(loo_df, "results/tables/prespecified/leave_one_out.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/leave_one_out.csv")

# ------------------------------------------------------------------------------
# STEP 6: PART 4 - Feature Importance Ranking (Random Forest & Linear SVM)
# ------------------------------------------------------------------------------
log_msg("\n--- Computing Feature Importance Rankings (RF and Linear SVM) on GSE32539 ---")

# Random Forest (500 trees, importance = TRUE)
df_rf <- as.data.frame(X_train_z)
df_rf$label <- factor(ifelse(y_train == 1, "IPF", "Control"), levels = c("Control", "IPF"))

set.seed(20261102)
rf_fit <- randomForest(label ~ ., data = df_rf, ntree = 500, importance = TRUE)
rf_imp <- importance(rf_fit)

rf_mda <- rf_imp[, "MeanDecreaseAccuracy"]
rf_mdg <- rf_imp[, "MeanDecreaseGini"]
rf_rank_mda <- rank(-rf_mda)
rf_rank_mdg <- rank(-rf_mdg)

# Linear SVM Feature Weights
# For linear SVM, w = sum(coef * support_vectors)
svm_fit <- svm(label ~ ., data = df_rf, kernel = "linear", scale = FALSE, cost = 1)
w_svm <- t(svm_fit$coefs) %*% svm_fit$SV
w_svm_abs <- abs(as.numeric(w_svm))
names(w_svm_abs) <- colnames(df_rf)[1:15]
svm_rank <- rank(-w_svm_abs)

# Ridge absolute regression coefficients
ridge_beta <- as.matrix(ridge_coefs)[-1, 1]
ridge_rank <- rank(-abs(ridge_beta))

feat_imp_df <- data.frame(
  gene = prespecified_15,
  RF_MeanDecreaseAccuracy = round(rf_mda[prespecified_15], 4),
  RF_MDA_rank = as.integer(rf_rank_mda[prespecified_15]),
  RF_MeanDecreaseGini = round(rf_mdg[prespecified_15], 4),
  RF_MDG_rank = as.integer(rf_rank_mdg[prespecified_15]),
  SVM_Linear_Weight = round(w_svm_abs[prespecified_15], 4),
  SVM_Rank = as.integer(svm_rank[prespecified_15]),
  Ridge_Beta = round(ridge_beta[prespecified_15], 4),
  Ridge_Rank = as.integer(ridge_rank[prespecified_15])
)

# Overall average rank across models
feat_imp_df$Mean_Importance_Rank <- round(
  rowMeans(feat_imp_df[, c("RF_MDA_rank", "RF_MDG_rank", "SVM_Rank", "Ridge_Rank")]), 2
)
feat_imp_df <- feat_imp_df[order(feat_imp_df$Mean_Importance_Rank), ]

write.csv(feat_imp_df, "results/tables/prespecified/feature_importance_15.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/feature_importance_15.csv")
log_msg("Top 5 genes by combined ML importance: ", paste(feat_imp_df$gene[1:5], collapse = ", "))

# ------------------------------------------------------------------------------
# STEP 7: Visualizations
# ------------------------------------------------------------------------------
log_msg("\n--- Generating PART 3 & PART 4 Figures ---")

# 7.1 Multi-Cohort ROC Curves
roc_df_list <- list()

add_roc_curve <- function(roc_obj, cohort_label, model_label) {
  data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Cohort = cohort_label,
    Model = model_label,
    AUC = round(as.numeric(pROC::auc(roc_obj)), 3)
  )
}

roc_df_list[[1]] <- add_roc_curve(cv_metrics_m1$roc, "GSE32539 (CV)", "Ridge Logistic")
roc_df_list[[2]] <- add_roc_curve(cv_metrics_m2$roc, "GSE32539 (CV)", "Unweighted Score")
roc_df_list[[3]] <- add_roc_curve(ext_roc_list[["GSE110147_Model1"]], "GSE110147", "Ridge Logistic")
roc_df_list[[4]] <- add_roc_curve(ext_roc_list[["GSE110147_Model2"]], "GSE110147", "Unweighted Score")
roc_df_list[[5]] <- add_roc_curve(ext_roc_list[["GSE24206_Model1"]], "GSE24206", "Ridge Logistic")
roc_df_list[[6]] <- add_roc_curve(ext_roc_list[["GSE24206_Model2"]], "GSE24206", "Unweighted Score")
roc_df_list[[7]] <- add_roc_curve(ext_roc_list[["GSE135065_Model1"]], "GSE135065 (Indep)", "Ridge Logistic")
roc_df_list[[8]] <- add_roc_curve(ext_roc_list[["GSE135065_Model2"]], "GSE135065 (Indep)", "Unweighted Score")

all_roc_df <- do.call(rbind, roc_df_list)
all_roc_df$Cohort_Label <- paste0(all_roc_df$Cohort, " (AUC = ", sprintf("%.2f", all_roc_df$AUC), ")")

p_roc <- ggplot(all_roc_df, aes(x = FPR, y = TPR, color = Cohort)) +
  geom_line(lwd = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Model) +
  scale_color_manual(values = c("GSE32539 (CV)" = "#1f77b4", "GSE110147" = "#2ca02c", 
                                "GSE24206" = "#ff7f0e", "GSE135065 (Indep)" = "#d62728")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Multi-Cohort ROC Curves: 15-Gene Diagnostic Classifier",
    subtitle = "Left: Ridge Logistic Regression | Right: Unweighted Directional Score",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme(legend.position = "bottom")

save_plot(p_roc, "results/figures/prespecified/classifier_roc_curves", width = 10, height = 6)

# 7.2 Null Distributions Comparison Plot
null_plot_df <- rbind(
  data.frame(Distribution = "1,000 Network Random Sets", AUC = auc_null_network),
  data.frame(Distribution = "1,000 Matrisome Random Sets", AUC = auc_null_matrisome),
  data.frame(Distribution = "1,000 Shared-45 DEG Sets", AUC = auc_null_deg45),
  data.frame(Distribution = "200 Label Permutations", AUC = auc_label_perm)
)

p_null <- ggplot(null_plot_df, aes(x = AUC, fill = Distribution)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", position = "identity") +
  geom_vline(xintercept = obs_auc_15, color = "red", linetype = "dashed", lwd = 1.1) +
  annotate("text", x = obs_auc_15 - 0.02, y = 120, label = paste0("15-Gene Observed\nAUC = ", round(obs_auc_15, 3)),
           color = "red", fontface = "bold", hjust = 1, size = 3.8) +
  facet_wrap(~ Distribution, scales = "free_y") +
  scale_fill_manual(values = c("#6baed6", "#74c476", "#fd8d3c", "#9e9ac8")) +
  theme_bw(base_size = 11) +
  labs(
    title = "Empirical Null Distributions vs Observed 15-Gene Classifier Performance",
    subtitle = "Dashed red line = Observed 15-gene directional score AUC on GSE32539",
    x = "Cross-Validated AUC",
    y = "Count"
  ) +
  theme(legend.position = "none")

save_plot(p_null, "results/figures/prespecified/classifier_null_distributions", width = 10, height = 7)

# 7.3 Leave-One-Gene-Out Delta AUC Bar Plot
loo_plot_df <- loo_df
loo_plot_df$omitted_gene <- factor(loo_plot_df$omitted_gene, levels = loo_plot_df$omitted_gene[order(loo_plot_df$delta_AUC_GSE135065)])

p_loo <- ggplot(loo_plot_df, aes(x = delta_AUC_GSE135065, y = omitted_gene, fill = delta_AUC_GSE135065 > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"),
                    labels = c("TRUE" = "Dropping improves AUC", "FALSE" = "Dropping hurts AUC")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Leave-One-Gene-Out Impact on Independent Validation Cohort (GSE135065)",
    subtitle = "Negative Delta AUC = Gene contributes positively to independent classification",
    x = "Delta AUC (AUC without gene - Baseline AUC)",
    y = "Omitted Gene",
    fill = "Impact"
  ) +
  theme(legend.position = "bottom")

save_plot(p_loo, "results/figures/prespecified/leave_one_out_delta_auc", width = 9, height = 6)

# 7.4 Feature Importance Ranking Bar Plot
feat_imp_df$gene <- factor(feat_imp_df$gene, levels = rev(feat_imp_df$gene))

p_imp_rf <- ggplot(feat_imp_df, aes(x = RF_MeanDecreaseAccuracy, y = gene)) +
  geom_col(fill = "#1f77b4", width = 0.7) +
  theme_bw(base_size = 11) +
  labs(title = "Random Forest (MDA)", x = "Mean Decrease Accuracy", y = "Gene")

p_imp_svm <- ggplot(feat_imp_df, aes(x = SVM_Linear_Weight, y = gene)) +
  geom_col(fill = "#ff7f0e", width = 0.7) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Linear SVM Feature Weight", x = "Absolute Weight (|w|)", y = "")

p_imp_ridge <- ggplot(feat_imp_df, aes(x = abs(Ridge_Beta), y = gene)) +
  geom_col(fill = "#2ca02c", width = 0.7) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Ridge Regression (|Beta|)", x = "|Standardized Beta|", y = "")

p_imp_combo <- p_imp_rf + p_imp_svm + p_imp_ridge + plot_layout(ncol = 3) +
  plot_annotation(title = "Multi-Model Feature Importance Ranking: 15 Pre-Specified Genes (GSE32539)")

save_plot(p_imp_combo, "results/figures/prespecified/feature_importance_comparison", width = 11, height = 6)

log_msg("=== PART 3 & PART 4 Complete Successfully ===")

