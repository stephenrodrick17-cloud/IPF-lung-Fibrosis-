#!/usr/bin/env Rscript
## ==============================================================================
## 07_machine_learning.R
## Machine Learning Diagnostic Biomarker Discovery for IPF Matrisome Study
## Models: Random Forest (RF) & Support Vector Machine Recursive Feature Elimination (SVM-RFE)
## Training Cohort: GSE32539 (119 IPF, 49 Control; N=168)
## External Validation 1: GSE110147 (22 IPF, 11 Control; N=33)
## External Validation 2: GSE24206 (9 IPF, 6 Control; N=15)
## Seed: 20261102 | R version: 4.6.2 (Windows)
## ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1: Setup & Environment
# ------------------------------------------------------------------------------
set.seed(20261102)
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(randomForest)
  library(e1071)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

# Directories
dirs <- c("logs", "results/tables/ml", "results/figures/ml")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- "logs/ml.log"
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("=== STEP 1: Initializing Machine Learning Pipeline ===")
log_msg("R version: ", R.version.string)
log_msg("RandomForest version: ", as.character(packageVersion("randomForest")))
log_msg("e1071 version: ", as.character(packageVersion("e1071")))
log_msg("pROC version: ", as.character(packageVersion("pROC")))
log_msg("Random seed set to 20261102")

# Save session info
sink("logs/sessionInfo_ml.txt")
cat("=== Machine Learning Session Info ===\nSeed: 20261102\n\n")
print(sessionInfo())
sink()

# Helper for dual format figures (PNG 300 dpi + PDF)
save_plot <- function(p, file_base, width = 8, height = 6) {
  ggsave(paste0(file_base, ".png"), plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

# ------------------------------------------------------------------------------
# STEP 2: Data Ingestion & Preprocessing
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 2: Loading & Aligning Data Matrices ===")

# Expression matrices (genes in rows)
e32_df <- read.csv("data/processed/GSE32539_expr_gene_level.csv", check.names = FALSE)
rownames(e32_df) <- e32_df$gene
e32_df$gene <- NULL

e110_df <- read.csv("data/processed/GSE110147_expr_gene_level.csv", check.names = FALSE)
rownames(e110_df) <- e110_df$gene
e110_df$gene <- NULL

e24_df <- read.csv("data/processed/GSE24206_expr_gene_level.csv", check.names = FALSE)
rownames(e24_df) <- e24_df$gene
e24_df$gene <- NULL

# Metadata labels
l32 <- read.csv("data/processed/ml_input/GSE32539_labels.csv")
l110 <- read.csv("data/processed/ml_input/GSE110147_labels.csv")
l24 <- read.csv("data/processed/ml_input/GSE24206_labels.csv")

# Remove outlier GSM806411 in GSE32539
keep32 <- setdiff(l32$sample_id, "GSM806411")
l32 <- l32[match(keep32, l32$sample_id), ]

log_msg("GSE32539 (Discovery): ", nrow(l32), " samples (", sum(l32$label == "IPF"), " IPF, ", sum(l32$label == "Control"), " Control)")
log_msg("GSE110147 (Validation 1): ", nrow(l110), " samples (", sum(l110$label == "IPF"), " IPF, ", sum(l110$label == "Control"), " Control)")
log_msg("GSE24206 (Validation 2): ", nrow(l24), " samples (", sum(l24$label == "IPF"), " IPF, ", sum(l24$label == "Control"), " Control)")

# Candidate feature space: Union of 45 shared matrisome DEGs and 60 WGCNA matrisome hubs
shared45 <- read.csv("results/tables/common_matrisome_DEGs_preview.csv")$gene
hubs_mat <- read.csv("results/tables/wgcna/hubs_matrisome.csv")$gene
cand_pool <- sort(unique(c(shared45, hubs_mat)))

# Exclude genes missing on any platform (e.g. SFTPA1 on GPL6244)
common_genes <- intersect(cand_pool, intersect(rownames(e32_df), intersect(rownames(e110_df), rownames(e24_df))))
log_msg("Total candidate matrisome pool: ", length(cand_pool), " genes")
log_msg("Candidate genes present across all 3 platforms: ", length(common_genes), " genes")

# Prepare training data: Samples in rows, Genes in columns
train_X <- t(as.matrix(e32_df[common_genes, l32$sample_id]))
train_y <- factor(l32$label, levels = c("Control", "IPF"))

# ------------------------------------------------------------------------------
# STEP 3: Random Forest Feature Selection
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 3: Random Forest Feature Selection ===")

# 3.1 Fit full Random Forest (1000 trees)
rf_model <- randomForest(
  x = train_X,
  y = train_y,
  ntree = 1000,
  importance = TRUE,
  proximity = TRUE,
  seed = 20261102
)

rf_imp <- importance(rf_model)
rf_imp_df <- data.frame(
  gene = rownames(rf_imp),
  MeanDecreaseAccuracy = round(rf_imp[, "MeanDecreaseAccuracy"], 4),
  MeanDecreaseGini = round(rf_imp[, "MeanDecreaseGini"], 4)
)
rf_imp_df <- rf_imp_df[order(-rf_imp_df$MeanDecreaseGini), ]
rf_imp_df$rank_RF <- 1:nrow(rf_imp_df)

write.csv(rf_imp_df, "results/tables/ml/rf_variable_importance.csv", row.names = FALSE)
log_msg("Saved results/tables/ml/rf_variable_importance.csv")

# 3.2 10-fold Cross-Validation with rfcv
log_msg("Running 10-fold cross-validation with rfcv to determine optimal feature count...")
set.seed(20261102)
rf_cv_res <- rfcv(train_X, train_y, cv.fold = 10, step = 0.8)

rf_cv_df <- data.frame(
  n_genes = as.numeric(names(rf_cv_res$error.cv)),
  cv_error = as.numeric(rf_cv_res$error.cv)
)
rf_cv_df <- rf_cv_df[order(rf_cv_df$n_genes), ]

best_rf_n <- rf_cv_df$n_genes[which.min(rf_cv_df$cv_error)[1]]
min_rf_err <- min(rf_cv_df$cv_error)
log_msg("Random Forest 10-fold CV: Minimum error = ", round(min_rf_err, 4), " achieved at n = ", best_rf_n, " genes.")

# Plot 1: RF CV Error Curve
p_rf_cv <- ggplot(rf_cv_df, aes(x = n_genes, y = cv_error)) +
  geom_line(color = "#1f78b4", lwd = 1.2) +
  geom_point(color = "#1f78b4", size = 3) +
  geom_vline(xintercept = best_rf_n, linetype = "dashed", color = "red", lwd = 1) +
  annotate("text", x = best_rf_n + 5, y = max(rf_cv_df$cv_error) * 0.9, 
           label = paste0("Optimal: n = ", best_rf_n, "\nError = ", round(min_rf_err, 4)), 
           color = "red", fontface = "bold", hjust = 0) +
  scale_x_continuous(trans = "log10", breaks = c(1, 2, 5, 10, 20, 40, 80)) +
  theme_bw(base_size = 13) +
  labs(
    title = "Random Forest: 10-Fold Cross-Validation Error Curve",
    x = "Number of Features (log10 scale)",
    y = "10-Fold CV Error Rate"
  )
save_plot(p_rf_cv, "results/figures/ml/rf_cv_error_curve", width = 8, height = 5.5)

# Plot 2: Top 25 RF Features by MeanDecreaseGini
top25_rf <- head(rf_imp_df, 25)
top25_rf$gene <- factor(top25_rf$gene, levels = rev(top25_rf$gene))
p_rf_imp <- ggplot(top25_rf, aes(x = MeanDecreaseGini, y = gene)) +
  geom_col(fill = "#33a02c", width = 0.7) +
  theme_bw(base_size = 12) +
  labs(
    title = "Top 25 Diagnostic Features by Random Forest Importance",
    x = "Mean Decrease Gini",
    y = "Gene Symbol"
  )
save_plot(p_rf_imp, "results/figures/ml/rf_variable_importance_top25", width = 7, height = 8)

# Top RF candidate set (using optimal CV size or top 20)
top_rf_genes <- head(rf_imp_df$gene, max(best_rf_n, 20))
log_msg("Selected top ", length(top_rf_genes), " Random Forest candidate features.")

# ------------------------------------------------------------------------------
# STEP 4: Support Vector Machine Recursive Feature Elimination (SVM-RFE)
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 4: SVM-RFE Feature Selection ===")

# SVM-RFE Ranking Function
svm_rfe_rank <- function(X, y) {
  n_genes <- ncol(X)
  genes_left <- colnames(X)
  ranked <- character(n_genes)
  y_num <- ifelse(y == "IPF", 1, -1)
  
  for (step in n_genes:1) {
    if (step == 1) { ranked[1] <- genes_left[1]; break }
    fit <- svm(x = X[, genes_left, drop = FALSE], y = y_num, type = "C-classification", 
               kernel = "linear", cost = 1, scale = FALSE)
    w <- t(fit$coefs) %*% fit$SV
    w2 <- as.vector(w^2)
    min_idx <- which.min(w2)
    ranked[step] <- genes_left[min_idx]
    genes_left <- genes_left[-min_idx]
  }
  return(ranked)
}

# 4.1 10-Fold Cross-Validation for SVM-RFE
log_msg("Running 10-fold cross-validation for SVM-RFE...")
set.seed(20261102)
k_folds <- 10
n_samples <- nrow(train_X)
n_total_genes <- ncol(train_X)
folds <- sample(rep(1:k_folds, length.out = n_samples))
y_num <- ifelse(train_y == "IPF", 1, -1)

subset_sizes <- c(1:15, seq(16, n_total_genes, by = 2))
cv_acc_matrix <- matrix(0, nrow = k_folds, ncol = length(subset_sizes))
colnames(cv_acc_matrix) <- subset_sizes

for (f in 1:k_folds) {
  train_idx <- which(folds != f)
  test_idx <- which(folds == f)
  
  X_tr <- train_X[train_idx, , drop = FALSE]
  y_tr <- y_num[train_idx]
  X_te <- train_X[test_idx, , drop = FALSE]
  y_te <- y_num[test_idx]
  
  # Rank features on training fold
  ranked_fold <- svm_rfe_rank(X_tr, factor(ifelse(y_tr == 1, "IPF", "Control"), levels = c("Control", "IPF")))
  
  # Evaluate test accuracy per subset size
  for (s_idx in seq_along(subset_sizes)) {
    s <- subset_sizes[s_idx]
    sub_genes <- ranked_fold[1:s]
    fit_sub <- svm(x = X_tr[, sub_genes, drop = FALSE], y = y_tr, type = "C-classification", 
                   kernel = "linear", cost = 1, scale = FALSE)
    preds <- predict(fit_sub, newdata = X_te[, sub_genes, drop = FALSE])
    cv_acc_matrix[f, s_idx] <- mean(preds == y_te)
  }
}

mean_cv_acc <- colMeans(cv_acc_matrix)
svm_cv_df <- data.frame(
  n_genes = as.numeric(colnames(cv_acc_matrix)),
  cv_accuracy = as.numeric(mean_cv_acc),
  cv_error = 1 - as.numeric(mean_cv_acc)
)

best_svm_n <- svm_cv_df$n_genes[which.max(svm_cv_df$cv_accuracy)[1]]
max_svm_acc <- max(svm_cv_df$cv_accuracy)
log_msg("SVM-RFE 10-fold CV: Maximum accuracy = ", round(max_svm_acc * 100, 2), "% achieved at n = ", best_svm_n, " genes.")

# Plot 3: SVM-RFE CV Accuracy Curve
p_svm_cv <- ggplot(svm_cv_df, aes(x = n_genes, y = cv_accuracy * 100)) +
  geom_line(color = "#e31a1c", lwd = 1.2) +
  geom_point(color = "#e31a1c", size = 3) +
  geom_vline(xintercept = best_svm_n, linetype = "dashed", color = "darkgreen", lwd = 1) +
  annotate("text", x = best_svm_n + 2, y = min(svm_cv_df$cv_accuracy * 100) + 1, 
           label = paste0("Optimal: n = ", best_svm_n, "\nAccuracy = ", round(max_svm_acc * 100, 1), "%"), 
           color = "darkgreen", fontface = "bold", hjust = 0) +
  theme_bw(base_size = 13) +
  labs(
    title = "SVM-RFE: 10-Fold Cross-Validation Accuracy Curve",
    x = "Number of Features",
    y = "10-Fold CV Accuracy (%)"
  )
save_plot(p_svm_cv, "results/figures/ml/svm_rfe_cv_accuracy", width = 8, height = 5.5)

# 4.2 Full SVM-RFE Ranking on GSE32539
log_msg("Computing final SVM-RFE feature ranking across all training samples...")
ranked_svm_genes <- svm_rfe_rank(train_X, train_y)
svm_ranking_df <- data.frame(
  rank_SVM = 1:length(ranked_svm_genes),
  gene = ranked_svm_genes
)
write.csv(svm_ranking_df, "results/tables/ml/svm_rfe_ranking.csv", row.names = FALSE)
log_msg("Saved results/tables/ml/svm_rfe_ranking.csv")

top_svm_genes <- head(ranked_svm_genes, max(best_svm_n, 20))
log_msg("Selected top ", length(top_svm_genes), " SVM-RFE candidate features.")

# ------------------------------------------------------------------------------
# STEP 5: Consensus Diagnostic Biomarker Identification
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 5: Consensus Diagnostic Biomarkers (RF & SVM-RFE Intersection) ===")

# Combine rankings
merged_ranking <- merge(rf_imp_df, svm_ranking_df, by = "gene")
merged_ranking$rank_product <- sqrt(merged_ranking$rank_RF * merged_ranking$rank_SVM)
merged_ranking <- merged_ranking[order(merged_ranking$rank_product), ]

# Primary consensus: Top 25 from RF and Top 25 from SVM-RFE
top_rf_set <- head(rf_imp_df$gene, 25)
top_svm_set <- head(ranked_svm_genes, 25)
consensus_biomarkers <- intersect(top_rf_set, top_svm_set)

# Ensure at least 8-12 robust biomarkers by taking top by rank product
if (length(consensus_biomarkers) < 8) {
  consensus_biomarkers <- head(merged_ranking$gene, 10)
}

log_msg("Identified ", length(consensus_biomarkers), " consensus diagnostic biomarker genes:")
log_msg("  ", paste(consensus_biomarkers, collapse = ", "))

# Save consensus biomarker table with full annotations
biomarker_table <- merged_ranking[merged_ranking$gene %in% consensus_biomarkers, ]
biomarker_table <- biomarker_table[order(biomarker_table$rank_product), ]
write.csv(biomarker_table, "results/tables/ml/diagnostic_biomarkers.csv", row.names = FALSE)
writeLines(consensus_biomarkers, "results/tables/ml/diagnostic_biomarkers.txt")
log_msg("Saved results/tables/ml/diagnostic_biomarkers.csv and .txt")

# Plot 4: Venn / Intersection Representation
venn_data <- data.frame(
  Category = c("Random Forest (Top 25)", "Consensus Overlap", "SVM-RFE (Top 25)"),
  Count = c(length(setdiff(top_rf_set, top_svm_set)), length(consensus_biomarkers), length(setdiff(top_svm_set, top_rf_set)))
)
p_overlap <- ggplot(venn_data, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Random Forest (Top 25)" = "#33a02c", "Consensus Overlap" = "#6a3d9a", "SVM-RFE (Top 25)" = "#e31a1c")) +
  theme_bw(base_size = 13) +
  ylim(0, max(venn_data$Count) + 5) +
  labs(
    title = "Machine Learning Feature Selection Concordance",
    subtitle = paste0(length(consensus_biomarkers), " Diagnostic Genes Identified at Intersection"),
    x = "",
    y = "Number of Genes"
  )
save_plot(p_overlap, "results/figures/ml/rf_svm_venn", width = 7, height = 5)

# ------------------------------------------------------------------------------
# STEP 6: Multi-Cohort Diagnostic Validation (ROC & AUC Analysis)
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 6: Multi-Cohort Validation (ROC & AUC Analysis) ===")

# Build standardized diagnostic risk scores across all cohorts:
# Direction flag: UP in IPF adds to score, DOWN subtracts
limma32 <- read.csv("results/tables/GSE32539_limma_full.csv")
rownames(limma32) <- limma32$gene
gene_direction <- ifelse(limma32[consensus_biomarkers, "logFC"] > 0, 1, -1)
names(gene_direction) <- consensus_biomarkers

calc_signature_score <- function(e_mat, sample_ids, genes, directions) {
  sub_mat <- t(as.matrix(e_mat[genes, sample_ids]))
  z_mat <- scale(sub_mat)
  # Weighted sum by direction
  score <- apply(z_mat, 1, function(row) sum(row * directions[genes]))
  return(score)
}

score32 <- calc_signature_score(e32_df, l32$sample_id, consensus_biomarkers, gene_direction)
score110 <- calc_signature_score(e110_df, l110$sample_id, consensus_biomarkers, gene_direction)
score24 <- calc_signature_score(e24_df, l24$sample_id, consensus_biomarkers, gene_direction)

# ROC calculations
roc32 <- roc(train_y, score32, levels = c("Control", "IPF"), direction = "<", quiet = TRUE)
roc110 <- roc(factor(l110$label, levels = c("Control", "IPF")), score110, levels = c("Control", "IPF"), direction = "<", quiet = TRUE)
roc24 <- roc(factor(l24$label, levels = c("Control", "IPF")), score24, levels = c("Control", "IPF"), direction = "<", quiet = TRUE)

ci32 <- ci.auc(roc32)
ci110 <- ci.auc(roc110)
ci24 <- ci.auc(roc24)

# Diagnostic metrics function
get_diag_metrics <- function(roc_obj, cohort_name, ci_obj) {
  coords_best <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", "accuracy", "ppv", "npv"))
  if (is.matrix(coords_best) || is.data.frame(coords_best)) coords_best <- coords_best[1, ]
  data.frame(
    Cohort = cohort_name,
    AUC = round(as.numeric(auc(roc_obj)), 4),
    CI_95_low = round(as.numeric(ci_obj[1]), 4),
    CI_95_high = round(as.numeric(ci_obj[3]), 4),
    Sensitivity = round(as.numeric(coords_best["sensitivity"]), 4),
    Specificity = round(as.numeric(coords_best["specificity"]), 4),
    Accuracy = round(as.numeric(coords_best["accuracy"]), 4),
    PPV = round(as.numeric(coords_best["ppv"]), 4),
    NPV = round(as.numeric(coords_best["npv"]), 4)
  )
}

perf_df <- rbind(
  get_diag_metrics(roc32, "GSE32539 (Discovery)", ci32),
  get_diag_metrics(roc110, "GSE110147 (Validation 1)", ci110),
  get_diag_metrics(roc24, "GSE24206 (Validation 2)", ci24)
)

write.csv(perf_df, "results/tables/ml/diagnostic_performance_summary.csv", row.names = FALSE)
log_msg("Saved results/tables/ml/diagnostic_performance_summary.csv")

log_msg("Diagnostic Signature Performance:")
for (i in 1:nrow(perf_df)) {
  log_msg(sprintf("  [%s] AUC = %.3f (95%% CI: %.3f - %.3f), Sens = %.1f%%, Spec = %.1f%%",
                  perf_df$Cohort[i], perf_df$AUC[i], perf_df$CI_95_low[i], perf_df$CI_95_high[i],
                  perf_df$Sensitivity[i] * 100, perf_df$Specificity[i] * 100))
}

# Plot 5: Multi-Cohort Overlayed ROC Curves
roc_df_all <- rbind(
  data.frame(Specificities = 1 - roc32$specificities, Sensitivities = roc32$sensitivities, Cohort = paste0("GSE32539 (Discovery: AUC = ", sprintf("%.3f", auc(roc32)), ")")),
  data.frame(Specificities = 1 - roc110$specificities, Sensitivities = roc110$sensitivities, Cohort = paste0("GSE110147 (Validation 1: AUC = ", sprintf("%.3f", auc(roc110)), ")")),
  data.frame(Specificities = 1 - roc24$specificities, Sensitivities = roc24$sensitivities, Cohort = paste0("GSE24206 (Validation 2: AUC = ", sprintf("%.3f", auc(roc24)), ")"))
)

p_roc <- ggplot(roc_df_all, aes(x = Specificities, y = Sensitivities, color = Cohort)) +
  geom_line(lwd = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c("#2b83ba", "#d7191c", "#fdae61")) +
  theme_bw(base_size = 13) +
  labs(
    title = "Multi-Cohort ROC Curves: Matrisome Diagnostic Signature",
    subtitle = paste0(length(consensus_biomarkers), "-Gene Composite Risk Score"),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Cohort"
  ) +
  theme(legend.position = c(0.68, 0.22), legend.background = element_rect(fill = alpha("white", 0.8), color = "gray80"))

save_plot(p_roc, "results/figures/ml/roc_multicohort_signature", width = 8, height = 6.5)

# ------------------------------------------------------------------------------
# STEP 7: Single-Gene Diagnostic Performance
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 7: Single-Gene Diagnostic Performance Analysis ===")

single_gene_perf <- list()

for (g in consensus_biomarkers) {
  # GSE32539
  r32 <- roc(train_y, as.numeric(e32_df[g, l32$sample_id]), levels = c("Control", "IPF"), quiet = TRUE)
  # GSE110147
  r110 <- roc(factor(l110$label, levels = c("Control", "IPF")), as.numeric(e110_df[g, l110$sample_id]), levels = c("Control", "IPF"), quiet = TRUE)
  # GSE24206
  r24 <- roc(factor(l24$label, levels = c("Control", "IPF")), as.numeric(e24_df[g, l24$sample_id]), levels = c("Control", "IPF"), quiet = TRUE)
  
  single_gene_perf[[g]] <- data.frame(
    gene = g,
    AUC_GSE32539 = round(as.numeric(auc(r32)), 4),
    AUC_GSE110147 = round(as.numeric(auc(r110)), 4),
    AUC_GSE24206 = round(as.numeric(auc(r24)), 4),
    Mean_AUC = round(mean(c(as.numeric(auc(r32)), as.numeric(auc(r110)), as.numeric(auc(r24)))), 4),
    Direction = ifelse(gene_direction[g] == 1, "Up", "Down")
  )
}

single_gene_df <- do.call(rbind, single_gene_perf)
single_gene_df <- single_gene_df[order(-single_gene_df$Mean_AUC), ]
write.csv(single_gene_df, "results/tables/ml/single_gene_diagnostic_performance.csv", row.names = FALSE)
log_msg("Saved results/tables/ml/single_gene_diagnostic_performance.csv")

log_msg("Single-Gene Performance Summary:")
for (i in 1:nrow(single_gene_df)) {
  log_msg(sprintf("  %-10s: Mean AUC = %.3f (GSE32539: %.3f, GSE110147: %.3f, GSE24206: %.3f) [%s]",
                  single_gene_df$gene[i], single_gene_df$Mean_AUC[i],
                  single_gene_df$AUC_GSE32539[i], single_gene_df$AUC_GSE110147[i], single_gene_df$AUC_GSE24206[i],
                  single_gene_df$Direction[i]))
}

# Plot 6: Single Gene AUC Dot Plot
single_gene_df$gene <- factor(as.character(single_gene_df$gene), levels = rev(as.character(single_gene_df$gene)))
p_single_auc <- ggplot(single_gene_df, aes(x = Mean_AUC, y = gene, color = Direction)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = pmin(AUC_GSE32539, AUC_GSE110147, AUC_GSE24206),
                    xmax = pmax(AUC_GSE32539, AUC_GSE110147, AUC_GSE24206)), 
                width = 0.2, lwd = 0.8, orientation = "y") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Up" = "#d7191c", "Down" = "#2b83ba")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Single-Gene Diagnostic Discrimination (Cross-Cohort AUC)",
    subtitle = "Point = Mean AUC across 3 cohorts; Bar = Range across cohorts",
    x = "Area Under the Curve (AUC)",
    y = "Diagnostic Gene Symbol"
  )
save_plot(p_single_auc, "results/figures/ml/single_gene_auc_summary", width = 8, height = 6)

# ------------------------------------------------------------------------------
# STEP 8: Gene Expression Boxplots across Cohorts
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 8: Generating Expression Boxplots across Cohorts ===")

# Build long-format expression data for top 6 consensus genes
top6_genes <- as.character(head(single_gene_df$gene, 6))

get_long_expr <- function(e_mat, l_df, cohort_name, genes) {
  res_list <- list()
  for (g in genes) {
    res_list[[g]] <- data.frame(
      sample_id = l_df$sample_id,
      cohort = cohort_name,
      gene = g,
      expression = as.numeric(unlist(e_mat[g, l_df$sample_id])),
      group = l_df$label
    )
  }
  do.call(rbind, res_list)
}

expr_long_all <- rbind(
  get_long_expr(e32_df, l32, "GSE32539", top6_genes),
  get_long_expr(e110_df, l110, "GSE110147", top6_genes),
  get_long_expr(e24_df, l24, "GSE24206", top6_genes)
)

expr_long_all$group <- factor(expr_long_all$group, levels = c("Control", "IPF"))
expr_long_all$cohort <- factor(expr_long_all$cohort, levels = c("GSE32539", "GSE110147", "GSE24206"))

p_boxplots <- ggplot(expr_long_all, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  facet_grid(gene ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#2b83ba", "IPF" = "#d7191c")) +
  theme_bw(base_size = 11) +
  labs(
    title = "Cross-Cohort Expression of Top Diagnostic Matrisome Genes",
    subtitle = "Consistent Dysregulation across Discovery (GSE32539) & Validation (GSE110147, GSE24206)",
    x = "Clinical Group",
    y = "Log2 Expression",
    fill = "Group"
  ) +
  theme(strip.text = element_text(face = "bold"))

save_plot(p_boxplots, "results/figures/ml/expression_boxplots_top_diagnostic_genes", width = 9, height = 11)

log_msg("\n=== Machine Learning Diagnostic Biomarker Pipeline Complete Successfully ===")
