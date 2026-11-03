#!/usr/bin/env Rscript
## ==============================================================================
## 09_prespecified_auc.R
## PART 2: Single-Gene Diagnostic Performance of the 15 Pre-Specified Genes
## Seed: 20261102 | R version: 4.6.2 (Windows)
## ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0: Pre-Specified Rules & Setup
# ------------------------------------------------------------------------------
# Significance criterion: Wilcoxon rank-sum p-value, Benjamini-Hochberg adj.p < 0.05
# Direction-free AUC = max(AUC, 1 - AUC), with recorded direction (+ or -)
# 95% Confidence Intervals via DeLong method (pROC)
# Independent Cohort: GSE135065 (strictly independent; NOT used in DEG selection)
# Discovery Cohorts: GSE32539, GSE110147, GSE24206 (used in prior DEG screening)
# ------------------------------------------------------------------------------

set.seed(20261102)
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(pROC)
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

log_msg("\n=== Initializing PART 2: Single-Gene AUC Testing (15 Genes) ===")

save_plot <- function(p, file_base, width = 10, height = 7) {
  ggsave(paste0(file_base, ".png"), plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

prespecified_15 <- c(
  "ASPN", "COL14A1", "COL15A1", "COL1A1", "COL3A1", "COL5A2", "COMP",
  "IL6", "MMP7", "MMP8", "MMP10", "OGN", "SFRP4", "TGFB3", "TIMP3"
)

# ------------------------------------------------------------------------------
# STEP 1: Load Datasets and Metadata
# ------------------------------------------------------------------------------
datasets <- c("GSE32539", "GSE110147", "GSE24206", "GSE135065")
data_list <- list()

for (gse in datasets) {
  expr_file <- paste0("data/processed/", gse, "_expr_gene_level.csv")
  meta_file <- paste0("data/processed/", gse, "_metadata.csv")
  limma_file <- paste0("results/tables/", gse, "_limma_full.csv")
  
  e <- read.csv(expr_file, check.names = FALSE)
  rownames(e) <- e$gene
  e$gene <- NULL
  
  meta <- read.csv(meta_file)
  
  # Normalize label/group column
  if ("group" %in% colnames(meta)) {
    group_col <- meta$group
  } else if ("disease_state" %in% colnames(meta)) {
    group_col <- meta$disease_state
  } else if ("label" %in% colnames(meta)) {
    group_col <- meta$label
  } else {
    stop("Unknown label column in metadata for ", gse)
  }
  
  names(group_col) <- meta$sample_id
  
  # For GSE32539, drop the outlier GSM806411 to match WGCNA 168 sample cohort
  if (gse == "GSE32539") {
    valid_samples <- intersect(colnames(e), setdiff(meta$sample_id, "GSM806411"))
  } else {
    valid_samples <- intersect(colnames(e), meta$sample_id)
  }
  
  group_sub <- group_col[valid_samples]
  # Filter to IPF and Control only
  keep <- group_sub %in% c("IPF", "Control")
  valid_samples <- valid_samples[keep]
  group_sub <- factor(group_sub[keep], levels = c("Control", "IPF"))
  
  limma_df <- read.csv(limma_file)
  rownames(limma_df) <- limma_df$gene
  
  data_list[[gse]] <- list(
    expr = e[, valid_samples],
    group = group_sub,
    limma = limma_df,
    is_independent = (gse == "GSE135065")
  )
  
  log_msg("Loaded ", gse, ": ", sum(group_sub == "IPF"), " IPF, ", sum(group_sub == "Control"), " Control. Independent: ", gse == "GSE135065")
}

# ------------------------------------------------------------------------------
# STEP 2: Compute AUC, CIs, Wilcoxon p-values per Gene per Dataset
# ------------------------------------------------------------------------------
log_msg("\n--- Computing AUC and statistical tests across datasets ---")

auc_results <- list()
long_plot_data <- list()

for (gse in datasets) {
  d <- data_list[[gse]]
  e <- d$expr
  grp <- d$group
  lim <- d$limma
  is_indep <- d$is_independent
  
  wilcox_p_vals <- numeric(length(prespecified_15))
  names(wilcox_p_vals) <- prespecified_15
  
  for (g in prespecified_15) {
    if (g %in% rownames(e)) {
      val <- as.numeric(e[g, names(grp)])
      wt <- wilcox.test(val ~ grp, exact = FALSE)
      wilcox_p_vals[g] <- wt$p.value
    } else {
      wilcox_p_vals[g] <- NA
    }
  }
  
  # Benjamini-Hochberg adjustment across 15 genes
  adj_wilcox <- p.adjust(wilcox_p_vals, method = "BH")
  
  for (g in prespecified_15) {
    if (g %in% rownames(e)) {
      val <- as.numeric(e[g, names(grp)])
      
      # Limma stats
      if (g %in% rownames(lim)) {
        l_fc <- lim[g, "logFC"]
        l_adj_p <- lim[g, "adj.P.Val"]
      } else {
        l_fc <- NA
        l_adj_p <- NA
      }
      
      direction <- ifelse(l_fc > 0, "UP", "DOWN")
      
      # ROC curve
      roc_obj <- suppressMessages(pROC::roc(response = grp, predictor = val, levels = c("Control", "IPF"), direction = "<"))
      raw_auc <- as.numeric(pROC::auc(roc_obj))
      
      # Direction-free AUC = max(AUC, 1-AUC)
      if (raw_auc >= 0.5) {
        auc_val <- raw_auc
        auc_dir <- "+"
        ci_obj <- suppressMessages(pROC::ci.auc(roc_obj, method = "delong"))
        ci_lower <- as.numeric(ci_obj[1])
        ci_upper <- as.numeric(ci_obj[3])
      } else {
        auc_val <- 1 - raw_auc
        auc_dir <- "-"
        # Flip predictor to calculate correct CI for reversed direction
        roc_flip <- suppressMessages(pROC::roc(response = grp, predictor = -val, levels = c("Control", "IPF"), direction = "<"))
        ci_obj <- suppressMessages(pROC::ci.auc(roc_flip, method = "delong"))
        ci_lower <- as.numeric(ci_obj[1])
        ci_upper <- as.numeric(ci_obj[3])
      }
      
      w_p <- wilcox_p_vals[g]
      w_adj_p <- adj_wilcox[g]
      
      auc_results[[paste(gse, g, sep = "_")]] <- data.frame(
        dataset = gse,
        cohort_type = ifelse(is_indep, "Independent_Validation", "Discovery_Cohort"),
        gene = g,
        direction_in_IPF = direction,
        logFC = round(l_fc, 3),
        limma_adj_p = signif(l_adj_p, 4),
        raw_AUC = round(raw_auc, 4),
        AUC = round(auc_val, 4),
        AUC_direction = auc_dir,
        AUC_CI_lower = round(ci_lower, 4),
        AUC_CI_upper = round(ci_upper, 4),
        Wilcoxon_p = signif(w_p, 4),
        Wilcoxon_adj_p = signif(w_adj_p, 4),
        is_significant = (w_adj_p < 0.05)
      )
      
      # Store expression values for boxplots
      long_plot_data[[paste(gse, g, sep = "_")]] <- data.frame(
        dataset = gse,
        gene = g,
        sample_id = names(grp),
        group = grp,
        expression = val
      )
      
    } else {
      auc_results[[paste(gse, g, sep = "_")]] <- data.frame(
        dataset = gse,
        cohort_type = ifelse(is_indep, "Independent_Validation", "Discovery_Cohort"),
        gene = g,
        direction_in_IPF = "Not_measured",
        logFC = NA,
        limma_adj_p = NA,
        raw_AUC = NA,
        AUC = NA,
        AUC_direction = NA,
        AUC_CI_lower = NA,
        AUC_CI_upper = NA,
        Wilcoxon_p = NA,
        Wilcoxon_adj_p = NA,
        is_significant = FALSE
      )
    }
  }
}

auc_df <- do.call(rbind, auc_results)
plot_expr_df <- do.call(rbind, long_plot_data)

# ------------------------------------------------------------------------------
# STEP 3: Concordance, Consistency, and Four-Cohort Significance Analysis
# ------------------------------------------------------------------------------
log_msg("\n--- Assessing Direction Concordance & 4-Cohort Significance ---")

concordance_list <- list()

for (g in prespecified_15) {
  sub_g <- auc_df[auc_df$gene == g, ]
  
  dirs <- sub_g$direction_in_IPF
  names(dirs) <- sub_g$dataset
  
  # Check if significant in all 4 datasets
  sig_all4 <- all(sub_g$is_significant)
  
  # Check direction consistency across the 4 datasets
  valid_dirs <- dirs[dirs %in% c("UP", "DOWN")]
  is_consistent <- length(unique(valid_dirs)) == 1
  
  # Independent cohort stats
  indep_row <- sub_g[sub_g$dataset == "GSE135065", ]
  indep_auc <- indep_row$AUC
  indep_sig <- indep_row$is_significant
  indep_dir <- indep_row$direction_in_IPF
  
  concordance_list[[g]] <- data.frame(
    gene = g,
    dir_GSE32539 = dirs["GSE32539"],
    dir_GSE110147 = dirs["GSE110147"],
    dir_GSE24206 = dirs["GSE24206"],
    dir_GSE135065 = dirs["GSE135065"],
    direction_consistent_all4 = ifelse(is_consistent, "yes", "no"),
    significant_all4 = ifelse(sig_all4, "yes", "no"),
    independent_AUC_GSE135065 = indep_auc,
    independent_significant = ifelse(indep_sig, "yes", "no"),
    thesis_discordance_check = ifelse(g %in% c("IL6", "MMP8", "TIMP3"), "Thesis_flagged_discordant", "Thesis_concordant_candidate")
  )
}

concord_df <- do.call(rbind, concordance_list)

# Merge concordance summary back into auc_df
auc_df_merged <- merge(auc_df, concord_df, by = "gene")
# Order nicely
auc_df_merged <- auc_df_merged[order(auc_df_merged$gene, auc_df_merged$dataset), ]

# Save Table 2
write.csv(auc_df_merged, "results/tables/prespecified/single_gene_auc_15.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/single_gene_auc_15.csv")

# Print Summary Table of Concordance and AUCs
log_msg("\nDirection Concordance & Independent Validation Summary:")
print(concord_df)

# Explicit Check on IL6, MMP8, TIMP3
log_msg("\nExplicit Verification of Thesis-Flagged Discordant Genes:")
for (tg in c("IL6", "MMP8", "TIMP3")) {
  sub_tg <- auc_df[auc_df$gene == tg, c("dataset", "direction_in_IPF", "logFC", "AUC", "Wilcoxon_adj_p")]
  cat("\nGene: ", tg, "\n")
  print(sub_tg)
}

# ------------------------------------------------------------------------------
# STEP 4: Visualizations
# ------------------------------------------------------------------------------
log_msg("\n--- Generating PART 2 Figures ---")

# 4.1 Forest-style Plot of AUC with 95% CIs across all 4 datasets
auc_clean <- auc_df[!is.na(auc_df$AUC), ]
auc_clean$dataset <- factor(auc_clean$dataset, levels = c("GSE32539", "GSE110147", "GSE24206", "GSE135065"))
auc_clean$gene <- factor(auc_clean$gene, levels = rev(prespecified_15))

p_forest <- ggplot(auc_clean, aes(x = AUC, y = gene, color = dataset, shape = cohort_type)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0.8, linetype = "dotted", color = "blue", alpha = 0.7) +
  geom_point(position = position_dodge(width = 0.7), size = 2.5) +
  geom_errorbar(aes(xmin = AUC_CI_lower, xmax = AUC_CI_upper), 
                position = position_dodge(width = 0.7), width = 0.3, orientation = "y") +
  scale_color_manual(values = c("GSE32539" = "#1f77b4", "GSE110147" = "#2ca02c", 
                                "GSE24206" = "#ff7f0e", "GSE135065" = "#d62728")) +
  scale_shape_manual(values = c("Discovery_Cohort" = 16, "Independent_Validation" = 17)) +
  scale_x_continuous(limits = c(0.4, 1.05), breaks = seq(0.4, 1.0, 0.1)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Single-Gene Diagnostic Performance (AUC with 95% DeLong CI)",
    subtitle = "Red triangle = Strictly Independent Validation (GSE135065); Circles = Discovery Cohorts",
    x = "Direction-Free AUC [max(AUC, 1 - AUC)]",
    y = "Pre-Specified Gene",
    color = "Dataset",
    shape = "Cohort Role"
  ) +
  theme(legend.position = "right")

save_plot(p_forest, "results/figures/prespecified/single_gene_auc_forest", width = 11, height = 8)
log_msg("Saved AUC forest plot (PNG & PDF)")

# 4.2 Boxplots of Gene Expression Across 4 Datasets
# Make 2 multi-panel figures: Part A (genes 1:8), Part B (genes 9:15)
plot_expr_df$dataset <- factor(plot_expr_df$dataset, levels = c("GSE32539", "GSE110147", "GSE24206", "GSE135065"))

genes_partA <- prespecified_15[1:8]
genes_partB <- prespecified_15[9:15]

p_box_A <- ggplot(subset(plot_expr_df, gene %in% genes_partA), aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_grid(gene ~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4575b4", "IPF" = "#d73027")) +
  theme_bw(base_size = 11) +
  labs(
    title = "Expression Distributions: Pre-Specified Genes (Part 1: ASPN to IL6)",
    x = "Phenotype",
    y = "Log2 Gene Expression",
    fill = "Group"
  ) +
  theme(legend.position = "bottom")

save_plot(p_box_A, "results/figures/prespecified/expression_boxplots_part1", width = 11, height = 12)

p_box_B <- ggplot(subset(plot_expr_df, gene %in% genes_partB), aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_grid(gene ~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4575b4", "IPF" = "#d73027")) +
  theme_bw(base_size = 11) +
  labs(
    title = "Expression Distributions: Pre-Specified Genes (Part 2: MMP7 to TIMP3)",
    x = "Phenotype",
    y = "Log2 Gene Expression",
    fill = "Group"
  ) +
  theme(legend.position = "bottom")

save_plot(p_box_B, "results/figures/prespecified/expression_boxplots_part2", width = 11, height = 11)

log_msg("Saved expression boxplots (PNG & PDF)")
log_msg("=== PART 2 Complete Successfully ===")

