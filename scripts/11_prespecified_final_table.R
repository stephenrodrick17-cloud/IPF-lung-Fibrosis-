#!/usr/bin/env Rscript
## ==============================================================================
## 11_prespecified_final_table.R
## PART 5: Final Synthesis Table for the Pre-Specified 15-Gene Evaluation
## Seed: 20261102 | R version: 4.6.2 (Windows)
## ==============================================================================

# ------------------------------------------------------------------------------
# Fixed Final Evidence Classification Rules:
# - "Strong":
#     (1) Strict WGCNA Hub ("Hub (strict)")
#     (2) Significant (adj.p < 0.05) in all 4 datasets
#     (3) Consistent direction in all 4 datasets
#     (4) Independent AUC >= 0.80 in GSE135065
#
# - "Moderate":
#     Strict WGCNA Hub OR (Significant in all 4 with consistent direction AND independent AUC >= 0.80),
#     but NOT both.
#
# - "Weak/not supported":
#     Fails both conditions (i.e., not a strict hub, and fails either 4-cohort significance,
#     direction consistency, or independent AUC >= 0.80).
# ------------------------------------------------------------------------------

set.seed(20261102)
options(stringsAsFactors = FALSE)

LOG_FILE <- "logs/pre_specified_test.log"
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("\n=== Initializing PART 5: Final Synthesis Table Assembly ===")

# Load previous outputs
wgcna_df <- read.csv("results/tables/prespecified/wgcna_status_15.csv")
auc_df <- read.csv("results/tables/prespecified/single_gene_auc_15.csv")
feat_df <- read.csv("results/tables/prespecified/feature_importance_15.csv")

prespecified_15 <- c(
  "ASPN", "COL14A1", "COL15A1", "COL1A1", "COL3A1", "COL5A2", "COMP",
  "IL6", "MMP7", "MMP8", "MMP10", "OGN", "SFRP4", "TGFB3", "TIMP3"
)

# Build synthesis rows
final_rows <- list()

for (g in prespecified_15) {
  w_row <- wgcna_df[wgcna_df$gene == g, ]
  f_row <- feat_df[feat_df$gene == g, ]
  
  # Cohort-specific stats from auc_df
  auc_32 <- auc_df[auc_df$gene == g & auc_df$dataset == "GSE32539", ]
  auc_110 <- auc_df[auc_df$gene == g & auc_df$dataset == "GSE110147", ]
  auc_24 <- auc_df[auc_df$gene == g & auc_df$dataset == "GSE24206", ]
  auc_135 <- auc_df[auc_df$gene == g & auc_df$dataset == "GSE135065", ]
  
  is_strict_hub <- (w_row$verdict_strict == "Hub (strict)")
  sig_all4 <- (auc_32$significant_all4 == "yes")
  dir_consistent <- (auc_32$direction_consistent_all4 == "yes")
  indep_auc_ge_80 <- (auc_135$AUC >= 0.80)
  
  condition_a <- is_strict_hub
  condition_b <- (sig_all4 && dir_consistent && indep_auc_ge_80)
  
  if (condition_a && condition_b) {
    tier <- "Strong"
    reason <- "Strict WGCNA hub + Significant & concordant in all 4 datasets + Independent AUC >= 0.80"
  } else if (condition_a || condition_b) {
    tier <- "Moderate"
    if (condition_a && !condition_b) {
      fail_parts <- c()
      if (!sig_all4) fail_parts <- c(fail_parts, "not significant in all 4 datasets")
      if (!dir_consistent) fail_parts <- c(fail_parts, "inconsistent direction in GSE135065")
      if (!indep_auc_ge_80) fail_parts <- c(fail_parts, paste0("independent AUC < 0.80 (", auc_135$AUC, ")"))
      reason <- paste0("Strict WGCNA hub, but ", paste(fail_parts, collapse = ", "))
    } else {
      reason <- "Valid in all 4 cohorts with AUC >= 0.80, but failed strict WGCNA hub cutoff"
    }
  } else {
    tier <- "Weak/not supported"
    fail_parts <- c()
    if (!is_strict_hub) fail_parts <- c(fail_parts, paste0("failed strict hub (", w_row$verdict_strict, ")"))
    if (!sig_all4) fail_parts <- c(fail_parts, "not significant in all 4 datasets")
    if (!dir_consistent) fail_parts <- c(fail_parts, "direction inconsistent")
    if (!indep_auc_ge_80) fail_parts <- c(fail_parts, paste0("independent AUC < 0.80 (", auc_135$AUC, ")"))
    reason <- paste(fail_parts, collapse = "; ")
  }
  
  final_rows[[g]] <- data.frame(
    gene = g,
    WGCNA_module = w_row$module_color,
    WGCNA_is_key_module = w_row$is_key_module,
    WGCNA_MM = w_row$MM,
    WGCNA_GS = w_row$GS,
    WGCNA_kWithin_rank = w_row$rank_kWithin_module,
    WGCNA_kWithin_pct = w_row$pct_kWithin_module,
    WGCNA_strict_verdict = w_row$verdict_strict,
    WGCNA_relaxed_verdict = w_row$verdict_exploratory,
    WGCNA_bootstrap_frequency = w_row$bootstrap_hub_frequency,
    logFC_GSE32539 = auc_32$logFC,
    AUC_GSE32539 = auc_32$AUC,
    logFC_GSE110147 = auc_110$logFC,
    AUC_GSE110147 = auc_110$AUC,
    logFC_GSE24206 = auc_24$logFC,
    AUC_GSE24206 = auc_24$AUC,
    logFC_GSE135065_indep = auc_135$logFC,
    AUC_GSE135065_indep = auc_135$AUC,
    Wilcoxon_adj_p_GSE135065 = auc_135$Wilcoxon_adj_p,
    direction_consistent_all4 = ifelse(dir_consistent, "yes", "no"),
    significant_all4 = ifelse(sig_all4, "yes", "no"),
    ML_importance_rank = f_row$Mean_Importance_Rank,
    final_evidence_tier = tier,
    tier_justification = reason
  )
}

final_df <- do.call(rbind, final_rows)
final_df <- final_df[order(match(final_df$final_evidence_tier, c("Strong", "Moderate", "Weak/not supported")), final_df$ML_importance_rank), ]

write.csv(final_df, "results/tables/prespecified/final_15_gene_table.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/final_15_gene_table.csv")

log_msg("\nFinal Evidence Tier Breakdown:")
print(table(final_df$final_evidence_tier))
print(final_df[, c("gene", "WGCNA_strict_verdict", "AUC_GSE135065_indep", "direction_consistent_all4", "final_evidence_tier")])

log_msg("=== PART 5 Assembly Complete Successfully ===")

