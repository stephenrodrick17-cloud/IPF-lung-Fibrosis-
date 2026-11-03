#!/usr/bin/env Rscript
## ==============================================================================
## 08_prespecified_wgcna.R
## PART 1: WGCNA Evaluation of the Pre-Specified 15-Gene List (from Thesis)
## Seed: 20261102 | R version: 4.6.2 (Windows)
## ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0: Fixed Criteria Defined in Advance
# ------------------------------------------------------------------------------
# Strict Hub Criterion:
# 1. Assigned to a Key Module (|r| >= 0.50 AND adj.p < 0.05 vs IPF trait)
# 2. Module Membership |MM| >= 0.80
# 3. Gene Significance |GS| >= 0.20
#
# Mutually Exclusive Strict Verdicts:
# - "Hub (strict)": Key module AND |MM| >= 0.8 AND |GS| >= 0.2
# - "Key module, below hub cutoff": Key module, but fails MM, fails GS, or fails both
# - "Non-key module": Module |r| < 0.50 or adj.p >= 0.05
# - "Grey / unassigned": Module 0 (grey)
# - "Not in network": Gene not in the top 8,000 MAD network
#
# Exploratory Relaxed Verdict (Reported separately, does not alter strict verdict):
# - "Hub (relaxed)": Key module AND (|MM| >= 0.70 OR kWithin in top 20% of its module)
# ------------------------------------------------------------------------------

set.seed(20261102)
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})
WGCNA::disableWGCNAThreads()
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)

# Directories
dirs <- c("logs", "results/tables/prespecified", "results/figures/prespecified")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- "logs/pre_specified_test.log"
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("=== Initializing PART 1: Pre-Specified 15-Gene WGCNA Testing ===")
log_msg("Random seed set to 20261102")

# Save session info
sink("logs/sessionInfo_prespecified.txt")
cat("=== Pre-Specified 15-Gene Testing Session Info ===\nSeed: 20261102\n\n")
print(sessionInfo())
sink()

save_plot <- function(p, file_base, width = 8, height = 6) {
  ggsave(paste0(file_base, ".png"), plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

# Pre-specified 15 genes
prespecified_15 <- c(
  "ASPN", "COL14A1", "COL15A1", "COL1A1", "COL3A1", "COL5A2", "COMP",
  "IL6", "MMP7", "MMP8", "MMP10", "OGN", "SFRP4", "TGFB3", "TIMP3"
)
log_msg("Pre-specified gene list: ", paste(prespecified_15, collapse = ", "))

# ------------------------------------------------------------------------------
# STEP 1: Load Existing WGCNA Artifacts
# ------------------------------------------------------------------------------
log_msg("\n--- Loading existing WGCNA results from GSE32539 ---")
mod_assign <- read.csv("results/tables/wgcna/module_assignment.csv")
mod_trait <- read.csv("results/tables/wgcna/module_trait.csv")
gene_mm_gs <- read.csv("results/tables/wgcna/gene_MM_GS.csv")

# Full expression data and labels for GSE32539
e32_full <- read.csv("data/processed/GSE32539_expr_gene_level.csv", check.names = FALSE)
rownames(e32_full) <- e32_full$gene
e32_full$gene <- NULL

l32 <- read.csv("data/processed/ml_input/GSE32539_labels.csv")
rownames(l32) <- l32$sample_id
keep32 <- setdiff(l32$sample_id, "GSM806411")
l32 <- l32[keep32, ]

net_genes <- mod_assign$gene
datExpr32 <- t(as.matrix(e32_full[net_genes, keep32]))
trait_vec <- ifelse(l32$label == "IPF", 1, 0)

key_modules <- mod_trait$color[mod_trait$is_key_module]
log_msg("Key modules (|r| >= 0.50, adj.p < 0.05): ", paste(key_modules, collapse = ", "))

# ------------------------------------------------------------------------------
# STEP 2: Characterize Each of the 15 Genes in GSE32539
# ------------------------------------------------------------------------------
log_msg("\n--- Characterizing 15 genes in the GSE32539 network ---")

# Rank in MAD across all genes in GSE32539
all_mads <- apply(e32_full[, keep32], 1, mad)
all_mad_ranks <- rank(-all_mads)

mod_trait_r_map <- setNames(mod_trait$cor_trait, mod_trait$color)
mod_trait_p_map <- setNames(mod_trait$adj_p_value, mod_trait$color)
mod_sizes <- table(mod_assign$module_color)

status_list <- list()

for (g in prespecified_15) {
  in_net <- g %in% net_genes
  mad_r <- all_mad_ranks[g]
  
  if (in_net) {
    row_g <- gene_mm_gs[gene_mm_gs$gene == g, ]
    m_col <- row_g$module_color
    m_num <- row_g$module_number
    m_r <- mod_trait_r_map[m_col]
    m_adj_p <- mod_trait_p_map[m_col]
    is_key <- m_col %in% key_modules
    
    mm_val <- row_g$MM
    mm_pval <- row_g$p_MM
    gs_val <- row_g$GS
    gs_pval <- row_g$p_GS
    k_w <- row_g$kWithin
    k_w_rank <- row_g$rank_in_module
    m_size <- as.numeric(mod_sizes[m_col])
    pct_kw <- round((1 - (k_w_rank - 1) / m_size) * 100, 2)
    
    k_tot_rank <- rank(-gene_mm_gs$kTotal)[which(gene_mm_gs$gene == g)]
    pct_ktot <- round((1 - (k_tot_rank - 1) / nrow(gene_mm_gs)) * 100, 2)
    
    # Strict verdict
    if (m_col == "grey") {
      verdict_strict <- "Grey / unassigned"
    } else if (!is_key) {
      verdict_strict <- "Non-key module"
    } else if (abs(mm_val) >= 0.80 && abs(gs_val) >= 0.20) {
      verdict_strict <- "Hub (strict)"
    } else if (abs(mm_val) < 0.80 && abs(gs_val) >= 0.20) {
      verdict_strict <- "Key module, below hub cutoff (fails MM)"
    } else if (abs(mm_val) >= 0.80 && abs(gs_val) < 0.20) {
      verdict_strict <- "Key module, below hub cutoff (fails GS)"
    } else {
      verdict_strict <- "Key module, below hub cutoff (fails both)"
    }
    
    # Exploratory relaxed verdict
    if (verdict_strict == "Hub (strict)") {
      verdict_exp <- "Hub (strict)"
    } else if (is_key && (abs(mm_val) >= 0.70 || pct_kw >= 80)) {
      verdict_exp <- "Hub (relaxed)"
    } else {
      verdict_exp <- "Non-hub"
    }
    
  } else {
    m_col <- "Not in network"
    m_num <- NA
    m_r <- NA
    m_adj_p <- NA
    is_key <- FALSE
    mm_val <- NA
    mm_pval <- NA
    gs_val <- NA
    gs_pval <- NA
    k_w <- NA
    k_w_rank <- NA
    pct_kw <- NA
    pct_ktot <- NA
    verdict_strict <- "Not in network"
    verdict_exp <- "Not in network"
  }
  
  status_list[[g]] <- data.frame(
    gene = g,
    in_network = ifelse(in_net, "yes", "no"),
    mad_rank = as.integer(mad_r),
    module_number = m_num,
    module_color = m_col,
    module_r = round(as.numeric(m_r), 4),
    module_adj_p = signif(as.numeric(m_adj_p), 4),
    is_key_module = ifelse(is_key, "yes", "no"),
    MM = round(as.numeric(mm_val), 4),
    p_MM = signif(as.numeric(mm_pval), 4),
    GS = round(as.numeric(gs_val), 4),
    p_GS = signif(as.numeric(gs_pval), 4),
    kWithin = round(as.numeric(k_w), 3),
    rank_kWithin_module = as.integer(k_w_rank),
    pct_kWithin_module = pct_kw,
    pct_kTotal_network = pct_ktot,
    verdict_strict = verdict_strict,
    verdict_exploratory = verdict_exp
  )
}

status_df <- do.call(rbind, status_list)

# ------------------------------------------------------------------------------
# STEP 3: Robustness Analysis (Bootstrap & Soft Power Sensitivity)
# ------------------------------------------------------------------------------
log_msg("\n--- Robustness Analysis: 100 Bootstrap Resampling ---")
# Resample samples 100 times, recompute MEs, recompute MM and GS for the 15 genes
B <- 100
n_samples <- nrow(datExpr32)
boot_hub_count <- setNames(rep(0, length(prespecified_15)), prespecified_15)

# Original module colors vector for all 8000 genes
colors_vec <- mod_assign$module_color

for (b in 1:B) {
  boot_idx <- sample(n_samples, replace = TRUE)
  dat_boot <- datExpr32[boot_idx, ]
  trait_boot <- trait_vec[boot_idx]
  
  # Recompute MEs
  me_boot <- moduleEigengenes(dat_boot, colors = colors_vec)$eigengenes
  colnames(me_boot) <- gsub("^ME", "", colnames(me_boot))
  
  for (g in prespecified_15) {
    m_col <- status_df$module_color[status_df$gene == g]
    if (m_col %in% colnames(me_boot) && status_df$is_key_module[status_df$gene == g] == "yes") {
      expr_g <- dat_boot[, g]
      # MM
      mm_b <- cor(expr_g, me_boot[, m_col], use = "p")
      # GS
      gs_b <- cor(expr_g, trait_boot, use = "p")
      if (abs(mm_b) >= 0.80 && abs(gs_b) >= 0.20) {
        boot_hub_count[g] <- boot_hub_count[g] + 1
      }
    }
  }
}

status_df$bootstrap_hub_frequency <- round(boot_hub_count[status_df$gene] / B, 3)
log_msg("Bootstrap hub stability frequencies computed across 100 resamplings.")

# Soft Power Sensitivity (Power 6 and Power 10)
log_msg("\n--- Soft Power Sensitivity Analysis (beta = 6 and beta = 10) ---")
# We test if intramodular connectivity and MM maintain the genes within the top tier at powers 6 and 10
# Calculate adjacency at power 6 and power 10 for the modules of the 15 genes
check_power_module_stability <- function(pwr) {
  # Test module preservation / connectivity consistency
  pwr_res <- character(length(prespecified_15))
  names(pwr_res) <- prespecified_15
  
  # Recompute adjacency on each gene's module genes
  for (g in prespecified_15) {
    m <- status_df$module_color[status_df$gene == g]
    m_genes <- mod_assign$gene[mod_assign$module_color == m]
    adj_m <- (0.5 * (1 + cor(datExpr32[, m_genes], use = "p")))^pwr
    k_w_pwr <- rowSums(adj_m) - 1
    rank_pwr <- rank(-k_w_pwr)[g]
    top20_pct <- rank_pwr <= (0.20 * length(m_genes))
    pwr_res[g] <- ifelse(top20_pct, "Top20%_preserved", "Moderate_connectivity")
  }
  return(pwr_res)
}

stab_p6 <- check_power_module_stability(6)
stab_p10 <- check_power_module_stability(10)

status_df$power_robustness_beta6 <- stab_p6[status_df$gene]
status_df$power_robustness_beta10 <- stab_p10[status_df$gene]
status_df$power_robustness <- ifelse(
  status_df$power_robustness_beta6 == "Top20%_preserved" & status_df$power_robustness_beta10 == "Top20%_preserved",
  "High_stability", "Moderate_stability"
)

# ------------------------------------------------------------------------------
# STEP 4: Replication in GSE110147 (N=33)
# ------------------------------------------------------------------------------
log_msg("\n--- Replicating MM and GS in GSE110147 (N=33; low power) ---")
e110_full <- read.csv("data/processed/GSE110147_expr_gene_level.csv", check.names = FALSE)
rownames(e110_full) <- e110_full$gene
e110_full$gene <- NULL

l110 <- read.csv("data/processed/ml_input/GSE110147_labels.csv")
trait110 <- ifelse(l110$label == "IPF", 1, 0)

# Overlapping genes
common_net110 <- intersect(net_genes, rownames(e110_full))
datExpr110 <- t(as.matrix(e110_full[common_net110, l110$sample_id]))
colors110 <- mod_assign$module_color[match(common_net110, mod_assign$gene)]

# Project GSE32539 module eigengenes in GSE110147
me110 <- moduleEigengenes(datExpr110, colors = colors110)$eigengenes
colnames(me110) <- gsub("^ME", "", colnames(me110))

mm110_list <- numeric(length(prespecified_15))
names(mm110_list) <- prespecified_15
gs110_list <- numeric(length(prespecified_15))
names(gs110_list) <- prespecified_15

for (g in prespecified_15) {
  if (g %in% rownames(e110_full)) {
    expr_g110 <- as.numeric(e110_full[g, l110$sample_id])
    m_col <- status_df$module_color[status_df$gene == g]
    if (m_col %in% colnames(me110)) {
      mm110_list[g] <- round(cor(expr_g110, me110[, m_col], use = "p"), 4)
    } else {
      mm110_list[g] <- NA
    }
    gs110_list[g] <- round(cor(expr_g110, trait110, use = "p"), 4)
  } else {
    mm110_list[g] <- NA
    gs110_list[g] <- NA
  }
}

status_df$MM_GSE110147 <- mm110_list[status_df$gene]
status_df$GS_GSE110147 <- gs110_list[status_df$gene]
status_df$replication_note_GSE110147 <- "Low_powered (N=33); evaluated for direction concordance"

# Save Table 1: wgcna_status_15.csv
write.csv(status_df, "results/tables/prespecified/wgcna_status_15.csv", row.names = FALSE)
log_msg("Saved results/tables/prespecified/wgcna_status_15.csv")

log_msg("\nWGCNA Strict Verdicts Summary:")
print(table(status_df$verdict_strict))

# ------------------------------------------------------------------------------
# STEP 5: Visualizations
# ------------------------------------------------------------------------------
log_msg("\n--- Generating PART 1 Figures ---")

# 5.1 MM vs GS Scatter with 15 Genes Highlighted for green and turquoise modules
plot_mm_gs_with_15 <- function(m_color) {
  mod_df <- gene_mm_gs[gene_mm_gs$module_color == m_color, ]
  mod_15 <- mod_df[mod_df$gene %in% prespecified_15, ]
  
  p <- ggplot(mod_df, aes(x = abs(MM), y = abs(GS))) +
    geom_point(color = "gray75", alpha = 0.6, size = 2) +
    geom_point(data = mod_15, aes(x = abs(MM), y = abs(GS)), color = "#e31a1c", size = 4) +
    geom_text(data = mod_15, aes(label = gene), color = "black", fontface = "bold", vjust = -0.8, size = 3.5) +
    geom_vline(xintercept = 0.8, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "blue") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Module Membership vs Gene Significance: ", toupper(m_color), " Module"),
      subtitle = "Red = Pre-specified thesis genes; Dashed blue = Hub thresholds (|MM| >= 0.8, |GS| >= 0.2)",
      x = "|Module Membership (MM)|",
      y = "|Gene Significance (GS)|"
    )
  return(p)
}

p_sc_green <- plot_mm_gs_with_15("green")
save_plot(p_sc_green, "results/figures/prespecified/MM_vs_GS_green_15genes", width = 8, height = 6)

p_sc_turq <- plot_mm_gs_with_15("turquoise")
save_plot(p_sc_turq, "results/figures/prespecified/MM_vs_GS_turquoise_15genes", width = 8, height = 6)

# 5.2 Rank plot of kWithin within green module with 15 genes marked
green_genes_df <- gene_mm_gs[gene_mm_gs$module_color == "green", ]
green_genes_df <- green_genes_df[order(green_genes_df$rank_in_module), ]
green_genes_df$is_prespecified <- green_genes_df$gene %in% prespecified_15

p_rank_k <- ggplot(green_genes_df, aes(x = rank_in_module, y = kWithin)) +
  geom_line(color = "darkgreen", lwd = 1) +
  geom_point(data = subset(green_genes_df, is_prespecified), color = "red", size = 3.5) +
  geom_text(data = subset(green_genes_df, is_prespecified), aes(label = gene), 
            vjust = -0.7, fontface = "bold", size = 3.2, check_overlap = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Intramodular Connectivity (kWithin) Distribution in GREEN Module",
    subtitle = "Red points = Pre-specified candidate genes (n = 10 in green module)",
    x = "Intramodular Rank (1 = Highest connectivity)",
    y = "Intramodular Connectivity (kWithin)"
  )
save_plot(p_rank_k, "results/figures/prespecified/kWithin_rank_green_15genes", width = 9, height = 6)

# 5.3 Heatmap of the 15 Genes (Module, MM, GS, kWithin, Bootstrap Hub Freq)
plot_df <- status_df
plot_df$gene <- factor(plot_df$gene, levels = rev(plot_df$gene))

p_heat1 <- ggplot(plot_df, aes(x = "Module Membership", y = gene, fill = MM)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", MM)), size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p_heat2 <- ggplot(plot_df, aes(x = "Gene Significance", y = gene, fill = GS)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", GS)), size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())

p_heat3 <- ggplot(plot_df, aes(x = "kWithin Percentile", y = gene, fill = pct_kWithin_module)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(pct_kWithin_module), "%")), size = 3.5) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 100)) +
  theme_bw(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())

p_heat4 <- ggplot(plot_df, aes(x = "Bootstrap Hub Freq", y = gene, fill = bootstrap_hub_frequency)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", bootstrap_hub_frequency)), size = 3.5) +
  scale_fill_gradient(low = "gray95", high = "darkgreen", limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())

p_heat_combo <- p_heat1 + p_heat2 + p_heat3 + p_heat4 + plot_layout(ncol = 4, widths = c(1.3, 1, 1, 1)) +
  plot_annotation(title = "WGCNA Characteristic Heatmap: 15 Pre-Specified Candidate Genes")

save_plot(p_heat_combo, "results/figures/prespecified/heatmap_15genes_wgcna", width = 10, height = 7)

log_msg("\n=== PART 1 Complete Successfully ===")

