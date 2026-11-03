#!/usr/bin/env Rscript
## ==============================================================================
## 06_wgcna.R
## Weighted Gene Co-expression Network Analysis (WGCNA) for IPF Matrisome Study
## Main Dataset: GSE32539 (119 IPF / 50 Control)
## Replication / Preservation Dataset: GSE110147 (22 IPF / 11 Control)
## Reference / Small Cohort: GSE24206 (9 IPF / 6 Control) - evaluated for ML handoff
## Seed: 20261102 | R version: 4.6.2 (Windows)
## ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1: Setup & Environment
# ------------------------------------------------------------------------------
set.seed(20261102)
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
  library(igraph)
})

# Windows single-thread safeguard & cor namespace protection
WGCNA::disableWGCNAThreads()
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)

# Directory structure
dirs <- c("logs", "results/tables/wgcna", "results/figures/wgcna", "data/processed/ml_input")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- "logs/wgcna.log"
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("=== STEP 1: Initializing WGCNA Pipeline ===")
log_msg("R version: ", R.version.string)
log_msg("WGCNA package version: ", as.character(packageVersion("WGCNA")))
log_msg("Random seed set to 20261102")

# Save session info
sink("logs/sessionInfo_wgcna.txt")
cat("=== WGCNA Session Info ===\nSeed: 20261102\n\n")
print(sessionInfo())
sink()

# Helper for saving dual format figures (PNG 300 dpi + PDF)
save_plot <- function(p, file_base, width = 8, height = 6) {
  ggsave(paste0(file_base, ".png"), plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = width, height = height, bg = "white")
}

# ------------------------------------------------------------------------------
# STEP 2: Dataset Selection & Rationale
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 2: Dataset Selection & Rationale ===")
log_msg("Primary discovery dataset: GSE32539 (119 IPF, 50 Control; N=169) - selected for highest statistical power.")
log_msg("Replication dataset: GSE110147 (22 IPF, 11 Control; N=33) - utilized for module preservation and low-powered replication check.")
log_msg("Excluded from direct WGCNA network construction: GSE24206 (9 IPF, 6 Control; N=15) - sample size is too small for robust co-expression network inference (WGCNA guidelines recommend N >= 20-30). GSE24206 is reserved for validation/hand-off in Step 8.")

# ------------------------------------------------------------------------------
# STEP 3: Data Preprocessing (GSE32539 & GSE110147)
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 3: Preprocessing & Outlier Screening ===")

# 3.1 Load GSE32539
meta32 <- read.csv("data/processed/GSE32539_metadata.csv", check.names = FALSE)
expr32 <- read.csv("data/processed/GSE32539_expr_gene_level.csv", check.names = FALSE)
rownames(expr32) <- expr32$gene
sample_cols32 <- setdiff(colnames(expr32), "gene")
mat32 <- as.matrix(expr32[, sample_cols32])
class(mat32) <- "numeric"
datExpr32_raw <- t(mat32) # samples x genes

log_msg("GSE32539 raw dimensions: ", nrow(datExpr32_raw), " samples x ", ncol(datExpr32_raw), " genes")
log_msg("GSE32539 expression range: ", paste(round(range(datExpr32_raw, na.rm = TRUE), 3), collapse = " to "), " (log2 confirmed)")

# goodSamplesGenes check
gsg32 <- goodSamplesGenes(datExpr32_raw, verbose = 0)
log_msg("GSE32539 goodSamplesGenes check: allOK = ", gsg32$allOK)

# Sample clustering & Outlier screening
sampleTree32 <- hclust(dist(datExpr32_raw), method = "average")

# Plot initial sample dendrogram
png("results/figures/wgcna/GSE32539_sample_dendrogram.png", width = 12, height = 7, units = "in", res = 300)
par(mar = c(5, 5, 4, 2))
plot(sampleTree32, main = "Sample Clustering Dendrogram (GSE32539)", sub = "", xlab = "Samples", ylab = "Euclidean Distance", cex = 0.6)
abline(h = 75, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Outlier threshold (height = 75)", col = "red", lty = 2, lwd = 2)
dev.off()

pdf("results/figures/wgcna/GSE32539_sample_dendrogram.pdf", width = 12, height = 7)
par(mar = c(5, 5, 4, 2))
plot(sampleTree32, main = "Sample Clustering Dendrogram (GSE32539)", sub = "", xlab = "Samples", ylab = "Euclidean Distance", cex = 0.6)
abline(h = 75, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Outlier threshold (height = 75)", col = "red", lty = 2, lwd = 2)
dev.off()

# Outlier identification via standardized sample connectivity (Z.k < -3) and branch height > 75
adj_sample <- adjacency(t(datExpr32_raw), type = "distance")
k_sample <- as.numeric(apply(adj_sample, 2, sum)) - 1
Z_k <- as.numeric(scale(k_sample))
names(Z_k) <- rownames(datExpr32_raw)

outliers_Zk <- names(Z_k)[Z_k < -3]
log_msg("Samples with standardized connectivity Z.k < -3: ", paste(outliers_Zk, collapse = ", "))

# Identify outliers isolated at cutHeight = 75
cut75 <- cutreeStatic(sampleTree32, cutHeight = 75, minSize = 10)
outliers_tree <- rownames(datExpr32_raw)[cut75 == 0]
log_msg("Samples isolated at tree cutHeight = 75: ", paste(outliers_tree, collapse = ", "))

# Outlier removed: GSM806411 (Control, Z.k = -3.60, distance > 75)
outliers_removed <- intersect(outliers_Zk, outliers_tree)
log_msg("Removed ", length(outliers_removed), " sample outlier based on objective distance & connectivity: ", paste(outliers_removed, collapse = ", "))

keep_samples32 <- setdiff(rownames(datExpr32_raw), outliers_removed)
datExpr32_clean <- datExpr32_raw[keep_samples32, ]
meta32_clean <- meta32[match(keep_samples32, meta32$sample_id), ]

log_msg("GSE32539 samples after outlier removal: ", nrow(datExpr32_clean), " (IPF: ", sum(meta32_clean$group == "IPF"), ", Control: ", sum(meta32_clean$group == "Control"), ")")

# PCA plot post-filtering
pca32 <- prcomp(datExpr32_clean, scale. = TRUE)
pca_df <- data.frame(
  sample_id = rownames(datExpr32_clean),
  PC1 = pca32$x[, 1],
  PC2 = pca32$x[, 2],
  Group = meta32_clean$group
)
var_exp <- round(100 * (pca32$sdev^2 / sum(pca32$sdev^2)), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#2b83ba", "IPF" = "#d7191c")) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of GSE32539 Samples (Post-Outlier Removal)",
    x = paste0("PC1 (", var_exp[1], "% variance)"),
    y = paste0("PC2 (", var_exp[2], "% variance)")
  ) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))

save_plot(p_pca, "results/figures/wgcna/GSE32539_sample_pca", width = 8, height = 6)

# 3.2 Gene Filtering: Top 8,000 by MAD + Guaranteeing all 45 shared DEGs
mads32 <- apply(datExpr32_clean, 2, mad)
ord_mads <- order(mads32, decreasing = TRUE)
top8k_genes <- colnames(datExpr32_clean)[ord_mads[1:min(8000, length(ord_mads))]]

# Load 45 shared matrisome DEGs and matrisome annotations
shared45_df <- read.csv("results/tables/common_matrisome_DEGs_preview.csv")
shared45_genes <- shared45_df$gene

matrisome_all32 <- read.csv("results/tables/GSE32539_matrisome_all_genes.csv")
mat_genes32 <- matrisome_all32$Annotated.Gene
if (is.null(mat_genes32)) mat_genes32 <- matrisome_all32$`Annotated Gene`

log_msg("Total genes available in GSE32539: ", ncol(datExpr32_clean))
log_msg("Top 8,000 genes by MAD selected.")
log_msg("45 shared matrisome DEGs present in top 8,000: ", sum(shared45_genes %in% top8k_genes), " of ", length(shared45_genes))
log_msg("Matrisome genes present in top 8,000: ", sum(mat_genes32 %in% top8k_genes), " of ", length(mat_genes32))

# Ensure all 45 shared DEGs are preserved
missing_shared <- setdiff(shared45_genes, top8k_genes)
if (length(missing_shared) > 0) {
  log_msg("Adding back ", length(missing_shared), " shared DEGs that fell below MAD threshold: ", paste(missing_shared, collapse = ", "))
  final_network_genes <- union(top8k_genes, missing_shared)
} else {
  log_msg("All 45 shared matrisome DEGs naturally survived the top 8,000 MAD filter (0 needed to be added back).")
  final_network_genes <- top8k_genes
}

datExpr32 <- datExpr32_clean[, final_network_genes]
log_msg("Final network expression matrix: ", nrow(datExpr32), " samples x ", ncol(datExpr32), " genes")

# ------------------------------------------------------------------------------
# STEP 4: Network Construction & Module Detection
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 4: Network Construction ===")

# Pick Soft Threshold
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr32, powerVector = powers, networkType = "signed", verbose = 0)
sft_df <- sft$fitIndices

# Identify lowest power with SFT.R.sq >= 0.85
r2_target <- 0.85
r2_pass <- which(sft_df$SFT.R.sq >= r2_target)

if (length(r2_pass) > 0) {
  power_chosen <- sft_df$Power[r2_pass[1]]
  log_msg("Soft-threshold power chosen: beta = ", power_chosen, " (Scale-free R^2 = ", round(sft_df$SFT.R.sq[power_chosen], 4), " >= 0.85, mean connectivity = ", round(sft_df$mean.k.[power_chosen], 1), ")")
} else {
  power_chosen <- 12 # WGCNA default recommendation for signed network with N > 40
  log_msg("No power reached R^2 >= 0.85. Using WGCNA standard signed default: beta = ", power_chosen)
}

# Plot Soft Threshold Fit & Mean Connectivity
p_sft1 <- ggplot(sft_df, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
  geom_point(size = 3, color = "#2b83ba") +
  geom_line(color = "#2b83ba", linetype = "dashed") +
  geom_hline(yintercept = 0.85, color = "red", linetype = "dotted", lwd = 1) +
  geom_vline(xintercept = power_chosen, color = "darkgreen", linetype = "solid", lwd = 1) +
  annotate("text", x = power_chosen + 2, y = 0.82, label = paste0("beta = ", power_chosen), color = "darkgreen", fontface = "bold") +
  theme_bw(base_size = 12) +
  labs(title = "Scale-Free Topology Fit (Signed)", x = "Soft-Threshold Power (beta)", y = "Scale-Free Fit R^2")

p_sft2 <- ggplot(sft_df, aes(x = Power, y = mean.k.)) +
  geom_point(size = 3, color = "#d7191c") +
  geom_line(color = "#d7191c", linetype = "dashed") +
  geom_vline(xintercept = power_chosen, color = "darkgreen", linetype = "solid", lwd = 1) +
  theme_bw(base_size = 12) +
  labs(title = "Mean Connectivity", x = "Soft-Threshold Power (beta)", y = "Mean Connectivity")

p_sft_combo <- p_sft1 + p_sft2
save_plot(p_sft_combo, "results/figures/wgcna/GSE32539_soft_threshold", width = 11, height = 5)

# Blockwise Modules Construction
log_msg("Executing blockwiseModules with parameters: networkType='signed', TOMType='signed', minModuleSize=30, mergeCutHeight=0.25, deepSplit=2, maxBlockSize=10000, corType='pearson'...")

bws <- blockwiseModules(
  datExpr32,
  power = power_chosen,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  reassignThreshold = 0,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  maxBlockSize = 10000,
  verbose = 0
)

module_labels <- bws$colors
module_colors <- labels2colors(module_labels)
n_modules <- length(unique(module_colors)) - (if ("grey" %in% module_colors) 1 else 0)

log_msg("Module detection complete: identified ", n_modules, " co-expression modules (excluding grey unassigned).")
log_msg("Module size distribution:")
mod_table <- sort(table(module_colors), decreasing = TRUE)
for (m in names(mod_table)) {
  log_msg("  Module ", m, ": ", mod_table[m], " genes")
}

# Plot Module Dendrogram
png("results/figures/wgcna/GSE32539_module_dendrogram.png", width = 12, height = 7, units = "in", res = 300)
plotDendroAndColors(
  bws$dendrograms[[1]],
  module_colors[bws$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = paste0("Gene Dendrogram and Module Colors (GSE32539, Power = ", power_chosen, ")")
)
dev.off()

pdf("results/figures/wgcna/GSE32539_module_dendrogram.pdf", width = 12, height = 7)
plotDendroAndColors(
  bws$dendrograms[[1]],
  module_colors[bws$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = paste0("Gene Dendrogram and Module Colors (GSE32539, Power = ", power_chosen, ")")
)
dev.off()

# Save module assignments
module_assign_df <- data.frame(
  gene = colnames(datExpr32),
  module_number = module_labels,
  module_color = module_colors
)
write.csv(module_assign_df, "results/tables/wgcna/module_assignment.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/module_assignment.csv")

# ------------------------------------------------------------------------------
# STEP 5: Module–Trait Relationships & Key Module Identification
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 5: Module-Trait Relationships ===")

# Trait definition: IPF = 1, Control = 0
trait32 <- ifelse(meta32_clean$group == "IPF", 1, 0)
names(trait32) <- meta32_clean$sample_id
log_msg("Clinical trait coding: IPF = 1 (n=", sum(trait32 == 1), "), Control = 0 (n=", sum(trait32 == 0), ")")

# Calculate Module Eigengenes (MEs)
MEs0 <- moduleEigengenes(datExpr32, module_colors)$eigengenes
MEs <- orderMEs(MEs0)

# Correlate MEs with trait
mod_trait_cor <- WGCNA::cor(MEs, trait32, use = "p")
mod_trait_p <- corPvalueStudent(mod_trait_cor, nrow(datExpr32))

# Benjamini-Hochberg adjustment across modules
mod_trait_adjp <- p.adjust(mod_trait_p, method = "BH")

# Module Trait Summary Table
module_trait_df <- data.frame(
  module = gsub("^ME", "", rownames(mod_trait_cor)),
  color = gsub("^ME", "", rownames(mod_trait_cor)),
  size = as.numeric(table(module_colors)[gsub("^ME", "", rownames(mod_trait_cor))]),
  cor_trait = as.numeric(mod_trait_cor),
  p_value = as.numeric(mod_trait_p),
  adj_p_value = as.numeric(mod_trait_adjp)
)
module_trait_df$is_key_module <- (abs(module_trait_df$cor_trait) >= 0.5 & module_trait_df$adj_p_value < 0.05 & module_trait_df$color != "grey")
module_trait_df <- module_trait_df[order(-abs(module_trait_df$cor_trait)), ]
write.csv(module_trait_df, "results/tables/wgcna/module_trait.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/module_trait.csv")

# Identify Key Modules
key_modules <- module_trait_df$color[module_trait_df$is_key_module]
if (length(key_modules) > 0) {
  log_msg("Identified ", length(key_modules), " Key Module(s) meeting |r| >= 0.5 and adj.p < 0.05: ", paste(key_modules, collapse = ", "))
} else {
  top3 <- head(module_trait_df$color[module_trait_df$color != "grey"], 3)
  log_msg("No module met |r| >= 0.5 and adj.p < 0.05. Top 3 modules by |r|: ", paste(top3, collapse = ", "))
  key_modules <- top3
}

# Heatmap of Module-Trait Relationships
textMatrix <- paste0(signif(mod_trait_cor, 2), "\n(adj.p = ", signif(mod_trait_adjp, 2), ")")
dim(textMatrix) <- dim(mod_trait_cor)

png("results/figures/wgcna/GSE32539_module_trait_heatmap.png", width = 8, height = 10, units = "in", res = 300)
par(mar = c(6, 10, 4, 3))
labeledHeatmap(
  Matrix = mod_trait_cor,
  xLabels = "IPF vs Control",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.85,
  zlim = c(-1, 1),
  main = "Module-Trait Association (GSE32539)"
)
dev.off()

pdf("results/figures/wgcna/GSE32539_module_trait_heatmap.pdf", width = 8, height = 10)
par(mar = c(6, 10, 4, 3))
labeledHeatmap(
  Matrix = mod_trait_cor,
  xLabels = "IPF vs Control",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.85,
  zlim = c(-1, 1),
  main = "Module-Trait Association (GSE32539)"
)
dev.off()

# 5.2 Gene Significance (GS) and Module Membership (MM)
geneModuleMembership <- as.data.frame(WGCNA::cor(datExpr32, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr32)))

geneTraitSignificance <- as.data.frame(WGCNA::cor(datExpr32, trait32, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr32)))
names(geneTraitSignificance) <- "GS"
names(GSPvalue) <- "p_GS"

# Intramodular Connectivity (kWithin)
adj_net <- adjacency(datExpr32, power = power_chosen, type = "signed")
intra_conn <- intramodularConnectivity(adj_net, module_colors, scaleByMax = FALSE)

# Build comprehensive gene_MM_GS table
genes_vec <- colnames(datExpr32)

gene_df_all <- data.frame(
  gene = genes_vec,
  module_number = module_labels,
  module_color = module_colors,
  kTotal = round(intra_conn$kTotal, 3),
  kWithin = round(intra_conn$kWithin, 3),
  kOut = round(intra_conn$kOut, 3),
  kDiff = round(intra_conn$kDiff, 3),
  GS = round(geneTraitSignificance$GS, 4),
  p_GS = GSPvalue$p_GS
)

# Extract MM for each gene's assigned module
assigned_me_names <- paste0("ME", gene_df_all$module_color)
gene_df_all$MM <- NA_real_
gene_df_all$p_MM <- NA_real_

for (i in 1:nrow(gene_df_all)) {
  me_col <- assigned_me_names[i]
  if (me_col %in% colnames(geneModuleMembership)) {
    gene_df_all$MM[i] <- round(geneModuleMembership[i, me_col], 4)
    gene_df_all$p_MM[i] <- MMPvalue[i, me_col]
  }
}

# Rank within module by kWithin (decreasing rank)
gene_df_all$rank_in_module <- ave(-gene_df_all$kWithin, gene_df_all$module_color, FUN = function(x) rank(x, ties.method = "min"))

write.csv(gene_df_all, "results/tables/wgcna/gene_MM_GS.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/gene_MM_GS.csv")

# MM vs GS Scatter plots for Key Modules
for (km in key_modules) {
  km_genes <- gene_df_all[gene_df_all$module_color == km, ]
  cor_val <- round(cor(km_genes$MM, km_genes$GS, use = "complete.obs"), 3)
  
  p_mm_gs <- ggplot(km_genes, aes(x = abs(MM), y = abs(GS))) +
    geom_point(color = km, size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("Module Membership vs Gene Significance: ", km, " module"),
      subtitle = paste0("Correlation: r = ", cor_val, " (n = ", nrow(km_genes), " genes)"),
      x = paste0("Module Membership |MM| in ", km),
      y = "Gene Significance |GS| for IPF"
    ) +
    geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "red")
  
  save_plot(p_mm_gs, paste0("results/figures/wgcna/MM_vs_GS_", km), width = 7, height = 6)
}

# ------------------------------------------------------------------------------
# STEP 6: Hub Gene Identification
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 6: Hub Gene Identification ===")

# Hub thresholds (configurable parameters)
CUTOFF_MM <- 0.8
CUTOFF_GS <- 0.2

# Define module hubs in key modules: key module AND |MM| >= 0.8 AND |GS| >= 0.2
hub_candidates <- gene_df_all[
  gene_df_all$module_color %in% key_modules &
  abs(gene_df_all$MM) >= CUTOFF_MM &
  abs(gene_df_all$GS) >= CUTOFF_GS,
]

log_msg("Total module hub candidates (|MM| >= ", CUTOFF_MM, ", |GS| >= ", CUTOFF_GS, " in key modules): ", nrow(hub_candidates))

# 3 Lists:
# 1. hubs_primary: hub genes that are also in the 45 shared matrisome DEGs
hubs_primary <- hub_candidates[hub_candidates$gene %in% shared45_genes, ]

# 2. hubs_matrisome: hub genes that are in the matrisome masterlist (excluding Non-matrisome)
mat_unique_annot <- read.csv("results/tables/GSE32539_matriannotated_unique_gene.csv")
mat_gene_names <- mat_unique_annot$Annotated.Gene
if (is.null(mat_gene_names)) mat_gene_names <- mat_unique_annot$`Annotated Gene`
mat_divs <- mat_unique_annot$Annotated.Matrisome.Division
if (is.null(mat_divs)) mat_divs <- mat_unique_annot$`Annotated Matrisome Division`
mat_cats <- mat_unique_annot$Annotated.Matrisome.Category
if (is.null(mat_cats)) mat_cats <- mat_unique_annot$`Annotated Matrisome Category`

mat_true_genes <- mat_gene_names[mat_divs != "Non-matrisome"]
mat_divisions <- setNames(mat_divs, mat_gene_names)
mat_categories <- setNames(mat_cats, mat_gene_names)

hubs_matrisome <- hub_candidates[hub_candidates$gene %in% mat_true_genes, ]

# 3. hubs_all: all module hub genes in key modules
hubs_all <- hub_candidates

log_msg("Hub Gene Counts:")
log_msg("  1. hubs_primary (in 45 shared matrisome DEGs): ", nrow(hubs_primary), " genes")
log_msg("  2. hubs_matrisome (matrisome hub genes in key modules): ", nrow(hubs_matrisome), " genes")
log_msg("  3. hubs_all (all hub genes in key modules): ", nrow(hubs_all), " genes")

# Save the three lists (CSV & plain TXT)
write.csv(hubs_primary, "results/tables/wgcna/hubs_primary.csv", row.names = FALSE)
writeLines(hubs_primary$gene, "results/tables/wgcna/hubs_primary.txt")

write.csv(hubs_matrisome, "results/tables/wgcna/hubs_matrisome.csv", row.names = FALSE)
writeLines(hubs_matrisome$gene, "results/tables/wgcna/hubs_matrisome.txt")

write.csv(hubs_all, "results/tables/wgcna/hubs_all.csv", row.names = FALSE)
writeLines(hubs_all$gene, "results/tables/wgcna/hubs_all.txt")

# Load limma results from all 3 datasets for final table annotation
limma24 <- read.csv("results/tables/GSE24206_limma_full.csv")
limma32 <- read.csv("results/tables/GSE32539_limma_full.csv")
limma110 <- read.csv("results/tables/GSE110147_limma_full.csv")

rownames(limma24) <- limma24$gene
rownames(limma32) <- limma32$gene
rownames(limma110) <- limma110$gene

# Thesis 15 hub genes reference list
thesis_15 <- c(
  "ASPN", "COL1A1", "COL3A1", "COL5A2", "COL14A1", "COL15A1", "COMP",
  "IL6", "MMP7", "MMP8", "MMP10", "OGN", "SFRP4", "TGFB3", "TIMP3"
)

# Build comprehensive final table: union of all hub candidates
all_hub_genes <- sort(unique(hubs_all$gene))

final_hub_table <- data.frame(
  gene = all_hub_genes,
  matrisome_division = ifelse(all_hub_genes %in% names(mat_divisions), mat_divisions[all_hub_genes], "Non-matrisome"),
  matrisome_category = ifelse(all_hub_genes %in% names(mat_categories), mat_categories[all_hub_genes], "Non-matrisome"),
  module = gene_df_all$module_color[match(all_hub_genes, gene_df_all$gene)],
  MM = gene_df_all$MM[match(all_hub_genes, gene_df_all$gene)],
  GS = gene_df_all$GS[match(all_hub_genes, gene_df_all$gene)],
  kWithin = gene_df_all$kWithin[match(all_hub_genes, gene_df_all$gene)],
  rank_in_module = gene_df_all$rank_in_module[match(all_hub_genes, gene_df_all$gene)],
  in_45_shared_DEGs = ifelse(all_hub_genes %in% shared45_genes, "yes", "no"),
  direction_GSE24206 = ifelse(limma24[all_hub_genes, "logFC"] > 0, "Up", "Down"),
  direction_GSE32539 = ifelse(limma32[all_hub_genes, "logFC"] > 0, "Up", "Down"),
  direction_GSE110147 = ifelse(limma110[all_hub_genes, "logFC"] > 0, "Up", "Down"),
  logFC_GSE24206 = round(limma24[all_hub_genes, "logFC"], 3),
  adj.P_GSE24206 = signif(limma24[all_hub_genes, "adj.P.Val"], 3),
  logFC_GSE32539 = round(limma32[all_hub_genes, "logFC"], 3),
  adj.P_GSE32539 = signif(limma32[all_hub_genes, "adj.P.Val"], 3),
  logFC_GSE110147 = round(limma110[all_hub_genes, "logFC"], 3),
  adj.P_GSE110147 = signif(limma110[all_hub_genes, "adj.P.Val"], 3),
  in_thesis_15_hubs = ifelse(all_hub_genes %in% thesis_15, "yes", "no")
)

# Concordant direction evaluation
final_hub_table$concordant <- apply(
  final_hub_table[, c("direction_GSE24206", "direction_GSE32539", "direction_GSE110147")], 1,
  function(x) {
    if (any(is.na(x))) return("no")
    if (length(unique(x)) == 1) "yes" else "no"
  }
)

# Flag known discordant genes
final_hub_table$concordant[final_hub_table$gene %in% c("CLEC4E", "MMP8", "PLXDC1")] <- "discordant_flag"

# Sort final hub table by in_45_shared_DEGs (yes first), then kWithin descending
final_hub_table <- final_hub_table[order(final_hub_table$in_45_shared_DEGs == "yes", final_hub_table$kWithin, decreasing = TRUE), ]
write.csv(final_hub_table, "results/tables/wgcna/hub_genes_final.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/hub_genes_final.csv (", nrow(final_hub_table), " total hub candidates)")

# Guidance on RF & SVM-RFE inputs if primary list has < 8 genes
if (nrow(hubs_primary) < 8) {
  log_msg("NOTE: hubs_primary has ", nrow(hubs_primary), " genes (< 8).")
  log_msg("RECOMMENDATION: Use 'hubs_matrisome' (n = ", nrow(hubs_matrisome), ") or the combined top-ranked candidates for RF/SVM-RFE feature selection to ensure adequate model stability and avoid underfitting.")
} else {
  log_msg("hubs_primary contains ", nrow(hubs_primary), " genes (>= 8). It is sufficient for ML feature selection, but hubs_matrisome offers a broader biological feature space.")
}

# ------------------------------------------------------------------------------
# STEP 7: Biological Validation & Module Preservation
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 7: Biological Validation & Module Preservation ===")

# 7.1 Fisher's Exact Test: Matrisome Enrichment per Module
mat_enrich_list <- list()
all_net_genes <- colnames(datExpr32)
is_mat_net <- all_net_genes %in% mat_true_genes

for (mod in unique(module_colors)) {
  in_mod <- module_colors == mod
  tbl <- table(factor(in_mod, levels = c(TRUE, FALSE)), factor(is_mat_net, levels = c(TRUE, FALSE)))
  
  fet <- fisher.test(tbl)
  mat_enrich_list[[mod]] <- data.frame(
    module = mod,
    module_size = sum(in_mod),
    matrisome_genes_in_mod = tbl[1, 1],
    odds_ratio = as.numeric(fet$estimate),
    p_value = fet$p.value
  )
}
mat_enrich_df <- do.call(rbind, mat_enrich_list)
mat_enrich_df$adj_p_value <- p.adjust(mat_enrich_df$p_value, method = "BH")
mat_enrich_df <- mat_enrich_df[order(mat_enrich_df$p_value), ]
write.csv(mat_enrich_df, "results/tables/wgcna/matrisome_module_enrichment.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/matrisome_module_enrichment.csv")

# 7.2 Module membership of 45 shared DEGs
deg_mod_df <- data.frame(
  gene = shared45_genes,
  module_number = module_labels[match(shared45_genes, colnames(datExpr32))],
  module_color = module_colors[match(shared45_genes, colnames(datExpr32))],
  MM = gene_df_all$MM[match(shared45_genes, gene_df_all$gene)],
  GS = gene_df_all$GS[match(shared45_genes, gene_df_all$gene)],
  kWithin = gene_df_all$kWithin[match(shared45_genes, gene_df_all$gene)],
  rank_in_module = gene_df_all$rank_in_module[match(shared45_genes, gene_df_all$gene)]
)
deg_mod_df <- deg_mod_df[order(deg_mod_df$module_color, -deg_mod_df$kWithin), ]
write.csv(deg_mod_df, "results/tables/wgcna/DEG_by_module.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/DEG_by_module.csv")

log_msg("Distribution of 45 shared DEGs across modules:")
print(table(deg_mod_df$module_color))

# 7.3 Functional Enrichment (GO BP & KEGG) for Key Modules
log_msg("Running functional enrichment analysis (GO BP & KEGG) for key module(s)...")

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

entrez_universe <- bitr(all_net_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

for (km in key_modules) {
  km_symbols <- gene_df_all$gene[gene_df_all$module_color == km]
  km_entrez_df <- bitr(km_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  km_entrez <- km_entrez_df$ENTREZID
  
  if (length(km_entrez) > 10) {
    ego <- enrichGO(
      gene = km_entrez,
      universe = entrez_universe,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      readable = TRUE
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      ego_res <- as.data.frame(ego)
      write.csv(ego_res, paste0("results/tables/wgcna/enrichment_", km, "_GOBP.csv"), row.names = FALSE)
      
      # Top 10 terms bubble plot
      top10 <- head(ego_res, 10)
      top10$GeneRatio_num <- sapply(strsplit(top10$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      top10$Description <- factor(top10$Description, levels = rev(top10$Description))
      
      p_bubble <- ggplot(top10, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
        geom_point(alpha = 0.85) +
        scale_color_gradient(low = "red", high = "blue") +
        theme_bw(base_size = 12) +
        labs(
          title = paste0("GO Biological Process: ", km, " Module"),
          x = "Gene Ratio",
          y = "",
          color = "Adjusted P",
          size = "Count"
        )
      save_plot(p_bubble, paste0("results/figures/wgcna/enrichment_bubble_", km), width = 9, height = 6)
      log_msg("Generated GO BP enrichment and bubble plot for ", km, " module (", nrow(ego_res), " enriched terms).")
    }
  }
}

# 7.4 Module Preservation in GSE110147
log_msg("Evaluating module preservation of GSE32539 modules in GSE110147 (N=33)...")

expr110 <- read.csv("data/processed/GSE110147_expr_gene_level.csv", check.names = FALSE)
rownames(expr110) <- expr110$gene
meta110 <- read.csv("data/processed/GSE110147_metadata.csv")

# Filter to identical network genes
common_genes <- colnames(datExpr32)
sample_cols110 <- meta110$sample_id

mat110 <- as.matrix(expr110[common_genes, sample_cols110])
class(mat110) <- "numeric"
datExpr110_pres <- t(mat110) # samples x genes

datExpr32_pres <- datExpr32[, common_genes]
colors_pres <- module_colors

multiData <- list(
  GSE32539 = list(data = datExpr32_pres),
  GSE110147 = list(data = datExpr110_pres)
)
multiColor <- list(GSE32539 = colors_pres)

# Run module preservation with 100 permutations as requested
log_msg("Starting modulePreservation calculation (100 permutations, maxModuleSize=1000)...")
pres <- modulePreservation(
  multiData,
  multiColor,
  referenceNetworks = 1,
  nPermutations = 100,
  maxModuleSize = 1000,
  randomSeed = 20261102,
  verbose = 2
)

# Extract preservation statistics
pres_stats <- pres$preservation$Z$ref.GSE32539$inColumnsAlsoPresentIn.GSE110147
pres_summary <- data.frame(
  module = rownames(pres_stats),
  module_size = pres$quality$observed$ref.GSE32539$inColumnsAlsoPresentIn.GSE110147$moduleSize,
  Zsummary = pres_stats$Zsummary.pres,
  medianRank = pres$preservation$observed$ref.GSE32539$inColumnsAlsoPresentIn.GSE110147$medianRank.pres
)
pres_summary$preservation_level <- ifelse(
  pres_summary$Zsummary > 10, "strong",
  ifelse(pres_summary$Zsummary >= 2, "moderate", "none/low")
)
pres_summary <- pres_summary[order(-pres_summary$Zsummary), ]
write.csv(pres_summary, "results/tables/wgcna/module_preservation_GSE110147.csv", row.names = FALSE)
log_msg("Saved results/tables/wgcna/module_preservation_GSE110147.csv")

log_msg("Module Preservation Results (GSE32539 in GSE110147):")
for (i in 1:nrow(pres_summary)) {
  m <- pres_summary$module[i]
  if (m != "gold") {
    log_msg("  Module ", m, ": Zsummary = ", round(pres_summary$Zsummary[i], 2), " (", pres_summary$preservation_level[i], ")")
  }
}
log_msg("Interpretation: GSE110147 has N=33 samples; moderate Zsummary (2-10) reflects limited statistical sample power rather than biological non-reproducibility.")

# Plot preservation Zsummary
p_pres <- ggplot(pres_summary[pres_summary$module != "gold", ], aes(x = module_size, y = Zsummary, label = module)) +
  geom_point(aes(color = preservation_level), size = 4) +
  geom_text(vjust = -0.7, size = 3.5) +
  geom_hline(yintercept = c(2, 10), linetype = "dashed", color = c("blue", "darkgreen")) +
  scale_color_manual(values = c("strong" = "darkgreen", "moderate" = "blue", "none/low" = "red")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Module Preservation in GSE110147",
    subtitle = "Reference: GSE32539 (nPerm = 100)",
    x = "Module Size",
    y = "Zsummary Preservation Statistic",
    color = "Preservation"
  )

save_plot(p_pres, "results/figures/wgcna/module_preservation_Zsummary", width = 8, height = 6)

# 7.5 Network Figure for Key Modules (Cytoscape Export & igraph Visualization)
log_msg("Exporting network graph and Cytoscape files for key module hub genes...")

for (km in key_modules) {
  km_hubs <- hubs_all[hubs_all$module_color == km, ]
  if (nrow(km_hubs) > 1) {
    # Top 30 hub genes ranked by kWithin
    hub_genes_km <- head(km_hubs$gene[order(-km_hubs$kWithin)], 30)
    
    sub_expr_km <- datExpr32[, hub_genes_km]
    adj_hubs <- adjacency(sub_expr_km, power = power_chosen, type = "signed")
    tom_hubs <- TOMsimilarity(adj_hubs)
    dimnames(tom_hubs) <- list(hub_genes_km, hub_genes_km)
    
    cyt <- exportNetworkToCytoscape(
      tom_hubs,
      edgeFile = paste0("results/tables/wgcna/network_edges_", km, ".txt"),
      nodeFile = paste0("results/tables/wgcna/network_nodes_", km, ".txt"),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = hub_genes_km,
      altNodeNames = hub_genes_km,
      nodeAttr = ifelse(hub_genes_km %in% shared45_genes, "Shared_DEG", "Module_Hub")
    )
    
    # igraph network diagram
    g_net <- graph_from_adjacency_matrix(adj_hubs, mode = "undirected", weighted = TRUE, diag = FALSE)
    q75 <- quantile(E(g_net)$weight, 0.75)
    g_sub <- delete_edges(g_net, E(g_net)[E(g_net)$weight < q75])
    
    # Node attributes
    V(g_sub)$is_shared <- V(g_sub)$name %in% shared45_genes
    V(g_sub)$color <- ifelse(V(g_sub)$is_shared, "#d7191c", "#2b83ba")
    
    png(paste0("results/figures/wgcna/network_graph_", km, ".png"), width = 8, height = 8, units = "in", res = 300)
    set.seed(20261102)
    plot(
      g_sub,
      vertex.label = V(g_sub)$name,
      vertex.label.cex = 0.8,
      vertex.label.color = "black",
      vertex.size = ifelse(V(g_sub)$is_shared, 12, 8),
      vertex.color = V(g_sub)$color,
      edge.color = "grey75",
      edge.width = E(g_sub)$weight * 3,
      layout = layout_with_fr(g_sub),
      main = paste0("Hub Gene Co-expression Network: ", km, " Module\n(Red = 45 Shared Matrisome DEG, Blue = Module Hub)")
    )
    dev.off()
    
    pdf(paste0("results/figures/wgcna/network_graph_", km, ".pdf"), width = 8, height = 8)
    set.seed(20261102)
    plot(
      g_sub,
      vertex.label = V(g_sub)$name,
      vertex.label.cex = 0.8,
      vertex.label.color = "black",
      vertex.size = ifelse(V(g_sub)$is_shared, 12, 8),
      vertex.color = V(g_sub)$color,
      edge.color = "grey75",
      edge.width = E(g_sub)$weight * 3,
      layout = layout_with_fr(g_sub),
      main = paste0("Hub Gene Co-expression Network: ", km, " Module\n(Red = 45 Shared Matrisome DEG, Blue = Module Hub)")
    )
    dev.off()
    log_msg("Saved network plots and Cytoscape edge/node files for ", km, " module.")
  }
}

# ------------------------------------------------------------------------------
# STEP 8: Hand-off Preparation for Downstream Machine Learning (RF / SVM-RFE)
# ------------------------------------------------------------------------------
log_msg("\n=== STEP 8: ML Hand-off Preparation ===")

# Candidate gene sets definition
union_candidates <- unique(c(hubs_primary$gene, hubs_matrisome$gene, hubs_all$gene))

cand_sets_df <- data.frame(gene = union_candidates)
cand_sets_df$in_primary <- ifelse(cand_sets_df$gene %in% hubs_primary$gene, "yes", "no")
cand_sets_df$in_matrisome <- ifelse(cand_sets_df$gene %in% hubs_matrisome$gene, "yes", "no")
cand_sets_df$in_all <- ifelse(cand_sets_df$gene %in% hubs_all$gene, "yes", "no")

cand_sets_df$set_membership <- apply(
  cand_sets_df[, c("in_primary", "in_matrisome", "in_all")], 1,
  function(x) {
    sets <- c()
    if (x[1] == "yes") sets <- c(sets, "primary")
    if (x[2] == "yes") sets <- c(sets, "matrisome")
    if (x[3] == "yes") sets <- c(sets, "all")
    paste(sets, collapse = ";")
  }
)
write.csv(cand_sets_df, "data/processed/ml_input/candidate_gene_sets.csv", row.names = FALSE)
log_msg("Saved data/processed/ml_input/candidate_gene_sets.csv (", nrow(cand_sets_df), " candidate genes)")

# Export samples x genes log2 matrices and labels for all 3 datasets
gse_list <- c("GSE32539", "GSE110147", "GSE24206")

for (g in gse_list) {
  g_expr_file <- sprintf("data/processed/%s_expr_gene_level.csv", g)
  g_meta_file <- sprintf("data/processed/%s_metadata.csv", g)
  
  g_expr <- read.csv(g_expr_file, check.names = FALSE)
  rownames(g_expr) <- g_expr$gene
  g_meta <- read.csv(g_meta_file)
  
  # Check which candidate genes are present
  avail_genes <- intersect(union_candidates, g_expr$gene)
  missing_genes <- setdiff(union_candidates, g_expr$gene)
  
  log_msg(g, ": ", length(avail_genes), " of ", length(union_candidates), " candidate genes present.")
  if (length(missing_genes) > 0) {
    log_msg("  Missing genes in ", g, ": ", paste(missing_genes, collapse = ", "))
  }
  
  # Matrix: samples in rows, candidate genes in columns
  g_sub <- as.matrix(g_expr[avail_genes, g_meta$sample_id])
  g_samples_x_genes <- t(g_sub)
  
  write.csv(g_samples_x_genes, sprintf("data/processed/ml_input/%s_expr_hubs_all_candidates.csv", g), row.names = TRUE)
  
  # Labels file: sample_id, label
  g_labels <- data.frame(
    sample_id = g_meta$sample_id,
    label = g_meta$group
  )
  write.csv(g_labels, sprintf("data/processed/ml_input/%s_labels.csv", g), row.names = FALSE)
  log_msg("Saved ML inputs for ", g, " (", nrow(g_samples_x_genes), " samples x ", ncol(g_samples_x_genes), " genes)")
}

log_msg("\n=== WGCNA Pipeline Complete Successfully ===")
