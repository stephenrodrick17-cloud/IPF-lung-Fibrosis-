suppressPackageStartupMessages({
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
})

set.seed(42)

log_msg <- function(msg) {
  cat(msg, "\n", sep = "")
  write(msg, file = "logs/pipeline.log", append = TRUE)
}

project_dir <- getwd()

for (gse in c("GSE24206", "GSE32539", "GSE110147")) {
  deg <- read.csv(file.path(project_dir, "results", "tables", sprintf("%s_limma_full.csv", gse)), stringsAsFactors = FALSE)
  mat <- read.csv(file.path(project_dir, "results", "tables", sprintf("%s_matrisome_DEGs.csv", gse)), stringsAsFactors = FALSE)
  if (nrow(mat) == 0) {
    log_msg(sprintf("[%s] no matrisome genes detected", gse))
    next
  }

  # Two threshold modes: A and B.
  modeA <- subset(mat, adj.P.Val < 0.05 & abs(logFC) > 0.585)
  modeB <- subset(mat, P.Value < 0.05 & abs(logFC) > 0.585)

  modeA$direction <- ifelse(modeA$logFC > 0, "up", "down")
  modeB$direction <- ifelse(modeB$logFC > 0, "up", "down")

  write.csv(modeA, file = file.path(project_dir, "results", "tables", sprintf("%s_matrisome_DEGs_modeA.csv", gse)), row.names = FALSE)
  write.csv(modeB, file = file.path(project_dir, "results", "tables", sprintf("%s_matrisome_DEGs_modeB.csv", gse)), row.names = FALSE)

  # Volcano plot using mode A as primary.
  p <- ggplot(deg, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.585, ifelse(logFC > 0, "Up", "Down"), "Not significant")), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Up" = "tomato", "Down" = "cornflowerblue", "Not significant" = "grey70")) +
    labs(color = "Status", x = "log2 fold change", y = "-log10(P)", title = sprintf("%s volcano plot", gse)) +
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed")
  ggsave(file.path(project_dir, "results", "figures", sprintf("%s_volcano.png", gse)), p, dpi = 300, width = 8, height = 6)
  ggsave(file.path(project_dir, "results", "figures", sprintf("%s_volcano.pdf", gse)), p, width = 8, height = 6)

  # Heatmap of matrisome DEGs.
  if (nrow(modeA) > 1) {
    heat <- reshape2::acast(modeA[, c("gene", "logFC")], gene ~ ., value.var = "logFC")
    if (is.null(dim(heat))) {
      next
    }
    pheat <- pheatmap(heat, cluster_rows = TRUE, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), main = sprintf("%s matrisome DEG heatmap", gse), angle_col = 45)
    save_pheatmap_pdf(pheat, file = file.path(project_dir, "results", "figures", sprintf("%s_matrisome_heatmap.pdf", gse)))
  }

  log_msg(sprintf("[%s] mode A count = %d; mode B count = %d; saved volcano and heatmap outputs", gse, nrow(modeA), nrow(modeB)))
}

log_msg("[STEP 4] DEG thresholding and visualization complete.")
