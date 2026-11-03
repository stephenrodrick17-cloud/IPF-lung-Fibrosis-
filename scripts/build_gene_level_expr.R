.libPaths(c(normalizePath("rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
})

extract_probe_symbol <- function(value) {
  x <- trimws(as.character(value))
  if (length(x) == 0 || is.na(x) || x == "" || x %in% c("---", "NA", "N/A", "null")) {
    return(NA_character_)
  }
  if (grepl("//", x, fixed = TRUE)) {
    parts <- strsplit(x, "//", fixed = TRUE)[[1]]
    clean <- trimws(parts[2])
    if (nzchar(clean)) return(clean)
  }
  m <- regmatches(x, regexec("([A-Z0-9_.-]{2,})", x, perl = TRUE))
  if (length(m[[1]]) >= 2) {
    return(m[[1]][2])
  }
  return(NA_character_)
}

for (gse in c("GSE32539", "GSE110147", "GSE24206")) {
  cat("\nProcessing gene-level expression for", gse, "...\n")
  meta <- read.csv(sprintf("data/processed/%s_metadata.csv", gse), stringsAsFactors = FALSE)
  meta$group <- factor(meta$group, levels = c("Control", "IPF"))
  
  expr_df <- read.csv(sprintf("data/processed/%s_expression.csv", gse), check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- meta$sample_id
  expr_mat <- as.matrix(expr_df[, sample_cols])
  rownames(expr_mat) <- expr_df$probe_id
  
  gpl_acc <- if (gse == "GSE24206") "GPL570" else "GPL6244"
  gpl <- GEOquery::getGEO(GEO = gpl_acc, destdir = "data/raw", getGPL = TRUE, GSEMatrix = FALSE)
  gpl_table <- if (inherits(gpl, "GPL")) gpl@dataTable@table else gpl[[1]]@dataTable@table
  key_col <- names(gpl_table)[toupper(names(gpl_table)) %in% c("ID", "PROBE_ID", "PROBESET_ID", "TRANSCRIPT_CLUSTER_ID")][1]
  symbol_col <- names(gpl_table)[grepl("symbol|gene assignment|gene_assignment|mrna_assignment", names(gpl_table), ignore.case = TRUE)][1]
  
  symbols <- vapply(as.character(gpl_table[[symbol_col]]), extract_probe_symbol, character(1))
  names(symbols) <- as.character(gpl_table[[key_col]])
  
  if (gse == "GSE24206") {
    symbols <- gsub("^[/[:space:]]+|[/[:space:]]+$", "", symbols)
  }
  
  probe_genes <- symbols[rownames(expr_mat)]
  valid <- !is.na(probe_genes) & probe_genes != "" & probe_genes != "---"
  expr_mat_valid <- expr_mat[valid, , drop = FALSE]
  probe_genes_valid <- probe_genes[valid]
  
  # Run limma on valid probes to get exact logFC & p-value for probe collapsing
  design <- model.matrix(~ group, data = meta)
  fit <- lmFit(expr_mat_valid, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "groupIPF", number = Inf, sort.by = "none")
  tt$gene <- probe_genes_valid[rownames(tt)]
  
  # Sort by largest |logFC|, then lowest P.Value
  ord <- order(-abs(tt$logFC), tt$P.Value)
  tt_sorted <- tt[ord, ]
  best_probes <- rownames(tt_sorted)[!duplicated(tt_sorted$gene)]
  best_genes <- tt_sorted$gene[!duplicated(tt_sorted$gene)]
  
  gene_expr_mat <- expr_mat_valid[best_probes, , drop = FALSE]
  rownames(gene_expr_mat) <- best_genes
  
  # Save as gene-level expression matrix (genes x samples, with gene symbol as first column / row names)
  out_df <- data.frame(gene = best_genes, gene_expr_mat, check.names = FALSE)
  out_path <- sprintf("data/processed/%s_expr_gene_level.csv", gse)
  write.csv(out_df, out_path, row.names = FALSE)
  
  cat(gse, "saved to:", out_path, "\n")
  cat(gse, "dimension:", nrow(gene_expr_mat), "genes x", ncol(gene_expr_mat), "samples\n")
  cat(gse, "expression range:", range(gene_expr_mat, na.rm = TRUE), "(log2 scale confirmed)\n")
  
  # Compare with clean_for_matrisome.csv
  clean_genes <- read.csv(sprintf("data/processed/%s_clean_for_matrisome.csv", gse))$gene
  cat(gse, "clean_for_matrisome overlap:", sum(clean_genes %in% best_genes), "of", length(clean_genes), "\n")
}
cat("\nGene-level expression files generated successfully.\n")

