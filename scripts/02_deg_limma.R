suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
})

set.seed(42)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
log_dir <- file.path(project_dir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_path <- file.path(log_dir, "pipeline.log")
if (file.exists(log_path)) file.remove(log_path)

log_msg <- function(msg) {
  cat(msg, "\n", sep = "")
  cat(sprintf("%s\n", msg), file = log_path, append = TRUE)
}

if (!dir.exists(file.path(project_dir, "results", "tables"))) {
  dir.create(file.path(project_dir, "results", "tables"), recursive = TRUE)
}

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

get_platform_annotation <- function(gpl_accession) {
  raw_dir <- file.path(project_dir, "data", "raw")
  if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

  gpl <- tryCatch({
    GEOquery::getGEO(GEO = gpl_accession, destdir = raw_dir, getGPL = TRUE, GSEMatrix = FALSE)
  }, error = function(e) NULL)

  if (is.null(gpl)) return(NULL)

  if (inherits(gpl, "GPL")) {
    table_df <- gpl@dataTable@table
  } else if (is.list(gpl) && length(gpl) > 0 && inherits(gpl[[1]], "GPL")) {
    table_df <- gpl[[1]]@dataTable@table
  } else {
    return(NULL)
  }

  if (is.null(table_df) || nrow(table_df) == 0) return(NULL)
  table_df
}

map_gene_symbols <- function(probes, gpl_accession) {
  gpl_table <- get_platform_annotation(gpl_accession)
  if (is.null(gpl_table)) {
    log_msg(sprintf("[%s] platform annotation unavailable; probe symbols cannot be mapped.", gpl_accession))
    return(setNames(rep(NA_character_, length(probes)), probes))
  }

  key_col <- names(gpl_table)[toupper(names(gpl_table)) %in% c("ID", "PROBE_ID", "PROBESET_ID", "TRANSCRIPT_CLUSTER_ID")]
  if (length(key_col) == 0) key_col <- "ID"

  symbol_col <- names(gpl_table)[grepl("symbol|gene assignment|gene_assignment|mrna_assignment", names(gpl_table), ignore.case = TRUE)]
  if (length(symbol_col) == 0) {
    symbol_col <- names(gpl_table)[grepl("gene", names(gpl_table), ignore.case = TRUE)]
  }
  if (length(symbol_col) == 0) {
    log_msg(sprintf("[%s] no gene_symbol-like column found in GPL annotation.", gpl_accession))
    return(setNames(rep(NA_character_, length(probes)), probes))
  }

  symbol_col <- symbol_col[1]
  gpl_table[[symbol_col]] <- as.character(gpl_table[[symbol_col]])
  symbol_vec <- vapply(gpl_table[[symbol_col]], extract_probe_symbol, character(1))
  names(symbol_vec) <- as.character(gpl_table[[key_col]])
  symbol_vec[is.na(symbol_vec)] <- NA_character_
  symbol_vec[probes]
}

for (gse in c("GSE24206", "GSE32539", "GSE110147")) {
  log_msg(sprintf("[%s] running limma DEG analysis", gse))
  metadata <- read.csv(file.path(project_dir, "data", "processed", sprintf("%s_metadata.csv", gse)), stringsAsFactors = FALSE)
  metadata$group <- factor(metadata$group, levels = c("Control", "IPF"))

  expr_path <- file.path(project_dir, "data", "processed", sprintf("%s_expression.csv", gse))
  expr_df <- read.csv(expr_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"probe_id" %in% names(expr_df)) stop(sprintf("[%s] expression file is missing a probe_id column.", gse))

  expr_mat <- as.matrix(expr_df[, setdiff(names(expr_df), "probe_id"), drop = FALSE])
  rownames(expr_mat) <- expr_df$probe_id
  expr_mat <- expr_mat[, metadata$sample_id, drop = FALSE]
  expr_mat <- expr_mat[, match(metadata$sample_id, colnames(expr_mat)), drop = FALSE]

  if (ncol(expr_mat) != nrow(metadata)) stop(sprintf("Sample count mismatch in %s", gse))

  design <- model.matrix(~ group, data = metadata)
  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = "groupIPF", number = Inf, sort.by = "none")
  res$probe <- rownames(res)

  gpl_accession <- if (gse == "GSE24206") "GPL570" else "GPL6244"
  gene_map <- map_gene_symbols(res$probe, gpl_accession)
  res$gene <- gene_map[res$probe]
  res <- res[!is.na(res$gene) & res$gene != "", , drop = FALSE]

  if (nrow(res) > 0) {
    # Keep the probe with the largest absolute logFC; if tied, keep the one with higher mean expression.
    res <- res[order(abs(res$logFC), -res$AveExpr, decreasing = TRUE), , drop = FALSE]
    res <- res[!duplicated(res$gene), , drop = FALSE]
  }

  res <- res[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  names(res) <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")

  write.csv(res, file = file.path(project_dir, "results", "tables", sprintf("%s_limma_full.csv", gse)), row.names = FALSE)
  log_msg(sprintf("[%s] saved %d genes with limma results", gse, nrow(res)))
}

log_msg("[STEP 2] Differential expression complete.")
