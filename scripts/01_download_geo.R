suppressPackageStartupMessages({
  library(GEOquery)
  library(ggplot2)
  library(limma)
})

set.seed(42)

log_msg <- function(msg) {
  cat(msg, "\n", sep = "")
  write(msg, file = "logs/pipeline.log", append = TRUE)
}

# Run from the project root, so relative paths work.
project_dir <- getwd()
if (!dir.exists(file.path(project_dir, "logs"))) dir.create(file.path(project_dir, "logs"), recursive = TRUE)
if (!dir.exists(file.path(project_dir, "data", "raw"))) dir.create(file.path(project_dir, "data", "raw"), recursive = TRUE)
if (!dir.exists(file.path(project_dir, "data", "processed"))) dir.create(file.path(project_dir, "data", "processed"), recursive = TRUE)
if (!dir.exists(file.path(project_dir, "results", "figures"))) dir.create(file.path(project_dir, "results", "figures"), recursive = TRUE)

prefixes <- c("GSE24206", "GSE32539", "GSE110147", "GSE135065")

make_group <- function(x) {
  s <- tolower(x)
  if (grepl("ipf|idiopathic pulmonary fibrosis", s)) return("IPF")
  if (grepl("control|normal|healthy|nonipf|non-ipf|non ipf|non-disease", s)) return("Control")
  return("Other")
}

save_boxplot <- function(expr_matrix, gse, title = NULL) {
  # Check whether data are already log2 or raw. If raw, apply log2 transform.
  if (max(expr_matrix, na.rm = TRUE) > 30) {
    expr_matrix <- log2(expr_matrix + 1)
    log_msg(sprintf("[%s] expression matrix was not log2 transformed; applied log2(x + 1)", gse))
  } else {
    log_msg(sprintf("[%s] expression matrix already appears log2 transformed", gse))
  }
  df <- as.data.frame(expr_matrix)
  df_long <- reshape2::melt(as.matrix(df))
  colnames(df_long) <- c("probe", "sample", "expression")
  p <- ggplot(df_long, aes(x = sample, y = expression)) +
    geom_boxplot(width = 0.5, outlier.alpha = 0.2) +
    theme_bw() +
    labs(title = if (is.null(title)) gse else title, x = "Sample", y = "Expression (log2)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  ggsave(file.path(project_dir, "results", "figures", sprintf("%s_sample_boxplot.png", gse)), p, dpi = 300, width = 12, height = 6)
  ggsave(file.path(project_dir, "results", "figures", sprintf("%s_sample_boxplot.pdf", gse)), p, width = 12, height = 6)
}

for (gse in prefixes) {
  log_msg(sprintf("[%s] downloading GEO dataset", gse))
  gset <- getGEO(GEO = gse, GSEMatrix = TRUE, destdir = file.path(project_dir, "data", "raw"), getGPL = TRUE)
  gset <- gset[[1]]

  pdata <- pData(gset)
  if (nrow(pdata) == 0) stop(sprintf("No phenotype data found for %s", gse))

  sample_id <- rownames(pdata)
  if ("geo_accession" %in% colnames(pdata)) sample_id <- pdata$geo_accession
  if ("title" %in% colnames(pdata)) sample_id <- pdata$title

  pdata$sample_id <- sample_id
  pdata$group_raw <- apply(pdata, 1, function(r) {
    vals <- r[grep("characteristics|source_name|title|description|disease|group", names(r), ignore.case = TRUE)]
    paste(vals[!is.na(vals)], collapse = "; ")
  })
  pdata$group <- vapply(pdata$group_raw, make_group, character(1))

  # Print platform and counts
  platform <- unique(pdata$platform_id)
  if (length(platform) == 0) platform <- "unknown"
  cat(sprintf("[%s] Platform(s): %s\n", gse, paste(platform, collapse = "; ")))
  print(table(pdata$group))

  # Apply filtering rules
  if (gse == "GSE24206") {
    filter_df <- pdata[pdata$group %in% c("IPF", "Control"), , drop = FALSE]
    log_msg(sprintf("[%s] GSE24206 kept IPF and Control only; early/advanced IPF labels were collapsed into IPF as they are not separate treatment groups in this workflow.", gse))
  } else if (gse == "GSE32539") {
    filter_df <- pdata[pdata$group %in% c("IPF", "Control"), , drop = FALSE]
    log_msg(sprintf("[%s] GSE32539 kept IPF and Control only; other disease classes were dropped.", gse))
  } else if (gse == "GSE110147") {
    keep_idx <- !grepl("nsip|mixed ipf-nsip|mixed.*nsip|ipf-nsip", pdata$group_raw, ignore.case = TRUE)
    filter_df <- pdata[keep_idx & pdata$group %in% c("IPF", "Control"), , drop = FALSE]
    log_msg(sprintf("[%s] GSE110147 removed NSIP and mixed IPF-NSIP samples, retaining only IPF vs Control.", gse))
  } else if (gse == "GSE135065") {
    filter_df <- pdata[pdata$group %in% c("IPF", "Control"), , drop = FALSE]
    log_msg(sprintf("[%s] validation dataset retained IPF and Control only.", gse))
  } else {
    filter_df <- pdata[pdata$group %in% c("IPF", "Control"), , drop = FALSE]
  }

  if (nrow(filter_df) == 0) stop(sprintf("No retained samples after filtering for %s", gse))
  keep_samples <- filter_df$sample_id
  filter_df <- filter_df[, c("sample_id", "group")]
  colnames(filter_df) <- c("sample_id", "group")

  write.csv(filter_df, file = file.path(project_dir, "data", "processed", sprintf("%s_metadata.csv", gse)), row.names = FALSE)
  log_msg(sprintf("[%s] saved metadata with %d samples (%d IPF, %d Control)", gse, nrow(filter_df), sum(filter_df$group == "IPF"), sum(filter_df$group == "Control")))
  print(table(filter_df$group))

  expr <- exprs(gset)
  expr <- expr[, colnames(expr) %in% keep_samples, drop = FALSE]
  expr <- expr[, match(keep_samples, colnames(expr)), drop = FALSE]
  save_boxplot(expr, gse)
}

log_msg("[STEP 1] GEO metadata retrieval complete.")
