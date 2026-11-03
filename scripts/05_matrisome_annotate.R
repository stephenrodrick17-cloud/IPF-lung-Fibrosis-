#!/usr/bin/env Rscript
## 05_matrisome_annotate.R
## Matrisome DEG identification with MatrisomeAnalyzeR (IPF vs Control; GSE24206, GSE32539, GSE110147)
## Run from the project root:  Rscript scripts/05_matrisome_annotate.R
## Everything is computed from results/tables/<GSE>_limma_full.csv and the package's own matrisome list.

set.seed(20261102)
stopifnot(dir.exists("scripts"), dir.exists("results"))
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))
for (d in c("logs","data/processed","results/tables","results/figures")) dir.create(d, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(MatrisomeAnalyzeR); library(ggplot2) })

## ---------------------------------------------------------------- logging
LOG <- "logs/matrisome_step.log"; if (file.exists(LOG)) file.remove(LOG)
lg <- function(...) { s <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste0(..., collapse = "")); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
WARN <- character(); warnf <- function(...) { s <- paste0(..., collapse = ""); WARN <<- c(WARN, s); lg("WARNING: ", s) }

## ---------------------------------------------------------------- parameters
GSES   <- c("GSE24206", "GSE32539", "GSE110147")
LFC    <- 0.585; PCUT <- 0.05
THESIS <- c(GSE24206 = 358, GSE32539 = 149, GSE110147 = 339)
DIVS   <- c("Core matrisome", "Matrisome-associated")
CATS   <- c("ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors")

## ---------------------------------------------------------------- STEP A: package + provenance
pv  <- as.character(packageVersion("MatrisomeAnalyzeR"))
pd  <- packageDescription("MatrisomeAnalyzeR")
sha <- if (!is.null(pd$GithubSHA1)) pd$GithubSHA1 else NA
lg("STEP A: MatrisomeAnalyzeR version ", pv, "; GitHub ", pd$GithubUsername, "/", pd$GithubRepo, " @ ", sha)
lg("exported: ", paste(ls("package:MatrisomeAnalyzeR"), collapse = ", "))
lg("R: ", R.version.string, "; seed = 20261102")
sink("logs/sessionInfo.txt"); cat("MatrisomeAnalyzeR version:", pv, "\nGitHub:", pd$GithubUsername, "/", pd$GithubRepo, "@", sha, "\nSeed: 20261102\n\n"); print(sessionInfo()); sink()

## built-in test with the bundled testfile (column is "Gene Symbol", verified by str(testfile))
tt <- matriannotate(data = testfile, gene.column = "Gene Symbol", species = "human")
lg("built-in test: testfile rows in = ", nrow(testfile), "; annotated rows out = ", nrow(tt),
   "; matrisome = ", sum(tt$`Annotated Matrisome Division` != "Non-matrisome"))
stopifnot(nrow(tt) == nrow(testfile), all(c("Annotated Gene","Annotated Matrisome Division","Annotated Matrisome Category") %in% colnames(tt)))
lg("built-in test OK. Output columns (first 3 are new): ", paste(head(colnames(tt), 3), collapse = " | "))
WF <- attr(tt, "workflow"); lg("annotation marker attr(,'workflow') = ", WF)
## helper: plotting functions in the package require this attribute; row/column subsetting may drop it
keepattr <- function(x) { attr(x, "workflow") <- WF; x }

## ---------------------------------------------------------------- INPUT: list, head, colnames, map
lg("INPUT: files in results/tables:"); print(list.files("results/tables"))
need <- c("gene","logFC","AveExpr","P.Value","adj.P.Val")
cand <- list(gene = c("gene","Gene","gene_symbol","Gene.symbol","Gene.Symbol","symbol","SYMBOL","Symbol","GeneSymbol"),
             logFC = c("logFC","log2FC","log2FoldChange"), AveExpr = c("AveExpr","aveexpr","baseMean"),
             P.Value = c("P.Value","pvalue","PValue","p_value"), adj.P.Val = c("adj.P.Val","padj","FDR","adj.P.value"))
read_de <- function(g) {
  f <- file.path("results/tables", paste0(g, "_limma_full.csv")); stopifnot(file.exists(f))
  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  cat("\n=====", g, "=====\ncolnames:", paste(colnames(x), collapse = ", "), "\n"); print(head(x, 4))
  mp <- c()
  for (k in need) { hit <- intersect(cand[[k]], colnames(x)); if (!length(hit)) stop("cannot map column for ", k, " in ", g); mp[k] <- hit[1] }
  ch <- mp[mp != names(mp)]
  lg(g, ": column mapping = ", paste(names(mp), "<-", mp, collapse = "; "), if (length(ch)) paste0("  [RENAMED: ", paste(ch, "->", names(ch), collapse = ", "), "]") else "  [no renaming needed]")
  y <- x[, mp]; colnames(y) <- names(mp)
  y$logFC <- as.numeric(y$logFC); y$AveExpr <- as.numeric(y$AveExpr); y$P.Value <- as.numeric(y$P.Value); y$adj.P.Val <- as.numeric(y$adj.P.Val)
  ## FULL-table check
  fr <- mean(y$adj.P.Val >= PCUT, na.rm = TRUE)
  lg(g, ": rows = ", nrow(y), "; fraction adj.P.Val >= 0.05 = ", round(fr, 3), "; max adj.P.Val = ", signif(max(y$adj.P.Val, na.rm = TRUE), 3))
  if (fr < 0.05) warnf(g, ": table looks SIGNIFICANT-ONLY (", round(100 * fr, 1), "% of rows non-significant). Percent-matrisome would be biased.")
  y
}

## ---------------------------------------------------------------- STEP B: clean
clean_symbols <- function(y, g) {
  st <- list(rows_in = nrow(y))
  y$gene <- trimws(as.character(y$gene))
  ## upstream artefact (02_deg_limma.R splits "A /// B" on "//"): symbols such as "/ KRT17" -> strip stray slashes/spaces at the ends
  slash <- !is.na(y$gene) & grepl("^[/[:space:]]+|[/[:space:]]+$", y$gene)
  st$slash_stripped <- sum(slash); y$gene <- gsub("^[/[:space:]]+|[/[:space:]]+$", "", y$gene)
  lg(g, " B0 stray leading/trailing '/' (upstream '///' split artefact) stripped in ", st$slash_stripped, " rows, e.g. ", paste(head(y$gene[slash], 3), collapse = ","))
  bad <- is.na(y$gene) | y$gene == "" | y$gene == "---"
  st$removed_empty <- sum(bad); y <- y[!bad, ]
  lg(g, " B1 empty/NA/'---' symbols: removed ", st$removed_empty, " (", nrow(y), " left)")
  nna <- !is.finite(y$logFC) | !is.finite(y$P.Value) | !is.finite(y$adj.P.Val)
  st$removed_na_stats <- sum(nna); y <- y[!nna, ]
  lg(g, " B1b rows with non-finite logFC/P/adjP: removed ", st$removed_na_stats)
  multi <- grepl("///", y$gene, fixed = TRUE); st$multi_rows <- sum(multi)
  sp <- strsplit(y$gene, "\\s*///\\s*")
  y <- y[rep(seq_len(nrow(y)), lengths(sp)), ]; y$gene <- trimws(unlist(sp))
  y <- y[!(y$gene == "" | y$gene == "---" | is.na(y$gene)), ]
  st$rows_after_split <- nrow(y)
  lg(g, " B2 multi-gene 'A///B' rule = SPLIT into separate rows; rows affected = ", st$multi_rows, "; rows after split = ", nrow(y))
  st$other_sep <- sum(grepl("[;,| ]", y$gene)); lg(g, " B2b symbols containing ; , | or space (left untouched): ", st$other_sep)
  st$lowercase <- sum(y$gene != toupper(y$gene)); lg(g, " B2c symbols not fully upper-case (matching is case-sensitive; left untouched): ", st$lowercase)
  y <- y[order(-abs(y$logFC), y$P.Value), ]           # largest |logFC|, then lowest P.Value
  nd <- sum(duplicated(y$gene)); y <- y[!duplicated(y$gene), ]
  st$dup_removed <- nd; st$rows_out <- nrow(y)
  lg(g, " B3 duplicate collapse (keep max |logFC|, tie -> min P.Value): removed ", nd, "; unique genes = ", nrow(y))
  lg(g, " STEP B total: rows in ", st$rows_in, " -> out ", st$rows_out, " (lost ", st$rows_in - st$rows_out, ")")
  rownames(y) <- NULL
  write.csv(y, file.path("data/processed", paste0(g, "_clean_for_matrisome.csv")), row.names = FALSE)
  list(df = y, st = st)
}

## ---------------------------------------------------------------- helpers: DEG calling, plots
call_degs <- function(a, mode) {
  p <- if (mode == "A") a$adj.P.Val else a$P.Value
  d <- a[!is.na(p) & p < PCUT & abs(a$logFC) > LFC, ]
  d$Direction <- ifelse(d$logFC > LFC, "Up", "Down")
  d
}
fmt_deg <- function(d) data.frame(gene = d$`Annotated Gene`, logFC = d$logFC, P.Value = d$P.Value, adj.P.Val = d$adj.P.Val,
                                  Direction = d$Direction, Division = d$`Annotated Matrisome Division`,
                                  Category = d$`Annotated Matrisome Category`, Other_annotations = d$Other_annotations, stringsAsFactors = FALSE)
save_plot <- function(p, base, w = 9, h = 7) {
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
}

## ---------------------------------------------------------------- main loop
CLEAN <- list(); ANN <- list(); MAT <- list(); DEG <- list(); SUMM <- list(); CLSTAT <- list(); ANSTAT <- list(); FIGLOG <- list()
for (g in GSES) {
  lg("================ ", g, " ================")
  y  <- read_de(g)
  cl <- clean_symbols(y, g); CLSTAT[[g]] <- cl$st; CLEAN[[g]] <- cl$df

  ## ---- STEP C
  a <- matriannotate(data = cl$df, gene.column = "gene", species = "human")
  lg(g, " C: matriannotate rows in ", nrow(cl$df), " -> rows out ", nrow(a), "; output columns: ", paste(colnames(a), collapse = " | "))
  if (!setequal(a$`Annotated Gene`, cl$df$gene)) warnf(g, ": annotated gene set differs from input genes")
  lg(g, " C: Division levels: ", paste(names(table(a$`Annotated Matrisome Division`)), collapse = ", "))
  lg(g, " C: Category levels: ", paste(names(table(a$`Annotated Matrisome Category`)), collapse = ", "))
  write.csv(a, sprintf("results/tables/%s_matriannotated_all.csv", g), row.names = FALSE)   # raw package output, unaltered
  ## The package returns one row per (gene, category): genes listed in 2 categories (e.g. CLC, NOV, MBP) appear twice.
  ## Keep ONE row per gene so counts are gene counts: Core matrisome takes priority over Matrisome-associated, otherwise first row returned.
  dupg <- unique(a$`Annotated Gene`[duplicated(a$`Annotated Gene`)])
  lg(g, " C: package returned ", nrow(a), " rows for ", nrow(cl$df), " input genes; ", length(dupg), " genes carry >1 annotation: ", paste(dupg, collapse = ","))
  a$.o <- seq_len(nrow(a)); a$.p <- ifelse(a$`Annotated Matrisome Division` == "Core matrisome", 0, 1)
  a <- a[order(a$`Annotated Gene`, a$.p, a$.o), ]
  oth <- tapply(paste0(a$`Annotated Matrisome Division`, ":", a$`Annotated Matrisome Category`), a$`Annotated Gene`, function(v) if (length(v) > 1) paste(v[-1], collapse = " | ") else "")
  a <- a[!duplicated(a$`Annotated Gene`), ]; a$Other_annotations <- unname(oth[a$`Annotated Gene`]); a$.o <- NULL; a$.p <- NULL
  rownames(a) <- NULL; a <- keepattr(a)
  lg(g, " C: after one-row-per-gene collapse: ", nrow(a), " genes (genes with >1 annotation collapsed = ", length(dupg), ")")
  write.csv(a, sprintf("results/tables/%s_matriannotated_unique_gene.csv", g), row.names = FALSE)
  m <- keepattr(a[a$`Annotated Matrisome Division` != "Non-matrisome", ])
  write.csv(m, sprintf("results/tables/%s_matrisome_all_genes.csv", g), row.names = FALSE)
  tdiv <- table(factor(m$`Annotated Matrisome Division`, levels = DIVS)); tcat <- table(factor(m$`Annotated Matrisome Category`, levels = CATS))
  other <- setdiff(unique(m$`Annotated Matrisome Category`), CATS); if (length(other)) warnf(g, ": unexpected categories ", paste(other, collapse = ","))
  lg(g, " C: total genes input = ", nrow(a), "; matrisome genes = ", nrow(m), " (", round(100 * nrow(m) / nrow(a), 2), "%); genes lost (non-matrisome) = ", nrow(a) - nrow(m))
  lg(g, " C: per Division: ", paste(names(tdiv), tdiv, sep = "=", collapse = "; "))
  lg(g, " C: per Category: ", paste(names(tcat), tcat, sep = "=", collapse = "; "))
  ANSTAT[[g]] <- list(n_in = nrow(a), n_mat = nrow(m), div = tdiv, cat = tcat)
  ANN[[g]] <- a; MAT[[g]] <- m

  ## ---- STEP D
  for (mode in c("A", "B")) {
    d  <- call_degs(m, mode); f <- fmt_deg(d)
    write.csv(f, sprintf("results/tables/%s_matrisome_DEGs_mode%s.csv", g, mode), row.names = FALSE)
    DEG[[paste(g, mode)]] <- f
    nu <- sum(f$Direction == "Up"); nd <- sum(f$Direction == "Down")
    nc <- sum(f$Division == "Core matrisome"); na_ <- sum(f$Division == "Matrisome-associated")
    SUMM[[paste(g, mode)]] <- data.frame(dataset = g, mode = mode, n_up = nu, n_down = nd, n_total = nrow(f), n_core = nc, n_associated = na_, stringsAsFactors = FALSE)
    lg(g, " D mode ", mode, ": matrisome rows ", nrow(m), " -> DEGs ", nrow(f), " (Up ", nu, ", Down ", nd, "; Core ", nc, ", Assoc ", na_, ")")
    if (mode == "A" && nrow(f) == 0) warnf(g, ": ZERO matrisome DEGs in mode A (adj.P.Val<0.05 & |logFC|>0.585). Continuing with mode B.")
  }
}

## ---------------------------------------------------------------- STEP D table
S <- do.call(rbind, SUMM); rownames(S) <- NULL
write.csv(S, "results/tables/matrisome_DEG_summary.csv", row.names = FALSE)
cmp <- data.frame(dataset = GSES,
                  modeA_up = S$n_up[S$mode == "A"][match(GSES, S$dataset[S$mode == "A"])], modeA_down = S$n_down[S$mode == "A"][match(GSES, S$dataset[S$mode == "A"])],
                  modeA_total = S$n_total[S$mode == "A"][match(GSES, S$dataset[S$mode == "A"])],
                  modeB_up = S$n_up[S$mode == "B"][match(GSES, S$dataset[S$mode == "B"])], modeB_down = S$n_down[S$mode == "B"][match(GSES, S$dataset[S$mode == "B"])],
                  modeB_total = S$n_total[S$mode == "B"][match(GSES, S$dataset[S$mode == "B"])], thesis = as.integer(THESIS[GSES]))
cmp$diff_modeB_vs_thesis <- cmp$modeB_total - cmp$thesis; cmp$diff_modeA_vs_thesis <- cmp$modeA_total - cmp$thesis
cat("\n=== dataset x mode x (Up, Down, Total) vs thesis ===\n"); print(cmp, row.names = FALSE)
lg("D: counts table:\n", paste(capture.output(print(cmp, row.names = FALSE)), collapse = "\n"))
write.csv(cmp, "results/tables/matrisome_DEG_counts_vs_thesis.csv", row.names = FALSE)

## ---------------------------------------------------------------- STEP E: plots + matrianalyze
for (g in GSES) {
  a <- ANN[[g]]; m <- MAT[[g]]
  for (mode in c("A", "B")) {
    degs <- DEG[[paste(g, mode)]]
    if (nrow(degs) == 0) { lg(g, " E mode ", mode, ": no DEGs, no figures"); next }
    ## same thresholds on ALL annotated genes -> provides the real Non-matrisome DEG count that matristar() needs
    pa <- if (mode == "A") a$adj.P.Val else a$P.Value
    aall <- a[!is.na(pa) & pa < PCUT & abs(a$logFC) > LFC, ]; aall$Direction <- ifelse(aall$logFC > LFC, "Up", "Down")
    for (grp in c("All", "Up", "Down")) {
      sel <- if (grp == "All") rep(TRUE, nrow(aall)) else aall$Direction == grp
      sub_all <- keepattr(aall[sel, ]); sub_mat <- keepattr(sub_all[sub_all$`Annotated Matrisome Division` != "Non-matrisome", ])
      tag <- sprintf("%s_mode%s_%s", g, mode, grp); ttl <- sprintf("%s | mode %s | %s DEGs (n=%d matrisome)", g, mode, grp, nrow(sub_mat))
      if (nrow(sub_mat) < 1) { lg(g, " E ", tag, ": 0 matrisome DEGs, skipped"); next }
      for (fn in c("matribar", "matriring", "matriflow", "matristar")) {
        dat <- if (fn == "matristar") sub_all else sub_mat   # matristar fails without Non-matrisome rows (inspected source)
        res <- tryCatch({ p <- get(fn)(dat, print.plot = FALSE); if (!inherits(p, "gg")) stop("did not return ggplot: ", class(p)[1])
                          save_plot(p + labs(title = ttl), file.path("results/figures", paste0(tag, "_", fn)), if (fn == "matristar") 10 else 9, if (fn == "matristar") 10 else 7); "ok" },
                        error = function(e) paste("FAILED:", conditionMessage(e)))
        lg(g, " E ", tag, " ", fn, " (input rows ", nrow(dat), "): ", res)
        FIGLOG[[length(FIGLOG) + 1]] <- data.frame(dataset = g, mode = mode, group = grp, plot = fn, n_input = nrow(dat), status = res)
        if (res != "ok") warnf(g, " ", tag, " ", fn, " ", res)
      }
      ## matrianalyze: package SUMS numeric columns per Division/Category; needs >=2 numeric columns and strips non-numeric chars
      ## (so p-values in scientific notation would be corrupted) -> use only logFC and AveExpr.
      ma <- sub_mat[, c("Annotated Gene", "Annotated Matrisome Division", "Annotated Matrisome Category", "logFC", "AveExpr")]
      ma <- keepattr(ma)
      r <- tryCatch(matrianalyze(ma), error = function(e) paste("FAILED:", conditionMessage(e)))
      if (is.data.frame(r)) { write.csv(r, sprintf("results/tables/%s_matrianalyze.csv", tag), row.names = FALSE); lg(g, " E ", tag, " matrianalyze: ", nrow(r), " rows saved") }
      else { lg(g, " E ", tag, " matrianalyze ", r); warnf(g, " ", tag, " matrianalyze ", r) }
    }
  }
}
if (length(FIGLOG)) write.csv(do.call(rbind, FIGLOG), "logs/figure_status.csv", row.names = FALSE)

## ---------------------------------------------------------------- STEP F
for (g in GSES) writeLines(DEG[[paste(g, "A")]]$gene, sprintf("results/tables/%s_matrisome_DEG_genes.txt", g))
lg("F: wrote mode-A gene lists (n = ", paste(sapply(GSES, function(g) nrow(DEG[[paste(g, "A")]])), collapse = ", "), ")")
preview <- function(mode) {
  L <- lapply(GSES, function(g) DEG[[paste(g, mode)]]); names(L) <- GSES
  com <- Reduce(intersect, lapply(L, `[[`, "gene"))
  out <- data.frame(gene = com, stringsAsFactors = FALSE)
  for (g in GSES) { i <- match(com, L[[g]]$gene); out[[paste0("logFC_", g)]] <- L[[g]]$logFC[i]; out[[paste0("Direction_", g)]] <- L[[g]]$Direction[i] }
  dc <- paste0("Direction_", GSES)
  out$Division <- L[[1]]$Division[match(com, L[[1]]$gene)]; out$Category <- L[[1]]$Category[match(com, L[[1]]$gene)]
  out$direction_consistent <- if (nrow(out)) apply(out[, dc, drop = FALSE], 1, function(r) length(unique(r)) == 1) else logical(0)
  out
}
pvA <- preview("A"); pvB <- preview("B")
write.csv(pvA, "results/tables/common_matrisome_DEGs_preview.csv", row.names = FALSE)
write.csv(pvB, "results/tables/common_matrisome_DEGs_preview_modeB.csv", row.names = FALSE)
lg("F: PREVIEW (not final) mode A intersection = ", nrow(pvA), "; direction-consistent = ", sum(pvA$direction_consistent))
lg("F: PREVIEW (not final) mode B intersection = ", nrow(pvB), "; direction-consistent = ", sum(pvB$direction_consistent))
if (any(S$n_total[S$mode == "A"] == 0)) warnf("At least one dataset has zero mode-A DEGs, so the mode-A intersection is empty by construction.")

## ---------------------------------------------------------------- MATRISOME_SUMMARY.md
md <- function(df) { df <- as.data.frame(df); c(paste0("| ", paste(colnames(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
                                                apply(df, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |"))) }
cl_tab <- do.call(rbind, lapply(GSES, function(g) { s <- CLSTAT[[g]]; data.frame(dataset = g, rows_in = s$rows_in, slash_artefact_fixed = s$slash_stripped, empty_or_NA_removed = s$removed_empty, nonfinite_stats_removed = s$removed_na_stats,
        multigene_rows_split = s$multi_rows, rows_after_split = s$rows_after_split, duplicates_collapsed = s$dup_removed, unique_genes_out = s$rows_out) }))
an_tab <- do.call(rbind, lapply(GSES, function(g) { s <- ANSTAT[[g]]; data.frame(dataset = g, genes_input = s$n_in, matrisome = s$n_mat, pct = round(100 * s$n_mat / s$n_in, 2),
        Core = unname(s$div["Core matrisome"]), Associated = unname(s$div["Matrisome-associated"]), Glycoproteins = unname(s$cat["ECM Glycoproteins"]), Collagens = unname(s$cat["Collagens"]),
        Proteoglycans = unname(s$cat["Proteoglycans"]), Affiliated = unname(s$cat["ECM-affiliated Proteins"]), Regulators = unname(s$cat["ECM Regulators"]), Secreted = unname(s$cat["Secreted Factors"])) }))
pvtxt <- function(p, nm) c(sprintf("- %s: **%d** genes common to all three datasets; direction consistent in **%d**.", nm, nrow(p), sum(p$direction_consistent)),
                           if (nrow(p)) paste0("  - ", p$gene, ": ", apply(p[, paste0("Direction_", GSES)], 1, paste, collapse = "/"), collapse = "\n"))
sm <- c("# Matrisome DEG summary (MatrisomeAnalyzeR)", "",
        sprintf("Generated by `scripts/05_matrisome_annotate.R` on %s. Package: MatrisomeAnalyzeR v%s (GitHub %s/%s @ %s). Seed 20261102. Species: human.", format(Sys.time(), "%Y-%m-%d %H:%M"), pv, pd$GithubUsername, pd$GithubRepo, sha), "",
        "Thresholds: |logFC| > 0.585; mode A = adj.P.Val < 0.05, mode B = P.Value < 0.05.", "",
        "## 1. Symbol cleaning (Step B)", md(cl_tab), "", "## 2. Matrisome annotation (Step C)", md(an_tab), "",
        "## 3. Matrisome DEGs vs thesis (Step D)", md(cmp), "",
        "## 4. Preview of intersection (not final)", pvtxt(pvA, "Mode A"), pvtxt(pvB, "Mode B"), "",
        "## 5. Warnings", if (length(WARN)) paste0("- ", WARN) else "- none", "",
        "## 6. Files", "- results/tables: `<GSE>_matriannotated_all.csv`, `<GSE>_matrisome_all_genes.csv`, `<GSE>_matrisome_DEGs_modeA/B.csv`, `<GSE>_matrisome_DEG_genes.txt`, `<GSE>_mode?_<All|Up|Down>_matrianalyze.csv`, `matrisome_DEG_summary.csv`, `matrisome_DEG_counts_vs_thesis.csv`, `common_matrisome_DEGs_preview*.csv`",
        "- results/figures: `<GSE>_mode?_<All|Up|Down>_<matribar|matriring|matriflow|matristar>.png/.pdf`", "- logs: `matrisome_step.log`, `sessionInfo.txt`, `figure_status.csv`")
writeLines(unlist(sm), "results/MATRISOME_SUMMARY.md")
lg("DONE. warnings: ", length(WARN))

