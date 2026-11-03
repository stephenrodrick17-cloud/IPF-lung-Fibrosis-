# Matrisome annotation step (05_matrisome_annotate.R)

Run from the project root (needs the package in `rlib/` or the user library):

    Rscript scripts/05_matrisome_annotate.R

Helper scripts `05a_install_check.R` (install, run with network), `05b_inspect*.R` and `05c_debug.R` (API inspection and debugging) are kept for traceability.
Package: MatrisomeAnalyzeR v1.0.1 from the `Matrisome/MatrisomeAnalyzeR` fork (SHA in `logs/sessionInfo.txt`). The fork installed fine, so the `izzilab` fallback was not used.
The first install attempt failed only because `remotes` was missing. I added it to the install list.

## Package facts found by inspection (not guessed)
- Exports: matrianalyze, matriannotate, matribar, matriflow, matriring, matrisome.list, matristar, testfile.
- `matriannotate` puts `Annotated Gene`, `Annotated Matrisome Division` and `Annotated Matrisome Category` first. It drops the original gene column and **re-sorts the rows**.
- Matching is exact and case-sensitive.
- The plot and analyze functions need `attr(,"workflow") == "matrisomeannotatoR"`. The script re-applies it after subsetting.

## Deviations / assumptions
1. **Upstream symbol artefact (GSE24206).** 1,738 symbols looked like `/ KRT17`. Script `02_deg_limma.R` splits `A /// B` on `//`, which leaves stray slashes, and these symbols cannot match the matrisome list. I stripped leading/trailing `/` and spaces in Step B, then collapsed duplicates (462 removed). I did not change the earlier script. The first symbol of those probes was lost upstream, and the stats are those of the probe.
2. **Multi-gene `A///B` rule.** Rows are split into separate rows. There were 0 such rows, so no effect.
3. **Dual annotations.** The package returns one row per (gene, category) pair, so CLC, EPO, HRG, MBP, NOV and PRG2 appear twice. I kept the raw output in `<GSE>_matriannotated_all.csv`.
   - Everything downstream uses one row per gene, saved as `<GSE>_matriannotated_unique_gene.csv`.
   - Core matrisome takes priority over Matrisome-associated, otherwise the first row returned.
   - The other annotation is kept in the `Other_annotations` column, which is an extra column in the matrisome and DEG CSVs.
4. **matristar.** It crashes if no `Non-matrisome` rows are present. It is therefore drawn from all DEGs (matrisome and non-matrisome, same thresholds). The other three plots use matrisome-only DEGs. The package scales the Non-matrisome bar down itself.
5. **matrianalyze.** It SUMS numeric columns per Division/Category and strips non-numeric characters. P-values in scientific notation would be corrupted, so only `logFC` and `AveExpr` are passed. The output is the sum of logFC and AveExpr per group, not a mean. It needs at least 2 numeric columns.
6. **Figures.** Made for modes A and B, and for All, Up and Down DEGs (72 PNG at 300 dpi and 72 PDF). Per-figure status is in `logs/figure_status.csv`.
7. **Mode A equals mode B** for GSE32539 and GSE110147. No matrisome gene has P < 0.05 and |logFC| > 0.585 but adj.P ≥ 0.05. The tables are full (26–34% of rows have adj.P ≥ 0.05), so this is not a significant-only table.
8. **Thesis comparison.** Counts differ from the thesis and I did not tune anything to match. Possible reasons: a different matrisome list or version (the thesis presumably used the Naba masterlist), and a different probe-collapse rule.
9. The seed is set (20261102). The pipeline is deterministic.

