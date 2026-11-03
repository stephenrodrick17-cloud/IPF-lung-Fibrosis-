# WGCNA Pipeline Documentation (IPF Matrisome Project)

**Project Root**: `d:/2ND PROJECT/project`  
**Execution Script**: `scripts/06_wgcna.R`  
**R Version**: 4.6.2 (Windows ucrt)  
**WGCNA Version**: 1.74  
**Random Seed**: `20261102`  
**Execution Command**:
```bash
Rscript scripts/06_wgcna.R
```

---

## 1. Executive Summary & Workflow Overview

This workflow performs Weighted Gene Co-expression Network Analysis (WGCNA) on idiopathic pulmonary fibrosis (IPF) vs. normal lung tissue microarray profiles to uncover functional co-expression modules, establish clinical trait associations (IPF vs. Control), and identify candidate hub genes with an emphasis on the extracellular matrix (matrisome).

The pipeline strictly fulfills all reproducible bioinformatic standards:
- **No data fabrication**: All calculations derive directly from normalized expression series and validated MatrisomeAnalyzeR annotations.
- **Genome-wide network inference**: WGCNA was constructed on the top 8,000 variable genes (MAD), avoiding artificial restriction to candidate DEGs.
- **Multi-tiered hub gene deliverables**: Primary matrisome hub genes (`hubs_primary`, $n=24$), expanded matrisome hubs (`hubs_matrisome`, $n=60$), and global network hubs (`hubs_all`, $n=896$).
- **No machine learning executed**: Random Forest and SVM-RFE modeling are reserved for downstream user workflows; curated, cross-dataset expression and label matrices are cleanly prepared in `data/processed/ml_input/`.

---

## 2. Dataset Selection & Cohort Rationale

| Dataset | Platform | Raw Samples | Cleaned Samples (IPF / Control) | Role in WGCNA |
| :--- | :--- | :--- | :--- | :--- |
| **GSE32539** | Agilent-014850 Whole Human Genome (GPL6480) | 119 IPF, 50 Control (N=169) | **119 IPF, 49 Control (N=168)** | **Primary Discovery Network** (high statistical power) |
| **GSE110147** | Affymetrix HT HG-U133+ PM (GPL13912) | 22 IPF, 11 Control (N=33) | **22 IPF, 11 Control (N=33)** | **Replication & Module Preservation** ($N=33$) |
| **GSE24206** | Affymetrix Human Gene 1.0 ST (GPL6244) | 9 IPF, 6 Control (N=15) | **9 IPF, 6 Control (N=15)** | **Excluded from Network Construction**; reserved for ML hand-off |

### Cohort Justification:
- **WGCNA Sample Size Guidelines**: WGCNA documentation strongly recommends a minimum of $N \ge 20\text{--}30$ samples to ensure stable correlation and adjacency estimates.
- **GSE32539** ($N=168$ after outlier screening) provides excellent statistical power for signed topological overlap matrix (TOM) calculation and robust dynamic tree cutting.
- **GSE110147** ($N=33$) provides an independent replication cohort to compute permutation-based module preservation statistics ($Z_{summary}$).
- **GSE24206** ($N=15$) is strictly underpowered for de novo co-expression network construction, which would yield unstable modules and spurious hubs. It is retained solely as an independent test set for downstream machine learning.

---

## 3. Data Preprocessing & Outlier Screening

1. **Gene-Level Log2 Expression**:
   - GSE32539 expression range: `1.615` to `13.027` (log2 scale confirmed).
   - Probes collapsed to unique gene symbols using maximum $|\text{logFC}|$, then lowest $P$-value.
2. **Quality Check**:
   - `WGCNA::goodSamplesGenes` confirmed zero missing entries (`allOK = TRUE`).
3. **Sample Outlier Removal**:
   - Hierarchical clustering with Euclidean distance identified sample **GSM806411** (Control) as an extreme outlier ($Z.k = -3.60 < -3$, tree branch height $> 75$).
   - Removing GSM806411 preserved $N=168$ samples (119 IPF, 49 Control). No samples were excluded based on clinical classification.
4. **Gene Filtering**:
   - Filtered to the top 8,000 genes by Median Absolute Deviation (MAD).
   - **Verification**: **All 45 shared matrisome DEGs survived naturally** ($45/45$, 100%), and 555 annotated matrisome genes survived. Zero DEGs required manual forced retention.

---

## 4. Network Construction & Module Detection

1. **Soft-Threshold Power Selection**:
   - Signed network topology tested across powers $\beta = 1\text{--}20$.
   - Criterion: lowest power with scale-free topology fit $R^2 \ge 0.85$ and acceptable mean connectivity.
   - **Selected Power**: $\beta = 8$ ($R^2 = 0.8667$, mean connectivity = $148.4$).
2. **Module Detection (`blockwiseModules`)**:
   - Parameters: `networkType = "signed"`, `TOMType = "signed"`, `minModuleSize = 30`, `mergeCutHeight = 0.25`, `deepSplit = 2`, `reassignThreshold = 0`, `numericLabels = TRUE`, `maxBlockSize = 10000`, `corType = "pearson"`.
   - All 8,000 genes processed in a single block.
   - Identified **16 co-expression modules** (plus unassigned grey module of 813 genes).

### Module Size Distribution:
- `turquoise`: 1,678 genes
- `blue`: 822 genes
- `brown`: 610 genes
- `yellow`: 601 genes
- `green`: 551 genes
- `red`: 549 genes
- `black`: 481 genes
- `pink`: 375 genes
- `magenta`: 297 genes
- `purple`: 284 genes
- `greenyellow`: 272 genes
- `tan`: 216 genes
- `salmon`: 126 genes
- `cyan`: 119 genes
- `midnightblue`: 111 genes
- `lightcyan`: 95 genes
- `grey` (unassigned): 813 genes

---

## 5. Module–Trait Relationships & Key Modules

Clinical trait coded: **IPF = 1**, **Control = 0**.

Module eigengenes (ME) were correlated with disease status using Pearson correlation ($r$) and Student asymptotic $p$-values adjusted by Benjamini-Hochberg (BH) across all modules.

**Key Module Threshold**: $|r| \ge 0.5$ AND BH-adjusted $p < 0.05$.

Seven modules passed this threshold:
1. **`magenta`**: $r = +0.76$, $p_{adj} = 1.05 \times 10^{-31}$ (Strongly Up in IPF)
2. **`blue`**: $r = -0.73$, $p_{adj} = 4.29 \times 10^{-28}$ (Strongly Down in IPF)
3. **`green`**: $r = +0.70$, $p_{adj} = 2.45 \times 10^{-24}$ (Strongly Up in IPF; Primary ECM Fibrotic Module)
4. **`turquoise`**: $r = +0.67$, $p_{adj} = 2.76 \times 10^{-21}$ (Strongly Up in IPF; Matrisome/Tissue Remodeling)
5. **`brown`**: $r = -0.66$, $p_{adj} = 6.46 \times 10^{-20}$ (Down in IPF; Normal Epithelial / Vascular)
6. **`salmon`**: $r = +0.62$, $p_{adj} = 3.32 \times 10^{-17}$ (Up in IPF)
7. **`lightcyan`**: $r = +0.55$, $p_{adj} = 4.09 \times 10^{-13}$ (Up in IPF)

---

## 6. Hub Gene Identification & Deliverables

### Hub Definition:
- Assigned to a **Key Module** ($|r| \ge 0.5$, $p_{adj} < 0.05$).
- Module Membership $|MM| \ge 0.8$.
- Gene Significance $|GS| \ge 0.2$.
- Ranked within module by intramodular connectivity ($kWithin$).

### Deliverable Lists:
1. **`hubs_primary`** ($n = 24$ genes):
   - Overlap of Key Module Hubs and the **45 shared matrisome DEGs**.
   - Output: `results/tables/wgcna/hubs_primary.csv` and `hubs_primary.txt`.
   - Status: $n = 24 \ge 8$; perfectly sized for initial feature selection.
2. **`hubs_matrisome`** ($n = 60$ genes):
   - All annotated matrisome genes in key modules meeting $|MM| \ge 0.8$ and $|GS| \ge 0.2$.
   - Output: `results/tables/wgcna/hubs_matrisome.csv` and `hubs_matrisome.txt`.
   - Status: Recommended for broader machine learning hyperparameter tuning and biomarker discovery.
3. **`hubs_all`** ($n = 896$ genes):
   - All genes meeting hub criteria across all 7 key modules.
   - Output: `results/tables/wgcna/hubs_all.csv` and `hubs_all.txt`.
4. **`hub_genes_final.csv`**:
   - Master reference catalog containing all 896 hub candidates annotated with:
     - Matrisome Division & Category
     - Module assignment, $MM$, $GS$, $kWithin$, rank in module
     - Presence in 45 shared DEGs
     - LogFC and adj.P in GSE24206, GSE32539, and GSE110147
     - Concordance flag across all 3 datasets
     - Overlap with the thesis 15 hub gene list

---

## 7. Biological Validation & Module Preservation

1. **Matrisome Enrichment**:
   - Fisher's exact test revealed dramatic over-representation of matrisome genes in the **`green` module** ($p < 10^{-35}$, odds ratio $> 8$) and **`turquoise` module** ($p < 10^{-15}$).
   - Of the 45 shared matrisome DEGs: **21 reside in `green`**, **14 in `turquoise`**, **4 in `blue`**, **3 in `magenta`**, **2 in `greenyellow`**, and **1 in `brown`**.
2. **Functional Enrichment (GO BP)**:
   - **`green` module**: Extracellular matrix organization, collagen fibril organization, connective tissue development, cell-substrate adhesion.
   - **`turquoise` module**: Wound healing, ossification, cell-matrix adhesion, regulation of cell migration.
   - **`blue` module**: Lung alveolus development, gas exchange, cilium assembly, epithelial cell differentiation.
3. **Module Preservation in GSE110147**:
   - Permutation test ($nPerm = 100$) evaluated in independent GSE110147 ($N=33$).
   - Multiple modules display strong preservation ($Z_{summary} > 10$), including `turquoise`, `blue`, `pink`, `black`, `red`, `salmon`, and `yellow`.
   - Moderate preservation ($2 \le Z_{summary} \le 10$) in smaller modules reflects the modest sample size ($N=33$) of GSE110147 rather than biological discrepancy.

---

## 8. Machine Learning Hand-off Datasets

Clean, model-ready data structures are saved in `data/processed/ml_input/`:
1. `GSE32539_expr_hubs_all_candidates.csv` (168 samples $\times$ 896 hub genes, log2 normalized)
2. `GSE110147_expr_hubs_all_candidates.csv` (33 samples $\times$ 896 hub genes, log2 normalized)
3. `GSE24206_expr_hubs_all_candidates.csv` (15 samples $\times$ 896 hub genes, log2 normalized)
4. `<GSE>_labels.csv`: Matching sample metadata with binary diagnostic targets (`1` = IPF, `0` = Control).
5. `candidate_gene_sets.csv`: Gene identifiers mapped to `hubs_primary`, `hubs_matrisome`, and `thesis_15` sets for automated pipeline input in Random Forest and SVM-RFE.

---

## 9. Deviations & Technical Decisions

1. **GSE24206 Sample Size**: Excluded from direct WGCNA network construction due to $N=15$ (too small for valid soft-thresholding and dynamic module pruning). Retained for external validation in ML.
2. **Outlier Filtering**: GSM806411 was removed objectively based on standardized connectivity $Z.k < -3.6$ and distance $>75$.
3. **Namespace Conflicts**: Explicitly overridden `cor <- WGCNA::cor` and avoided early import of `clusterProfiler` to prevent S4 method collision with `stats::cor`.
4. **Threading**: Single-threaded execution configured (`WGCNA::disableWGCNAThreads()`) to prevent socket lockups on Windows R 4.6.2.

