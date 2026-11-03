# Pre-Specified 15-Gene List Evaluation Protocol & Documentation

**Project**: Idiopathic Pulmonary Fibrosis (IPF) Matrisome Analysis  
**R Version**: 4.6.2 (Windows ucrt)  
**Fixed Random Seed**: `20261102`  
**Execution Date**: 2026-11-03  
**Log File**: `logs/pre_specified_test.log`  
**Session Info**: `logs/sessionInfo_prespecified.txt`  

---

## 1. Candidate Gene List (Pre-Specified from Thesis)
The following 15 candidate matrisome genes were fixed in advance and evaluated without any post-hoc additions, exclusions, or criteria tuning:
1. `ASPN` (Asporin)
2. `COL14A1` (Collagen Type XIV Alpha 1 Chain)
3. `COL15A1` (Collagen Type XV Alpha 1 Chain)
4. `COL1A1` (Collagen Type I Alpha 1 Chain)
5. `COL3A1` (Collagen Type III Alpha 1 Chain)
6. `COL5A2` (Collagen Type V Alpha 2 Chain)
7. `COMP` (Cartilage Oligomeric Matrix Protein)
8. `IL6` (Interleukin 6)
9. `MMP7` (Matrix Metallopeptidase 7)
10. `MMP8` (Matrix Metallopeptidase 8)
11. `MMP10` (Matrix Metallopeptidase 10)
12. `OGN` (Osteoglycin)
13. `SFRP4` (Secreted Frizzled Related Protein 4)
14. `TGFB3` (Transforming Growth Factor Beta 3)
15. `TIMP3` (TIMP Metallopeptidase Inhibitor 3)

---

## 2. Pre-Specified Evaluation Criteria (Fixed in Advance)

### Part 1: WGCNA Status Criteria
1. **Network**: Top 8,000 MAD genes from GSE32539 ($N=168$; outlier `GSM806411` removed). Signed network, soft threshold power $\beta = 8$, minimum module size = 30, merge cut height = 0.25.
2. **Key Module Definition**: Module–trait correlation $|r| \ge 0.50$ AND Benjamini–Hochberg adjusted $P < 0.05$ with IPF phenotype ($1 = \text{IPF}, 0 = \text{Control}$).
3. **Strict Hub Definition**: Assigned to a Key Module AND Module Membership $|MM| \ge 0.80$ AND Gene Significance $|GS| \ge 0.20$.
4. **Mutually Exclusive Verdicts**:
   - `Hub (strict)`: Key module AND $|MM| \ge 0.80$ AND $|GS| \ge 0.20$.
   - `Key module, below hub cutoff`: Key module, but fails $|MM| < 0.80$ or $|GS| < 0.20$.
   - `Non-key module`: Assigned to a module with $|r| < 0.50$ or adj. $P \ge 0.05$.
   - `Grey / unassigned`: Module 0 (unassigned genes).
   - `Not in network`: Not present in top 8,000 MAD gene set.
5. **Exploratory Relaxed Verdict** (Reported separately; does not alter strict verdict):
   - `Hub (relaxed)`: Key module AND ($|MM| \ge 0.70$ OR intramodular connectivity $kWithin$ in top 20% of its module).

### Part 2: Single-Gene Diagnostic Performance Criteria
1. **Datasets**:
   - `GSE32539`: Discovery training cohort ($N=168$: 119 IPF, 49 Control; whole lung tissue).
   - `GSE110147`: Discovery replication cohort ($N=33$: 22 IPF, 11 Control; whole lung tissue).
   - `GSE24206`: Discovery cohort ($N=15$: 9 IPF, 6 Control; lung tissue).
   - `GSE135065`: Strictly independent validation cohort ($N=18$: 9 IPF, 9 Control; Affymetrix HuGene-2.1 ST array, primary human lung fibroblasts).
2. **Metrics**:
   - Direction in IPF: `UP` ($\log FC > 0$) or `DOWN` ($\log FC < 0$).
   - Direction-free AUC: $\max(AUC, 1 - AUC)$ with signed direction ($+$ or $-$).
   - DeLong 95% Confidence Intervals (`pROC`).
   - Wilcoxon rank-sum test $P$-value with Benjamini–Hochberg FDR adjustment across the 15 genes.
   - 4-Cohort significance: Significant ($adj.P < 0.05$) in all 4 datasets.
   - Direction consistency: Identical direction across all 4 datasets.

### Part 3 & 4: Multi-Gene Classifier & Feature Importance Criteria
1. **Standardization**: Within-cohort $z$-score normalization to preserve relative distribution without inter-cohort leakage.
2. **Model 1 (Ridge Logistic Regression)**: L2-regularized binomial regression (`glmnet`, $\alpha = 0$), inner 5-fold cross-validation for optimal $\lambda$.
3. **Model 2 (Unweighted Directional Score)**: Sum of standardized expression weighted by the training logFC sign: $Score_i = \sum_{g=1}^{15} \text{sign}(\log FC_g) \cdot z_{gi}$.
4. **Internal Validation**: Repeated stratified 5-fold cross-validation $\times$ 20 repeats (100 folds) on GSE32539.
5. **External Evaluation**: Models trained on full GSE32539 ($N=168$) and tested on GSE110147 ($N=33$), GSE24206 ($N=15$), GSE135065 ($N=18$), and GSE135065 untreated baseline ($N=6$).
6. **Empirical Baselines (1,000 sets each)**:
   - 1,000 random sets of 15 genes from Network (8,000 genes).
   - 1,000 random sets of 15 genes from Matrisome (639 genes).
   - 1,000 random sets of 15 genes from Shared 45 DEGs.
   - 200 label permutations of GSE32539 phenotype.
7. **Leave-One-Gene-Out Sensitivity**: Marginal impact $\Delta AUC = AUC_{-g} - AUC_{all15}$ on GSE135065.
8. **Feature Importance**: Random Forest (500 trees, Mean Decrease Accuracy and Mean Decrease Gini) and Linear SVM feature weights ($|w|$).

### Part 5: Final Pre-Specified Evidence Tiers
- **Strong**: Strict WGCNA Hub AND Significant in all 4 datasets AND Consistent direction in all 4 datasets AND Independent AUC $\ge 0.80$.
- **Moderate**: Strict WGCNA Hub OR (Significant in all 4 datasets with consistent direction AND Independent AUC $\ge 0.80$), but NOT both.
- **Weak / Not Supported**: Fails both conditions.

---

## 3. Directory Structure and Deliverables

```
project/
├── logs/
│   ├── pre_specified_test.log           # Complete execution log with timestamps
│   └── sessionInfo_prespecified.txt     # R sessionInfo()
├── scripts/
│   ├── 08_prespecified_wgcna.R          # Part 1: WGCNA status & robustness
│   ├── 09_prespecified_auc.R            # Part 2: Single-gene AUC & boxplots
│   ├── 10_prespecified_classifier.R     # Part 3 & 4: Multi-gene ML & baselines
│   └── 11_prespecified_final_table.R    # Part 5: Final synthesis table assembly
├── results/
│   ├── PRESPECIFIED_SUMMARY.md          # Comprehensive final narrative report
│   ├── tables/prespecified/
│   │   ├── wgcna_status_15.csv          # WGCNA metrics, module membership, verdicts
│   │   ├── single_gene_auc_15.csv       # Single-gene AUCs, CIs, Wilcoxon p per dataset
│   │   ├── classifier_performance.csv   # Multi-gene classifier CV & external metrics
│   │   ├── random_set_comparison.csv    # 1,000 random sets & permutation p-values
│   │   ├── leave_one_out.csv            # Leave-one-gene-out delta AUCs
│   │   ├── feature_importance_15.csv    # RF, SVM, and Ridge importance rankings
│   │   └── final_15_gene_table.csv      # Unified 15-gene synthesis & evidence tiers
│   └── figures/prespecified/
│       ├── MM_vs_GS_green_15genes.png / .pdf
│       ├── MM_vs_GS_turquoise_15genes.png / .pdf
│       ├── kWithin_rank_green_15genes.png / .pdf
│       ├── heatmap_15genes_wgcna.png / .pdf
│       ├── single_gene_auc_forest.png / .pdf
│       ├── expression_boxplots_part1.png / .pdf
│       ├── expression_boxplots_part2.png / .pdf
│       ├── classifier_roc_curves.png / .pdf
│       ├── classifier_null_distributions.png / .pdf
│       ├── leave_one_out_delta_auc.png / .pdf
│       └── feature_importance_comparison.png / .pdf
```

---

## 4. Notable Biological and Methodological Observations
1. **Independent Cohort Biology (GSE135065)**:
   - GSE135065 is composed of primary human lung fibroblasts cultured in vitro, including untreated baseline cells and cells treated with TGF-$\beta$ (24h and 5d).
   - In untreated baseline conditions ($N=6$: 3 IPF vs 3 Control), the 15-gene classifier achieved perfect separation (AUC = 1.000, Balanced Accuracy = 1.000).
   - Across all 18 samples, normal control fibroblasts stimulated with TGF-$\beta$ acquired an intense myofibroblast/fibrotic gene expression state, altering unadjusted whole-population AUCs.
2. **Direction Inconsistency in In Vitro vs In Vivo Fibrosis**:
   - Collagen and proteoglycan genes (`COL1A1`, `COL3A1`, `COL5A2`, `COL14A1`, `ASPN`, `OGN`) are robustly upregulated in IPF whole lung tissue biopsies (due to massive ECM accumulation and myofibroblast infiltration). In cultured 2D fibroblasts, normal fibroblasts express higher basal levels or different response kinetics on plastic.
   - Conversely, `COL15A1`, `MMP7`, `MMP10`, `SFRP4`, and `TIMP3` preserved concordant expression direction between tissue biopsies and cultured fibroblasts.

