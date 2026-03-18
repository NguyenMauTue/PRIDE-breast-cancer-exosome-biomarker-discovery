# PRIDE Breast Cancer Exosome Biomarker Discovery

Bioinformatics pipeline for prioritizing candidate protein biomarkers 
in breast cancer–derived exosomes using public LFQ proteomics data.

---

## Biological Background

Extracellular vesicles (exosomes) play central roles in cancer biology,
mediating intercellular communication, metastasis, and immune modulation.
Proteins packaged into tumor-derived exosomes may serve as **non-invasive 
biomarkers detectable through liquid biopsy**.

This project investigates exosomal proteomic differences between breast 
cancer and normal breast epithelial cells, and prioritizes candidate 
proteins using a multi-criteria scoring approach. Given the exploratory 
nature of the dataset (n=3 per condition), results are framed as 
**hypothesis-generating candidate discovery** pending independent validation.

---

## Dataset

| Field | Details |
|---|---|
| Dataset ID | PXD056161 |
| PRIDE URL | https://www.ebi.ac.uk/pride/archive/projects/PXD056161 |
| DOI | 10.1021/acs.jproteome.4c00795 |

**Cell models:**
- MCF-10A — non-tumorigenic breast epithelial cells (Normal)
- MDA-MB-231 — triple-negative breast cancer cells (Tumor)

**Biological material:** Extracellular vesicles (exosomes)

**Proteomics workflow:**
- Bottom-up proteomics, label-free quantification (LFQ)
- Orbitrap Fusion mass spectrometer
- MaxQuant/Andromeda search engine, Human UniProt reference proteome

---

## Analysis Workflow
```
01_load_and_filter_data.R         # MaxQuant QC filtering
02_metadata_and_matrix.R          # LFQ matrix construction, log2 transform
03_quality_control.R              # Technical replicate PCA, correlation QC
04_missingness_analysis.R         # MNAR characterization, logistic regression
05_imputation.R                   # Left-shifted Gaussian imputation
06_differential_expression.R      # limma DE analysis
07_robustness_analysis.R          # Imputation sensitivity, complete-case comparison
08_pathway_enrichment_analysis.R  # Reactome GSEA
09_pathway_thematic_grouping.R    # Unbiased gene list export for STRING
10_protein_protein_interaction.R  # STRING network construction (score > 0.7)
11_network_expression_integration.R  # Merge DE results with network topology
12_CDS_via_AHP.R                  # AHP-weighted Composite Driver Score + sensitivity analysis
13_functional_annotation.R        # biomaRt GO annotation
14_biological_theme_classification.R  # Post-hoc theme annotation, localization, PubMed mining
15_final_tables_and_visualization.R   # Final candidates, plots, Excel export
```

---

## Methods Summary

### Composite Driver Score (CDS)
Candidate proteins were ranked using an AHP-weighted composite score 
integrating four criteria:

| Criterion | Weight | Rationale |
|---|---|---|
| Fold Change (FC) | ~0.52 | Primary signal for clinical detectability |
| FDR | ~0.16 | Downweighted due to limited statistical power (n=3) |
| Betweenness centrality | ~0.08 | Network topology — supporting evidence |
| Degree | ~0.24 | Network connectivity — supporting evidence |

Consistency Ratio (CR) < 0.1, confirming internal consistency of the 
pairwise judgment matrix.

All features were normalized to [0,1] prior to scoring. No statistical 
pre-filtering was applied given the n=3 sample size per condition.

### Sensitivity Analysis
AHP weights were perturbed ±20% across 5 levels per criterion. Proteins 
with stable rank across all perturbations are labeled **robust_candidate**; 
proteins with high rank variance (rank_sd > 1.5) are labeled 
**weight_sensitive_candidate**.

### Localization Filter
Proteins with exclusively intracellular or nuclear localization (GO_CC) 
were excluded as likely co-isolation contaminants, consistent with 
exosome cargo biology. Histone proteins were additionally excluded based 
on protein description annotation.

---

## Key Findings

The prioritization pipeline identified **47 candidate proteins** with 
exosome-compatible localization profiles (extracellular, membrane, 
vesicle/exosome compartments).

**Top candidates by CDS:**

| Rank | Symbol | Localization | Robustness |
|---|---|---|---|
| 1 | CD47 | membrane | robust |
| 2 | SDC1 | extracellular | robust |
| 3 | NRAS | extracellular | robust |
| 4 | LAMB2 | extracellular | robust |
| 5 | CD36 | extracellular | robust |
| 6 | TUBA1A | membrane | robust |
| 7 | RAB7A | membrane | robust |
| 8 | ACTB | extracellular | robust |
| 9 | GRB2 | extracellular | robust |
| 10 | FN1 | extracellular | robust |

**Novel high-priority candidates** (high CDS, low PubMed hits):
- TUBA1A
- LAMB2

CD47 (rank #1) has independent clinical validation in breast cancer 
exosome cohorts, supporting the validity of the CDS prioritization approach.

---

## Software and Tools

Analysis performed in R (v4.5.2) using:

| Package | Purpose |
|---|---|
| data.table, dplyr | Data manipulation |
| ggplot2, ggrepel | Visualization |
| limma | Differential expression |
| igraph | Network analysis |
| clusterProfiler, ReactomePA | Pathway enrichment |
| biomaRt | Functional annotation |
| rentrez | PubMed mining |
| purrr | Sensitivity analysis |
| openxlsx | Excel export |

---

## Reproducibility

All scripts are available in the `scripts/` directory and are designed 
to be run sequentially. Raw proteomics data can be downloaded directly 
from PRIDE using dataset identifier **PXD056161** — script 01 loads 
the data automatically via FTP.

---

## Limitations

- Sample size: n=3 per condition limits statistical power; 
  adj.P.Val thresholds were not applied as pre-filters
- Cell line model: MDA-MB-231 represents aggressive/metastatic 
  TNBC; findings require validation in early-stage clinical specimens
- Results are hypothesis-generating and require independent 
  experimental validation

---

## Author

Nguyen Mau Tue  
*Documentation assistance provided with the help of AI tools.*
