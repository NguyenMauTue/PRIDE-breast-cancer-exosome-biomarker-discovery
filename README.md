# PRIDE Breast Cancer Exosome Biomarker Discovery

### A reproducible pipeline implementing the Composite Driver Score (CDS) framework for exosomal protein biomarker prioritisation in triple-negative breast cancer

---

> **Preprint:** [bioRxiv — link TBD]  
> **Manuscript PDF:** [link TBD]  
> **Data:** [PRIDE PXD056161](https://www.ebi.ac.uk/pride/archive/projects/PXD056161) · [PRIDE PXD012162](https://www.ebi.ac.uk/pride/archive/projects/PXD012162)

---

## 🔑 Key Idea

Most exosomal proteomics studies rank candidates by fold-change or statistical significance alone. This pipeline implements a **Composite Driver Score (CDS)** that integrates expression magnitude with network topology using an **Analytic Hierarchy Process (AHP)**, enabling prioritisation of proteins that are both differentially abundant and biologically embedded in functionally coherent networks.

Rather than filtering by arbitrary FDR thresholds — inappropriate given n=3 per condition — CDS weights each criterion by its relative contribution to clinical detectability and biological plausibility. The result is a ranked, robustness-tested candidate list designed to generate experimentally testable hypotheses.

---

## 📦 Repository Structure

```
├── scripts/          # Full analysis pipeline (01–16, run sequentially)
├── data/             # Input files (MaxQuant output, STRING network TSV)
└── results/          # Final tables, DE results, GSEA objects and all manuscript and supplementary figures
```

## ⚙️ How to Run

```r
# 1. Clone the repository
# git clone https://github.com/NguyenMauTue/PRIDE-breast-cancer-exosome-biomarker-discovery

# 2. Install dependencies (R >= 4.4)
install.packages(c("dplyr", "ggplot2", "ggrepel", "limma", "igraph",
                   "clusterProfiler", "ReactomePA", "biomaRt",
                   "rentrez", "openxlsx", "patchwork", "circlize"))
BiocManager::install(c("limma", "clusterProfiler", "ReactomePA",
                       "biomaRt", "org.Hs.eg.db", "reactome.db"))

# 3. Download raw data from PRIDE (PXD056161) — script 01 loads automatically via FTP

# 4. Run pipeline sequentially
source("scripts/run_all.R")
# or manually: 01 → 02 → ... → 16
```

> A `sessionInfo()` export is saved to `results/sessionInfo.txt` after each full run.  
> An `renv.lock` file is provided for full environment reproducibility.

---

## 📊 Outputs

| File | Description |
|---|---|
| `results/BiomarkerCandidates_final.xlsx` | Full ranked candidate table with CDS, localization, PubMed hits |
| `results/limma_network_table.csv` | DE results merged with STRING network topology |
| `results/agrn_pathway_dir_summary.csv` | AGRN satellite protein directionality by Reactome pathway |
| `figures/volcano.png` | Volcano plot — differential expression |
| `figures/DriverScore_landscape.png` | CDS landscape plot |
| `figures/PPI_network.png` | STRING network with DE integration |
| `figures/GO_enrichment.png` | GO enrichment dot plot |
| `figures/16_AGRN_satellite_directionality.png` | AGRN pathway co-protein concordance |

---

## ⚗️ Methods Summary

### Composite Driver Score (CDS)

Candidates ranked by AHP-weighted composite of four criteria:

| Criterion | Weight | Rationale |
|---|---|---|
| Fold Change | ~0.52 | Primary signal for clinical detectability |
| FDR | ~0.16 | Downweighted — limited power at n=3 |
| Degree | ~0.24 | Network connectivity |
| Betweenness centrality | ~0.08 | Network topology |

Consistency Ratio (CR) < 0.1. All features normalised to [0,1] prior to scoring.

### Sensitivity Analysis

AHP weights perturbed ±20% across 5 levels per criterion. Candidates labelled `robust_candidate` if rank is stable across all perturbations; `weight_sensitive_candidate` if rank SD > 1.5.

---

## 🔬 Key Results

CDS prioritisation converges on a coordinated **ECM/adhesion module** — integrins (ITGA2, ITGB1, ITGA3, ITGAV, ITGB4), fibronectin (FN1), and **AGRN** — as a coherent exosomal program in TNBC-derived vesicles.

**AGRN** (rank 9, CDS = 0.505) emerges as the highest-priority novel candidate: detected across all 6 samples with zero missing values, cross-dataset logFC concordance confirmed (PXD056161: +2.98; PXD012162: +3.43), and embedded within ECM proteoglycan and integrin interaction pathways showing 100% directional concordance with co-occurring proteins.

**FN1** (rank 1, CDS = 0.907) serves as a pipeline validation anchor — well-established in the exosome literature with 48 PubMed hits.

---

## ⚠️ Limitations

- **Sample size:** n=3 per condition limits statistical power; FDR thresholds not applied as pre-filters
- **Cell line model:** MDA-MB-231 represents aggressive/metastatic TNBC; findings require validation in early-stage clinical specimens
- **Single dataset primary:** Cross-dataset concordance (PXD012162) supports robustness but does not substitute for prospective validation
- All results are **hypothesis-generating** and require independent experimental validation

---

## 🗂️ Biological Background

Tumor-derived exosomes mediate intercellular communication, ECM remodelling, and immune modulation. Proteins packaged into these vesicles are detectable non-invasively via liquid biopsy, making them attractive early-detection biomarker candidates. This project reanalyses publicly available LFQ proteomics data comparing MDA-MB-231 (TNBC) vs. MCF-10A (normal breast epithelial) exosomes.

---

## 👤 Author

**Nguyễn Mậu Tuệ**  
University of Science, Vietnam National University Hanoi (K70)  
*Analysis pipeline developed independently. Documentation assistance provided with the help of AI tools.*
