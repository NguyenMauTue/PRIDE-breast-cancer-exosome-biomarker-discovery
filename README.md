# PRIDE Breast Cancer Exosome Biomarker Discovery

### A reproducible pipeline implementing the Composite Driver Score (CDS) framework for exosomal protein biomarker prioritisation in triple-negative breast cancer

---

> **Preprint:** [https://doi.org/10.64898/2026.05.14.725271](https://doi.org/10.64898/2026.05.14.725271)
> 
> **Manuscript PDF:** [bioRxiv PDF](https://www.biorxiv.org/content/10.64898/2026.05.14.725271v1.full.pdf)
> 
> **Data:** [PRIDE PXD056161](https://www.ebi.ac.uk/pride/archive/projects/PXD056161) · [PRIDE PXD012162](https://www.ebi.ac.uk/pride/archive/projects/PXD012162)

---

## 🔑 Key Idea

Most exosomal proteomics studies rank candidates by fold-change or statistical significance alone. This pipeline implements a **Composite Driver Score (CDS)** that integrates expression magnitude with network topology using an **Analytic Hierarchy Process (AHP)**, enabling prioritisation of proteins that are both differentially abundant and biologically embedded in functionally coherent networks.

Rather than filtering by arbitrary FDR thresholds — inappropriate given n=3 per condition — CDS weights each criterion by its relative contribution to clinical detectability and biological plausibility. The result is a ranked, robustness-tested candidate list designed to generate experimentally testable hypotheses.

---

## 📦 Repository Structure

```
├── scripts/              # Main analysis pipeline (01–16, run sequentially)
├── data/                 # Input files (MaxQuant output, STRING network TSV)
├── results/              # Final tables, DE results, GSEA objects and all figures
├── Cross checking/       # Cross-dataset validation pipeline (PXD012162)
│   ├── scripts/          # Mirrored pipeline for validation dataset
│   └── results/          # Cross-validation outputs
└── Master run_all.R      # Single entry point — runs both pipelines end-to-end
```

## ⚙️ How to Run

```r
# 1. Clone the repository
# git clone https://github.com/NguyenMauTue/PRIDE-breast-cancer-exosome-biomarker-discovery

# 2. Install dependencies (R >= 4.4)
install.packages(c("dplyr", "ggplot2", "ggrepel", "openxlsx",
                   "patchwork", "circlize", "rentrez", "xtable"))
BiocManager::install(c("limma", "clusterProfiler", "ReactomePA",
                       "biomaRt", "org.Hs.eg.db", "reactome.db",
                       "GO.db", "igraph"))

# 3. Download raw data from PRIDE (PXD056161) — script 01 loads automatically via FTP

# 4. Run full pipeline (main + cross-validation)
source("Master run_all.R")

# or run main pipeline only
source("scripts/run_all.R")
```

> A `sessionInfo()` export is saved to `results/sessionInfo.txt` after each full run.  
> An `renv.lock` file is provided for full environment reproducibility.

## 📊 Outputs

| File | Description |
|---|---|
| `results/Module_tables.xlsx` | Full ranked candidate table with CDS, biological module, robustness, PubMed hits |
| `results/Supplementary_Table1.tex` | LaTeX-ready longtable for manuscript Supplementary Material |
| `results/limma_network_table.csv` | DE results merged with STRING network topology |
| `results/agrn_pathway_dir_summary.csv` | AGRN satellite protein directionality by Reactome pathway |
| `results/volcano.png` | Volcano plot — differential expression |
| `results/DriverScore_landscape.png` | CDS landscape plot |
| `results/PPI_network.png` | STRING network with DE integration |
| `results/GO_enrichment.png` | GO enrichment dot plot |
| `results/16_AGRN_satellite_directionality.png` | AGRN pathway co-protein concordance |


## ⚗️ Methods Summary

### Biological Module Classification

Proteins assigned to biological modules via GO term hierarchy using `GO.db` (Bioconductor). Module membership determined by matching annotated GO term IDs against programmatically derived offspring sets of predefined root terms:

| Module | Root Terms |
|---|---|
| ECM & Cell Adhesion | GO:0031012, GO:0030198, GO:0007160 |
| Motility & Signaling | GO:0016477, GO:0000165, GO:0007265, GO:0035023 |
| Vesicle Trafficking | GO:0016192, GO:0036258 |

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
