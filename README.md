# Breast Cancer Exosome Biomarker Discovery

Bioinformatics analysis pipeline for identifying candidate protein biomarkers in breast cancer–derived exosomes using public proteomics data.

---

# Biological Background

Extracellular vesicles (exosomes) play important roles in cancer biology, mediating intercellular communication, metastasis, and immune modulation.

Proteins packaged into tumor-derived exosomes may serve as **non-invasive biomarkers detectable through liquid biopsy**.

This project investigates whether **exosomal proteomic profiles can distinguish breast cancer–associated samples from normal controls**, and prioritizes candidate proteins potentially involved in tumor-related vesicle biology. 

---

# Dataset

Proteomics data were obtained from the PRIDE Archive.

Dataset ID
PXD056161
PRIDE dataset: https://www.ebi.ac.uk/pride/archive/projects/PXD056161

Original study:

Metabolomic and Proteomic Analysis of Breast Cancer Exosomes Reveals Alterations in Purine and Carnitine Metabolism  
DOI: 10.1021/acs.jproteome.4c00795

Experimental system

Cell models:

* MCF-10A (non-tumorigenic breast epithelial cells)
* MDA-MB-231 (tumorigenic breast cancer cells)

Biological material:

* Extracellular vesicles (exosomes)

Proteomics workflow:

* Bottom-up proteomics
* Label-free quantification
* CE-MS with electrospray ionization

Instrumentation:

* Orbitrap Fusion mass spectrometer

Original data processing:

* MaxQuant with Andromeda search engine
* Human UniProt reference proteome

---

# Analysis Workflow

The bioinformatics pipeline consists of several stages:

01. Data loading and preprocessing
02. Quality control and replicate validation
03. Protein filtering (peptide count and sequence coverage)
04. Log2 transformation of LFQ intensities
05. Missing value modeling (MNAR detection)
06. Left-shifted Gaussian imputation
07. Median normalization
08. Differential expression analysis using limma
09. Principal Component Analysis
10. Protein–protein interaction network construction using the STRING database
11. Network topology analysis (degree, betweenness centrality)
12. PCA-based multi-feature integration
13. Composite driver score calculation
14. Functional annotation and literature mining

---

# Key Findings

The prioritization pipeline identified **21 candidate driver proteins** associated with breast cancer exosome proteomic signatures.

Top driver proteins include:

* MVB12A
* CHMP5
* CHMP4B
* CHMP2B
* VPS37B
* VPS4B
* TSG101
* CHMP2A
* VPS28

Many of these proteins belong to the **ESCRT machinery**, which is responsible for membrane remodeling and exosome biogenesis.

These results suggest that **alterations in ESCRT-associated pathways may contribute to tumor-related exosome production in breast cancer cells**.

---

# Software and Tools

Analysis performed in R using:

* data.table
* ggplot2
* limma
* igraph
* ReactomePA
* clusterProfiler

---

# Reproducibility

All scripts required to reproduce the analysis are available in this repository.

Raw mass spectrometry data can be downloaded from the PRIDE repository using dataset identifier **PXD056161**.

---

# Citation

If you use this repository in your work, please cite the project using the information provided in `CITATION.cff`.

---

# Author

Nguyen Mau Tue
Documentation assistance provided with the help of AI tools.
