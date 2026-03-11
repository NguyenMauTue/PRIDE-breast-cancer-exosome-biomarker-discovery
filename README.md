# Breast-cancer-exosome-biomarker-discovery
## Research question
Can heparin-binding proteins distinguish cancer cells from normal cells based on proteomics profiles?

## Dataset
PRIDE dataset: PXD056161
Exosome proteomics data from breast tissue of breast cancer patients and healthy individuals.
Technique: CE-MS (Capillary Electrophoresis Mass Spectrometry)

## Methods
- Data filtering (peptide count, sequence coverage)
- Log2 transformation
- Missing value imputation (left-shifted Gaussian)
- Median normalization
- Principal Component Analysis (PCA)
- Statistical testing (t-test)
- Visualization (ggplot2)

## Results
PCA analysis reveals clear separation between tumor and normal samples.
Statistical testing confirms significant association between proteomic profile and disease state.

## Tools
- R
- data.table
- ggplot2

## Learning purpose
This project was developed as part of my learning process in proteomics data analysis.
I explored publicly available datasets and implemented statistical workflows to understand how protein expression differs between biological conditions.
## Author
Nguyen Mau Tue
