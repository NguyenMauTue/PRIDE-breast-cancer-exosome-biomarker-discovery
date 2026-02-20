# Script: 01_download_and_preprocess_data.R
# Project: PRIDE proteomics differential expression analysis
# Description:
#   - Download proteomics data from PRIDE (PXD056161)
#   - Filter low-quality proteins
#   - Handle missing values
#   - Perform log2 transformation
#   - Impute missing values (left-shifted Gaussian)
#   - Normalize expression data
# Output:
#   data/processed/proteomics_processed.RData


#Data extracting
library(data.table)
proteinGroup = fread("https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD056161/proteinGroups.txt")

#Filtering
proteinGroup = proteinGroup[proteinGroup$Peptides >= 2]
proteinGroup = proteinGroup[proteinGroup$`Sequence coverage [%]` >= 5]

#Extracting LFQ intensity matrix
lfqData = grep("LFQ intensity", colnames(proteinGroup), value= T)
expMatrix = as.matrix(proteinGroup[, ..lfqData])
rownames(expMatrix) = make.unique(proteinGroup$`Protein IDs`)
expMatrix[expMatrix == "" | expMatrix == 0] = NA

#Log2 transformation
logExpMatrix = log2(expMatrix)
expMatrixClean = logExpMatrix[rowSums(!is.na(logExpMatrix)) >= 6, ]
group = factor(c(rep("Normal", 6), rep("Tumor", 6)), levels = c("Tumor", "Normal"))
metadata = data.frame(sample = sub("LFQ intensity ", "", colnames(expMatrixClean)),
                      group = group)

#Imputation NA (left-shifted gaussian)
flatexpMatrixClean = as.vector(as.matrix(expMatrixClean))
mu = mean(flatexpMatrixClean, na.rm = T)
sigma = sd(flatexpMatrixClean, na.rm = T)
expMatrixClean = as.matrix(expMatrixClean)
expMatrixClean[is.na(expMatrixClean)] = rnorm(sum(is.na(expMatrixClean)), mean = mu - 1.8 * sigma, sd = 0.3 * sigma)

#Normalisation
medians = apply(expMatrixClean, 2, median)
expMatrixNorm = sweep(expMatrixClean, 2, medians, "-")

#Output
save(expMatrixNorm, metadata, file = "data/processed/proteomics_processed.RData")
