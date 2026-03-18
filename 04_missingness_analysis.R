############################################################
# Missingness analysis
############################################################

library(reshape2)
library(ggplot2)

############################################################
# Load input data
############################################################

# Load LFQ intensity matrix (log transformed)
log_lfq_matrix <- readRDS("data/log_lfq_matrix.rds")
metadata <- read.csv("data/sample_metadata.csv")

############################################################
# Collapsing
############################################################
unique_bio = unique(metadata$BioRep_full)
collapsed_matrix = sapply(unique_bio, function(bio) {
  cols = metadata$Sample[metadata$BioRep_full == bio]
  rowMeans(log_lfq_matrix[, cols], na.rm = TRUE)
})
colnames(collapsed_matrix) = unique_bio

############################################################
# Filtering proteins
############################################################

# Identify sample columns for each condition
nor_col = grep("^Normal", colnames(collapsed_matrix))
tum_col = grep("^Tumor", colnames(collapsed_matrix))

saveRDS(collapsed_matrix,
        "../data/collapsed_matrix.rds")
# Retain proteins detected in at least two samples
# within either biological condition
filtered_matrix = collapsed_matrix[
  rowSums(!is.na(collapsed_matrix[, nor_col])) >= 2 |
    rowSums(!is.na(collapsed_matrix[, tum_col])) >= 2, ]

############################################################
# Evaluate missingness mechanism
############################################################

# Convert LFQ matrix to long format
long_df <- melt(log_lfq_matrix)
colnames(long_df) <- c("Protein", "Sample", "Intensity")

# Detection status
long_df$Detected <- !is.na(long_df$Intensity)

# Mean abundance per protein
protein_mean <- rowMeans(log_lfq_matrix, na.rm = TRUE)
long_df$MeanIntensity <- protein_mean[long_df$Protein]

############################################################
# Logistic regression: detection probability vs abundance
############################################################

fit_det <- glm(
  Detected ~ MeanIntensity,
  data = long_df,
  family = binomial
)

summary(fit_det)

############################################################
# Visualization of detection bias
############################################################

# Convert logical detection to numeric
long_df$Detected_Number = 0
long_df$Detected_Number[long_df$Detected == TRUE] = 1

ggplot(long_df, aes(x = MeanIntensity, y = Detected_Number)) +
  geom_jitter(height = 0.05, width = 0.1, alpha = 0.1) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE,
              color = "red",
              size = 1.2) +
  theme_minimal()

# Interpretation example:
# Logistic regression analysis demonstrated a strong positive
# association between protein abundance and detection
# probability (β = 0.43, p < 2e-16), indicating abundance-
# dependent missingness consistent with a limit-of-detection
# mechanism.

############################################################
# Quantification of extreme missingness patterns (3v0)
############################################################

normal_count = rowSums(!is.na(filtered_matrix[, nor_col]))
tumor_count  = rowSums(!is.na(filtered_matrix[, tum_col]))

n3v0_tumor   = sum(tumor_count == 3 & normal_count == 0)
n3v0_normal  = sum(tumor_count == 0 & normal_count == 3)

n3v0_count   = n3v0_tumor + n3v0_normal
n_protein    = nrow(filtered_matrix)

prop_3v0 = n3v0_count / n_protein * 100
prop_3v0

# Approximately 17.7% of proteins exhibited complete 3v0
# missingness patterns. Given the previously demonstrated
# abundance-dependent detection mechanism, these cases
# likely reflect censoring effects rather than true
# biological absence.

############################################################
# Missingness pattern classification
############################################################

pattern_df = data.frame(
  row.names = rownames(filtered_matrix),
  Normal_count = normal_count,
  Tumor_count  = tumor_count
)

pattern_df$Pattern = "Other"

pattern_df$Pattern[normal_count == 3 & tumor_count == 0] =
  "Normal_only_3v0"

pattern_df$Pattern[tumor_count == 3 & normal_count == 0] =
  "Tumor_only_3v0"

############################################################
# Export results
############################################################

exclusive_df = subset(
  pattern_df,
  pattern_df$Pattern == "Normal_only_3v0" |
    pattern_df$Pattern == "Tumor_only_3v0"
)

write.csv(pattern_df,
          "../data/missing_pattern_summary.csv",
          row.names = TRUE)

write.csv(exclusive_df,
          "../data/exclusive_3v0_proteins.csv",
          row.names = TRUE)

