############################################################
# Missingness analysis
############################################################

library(reshape2)
library(ggplot2)

############################################################
# Load input data
############################################################

# Load LFQ intensity matrix (log transformed)
log_lfq_matrix <- readRDS(here::here("PXD012162", "data", "log_lfq_matrix.rds"))
metadata <- read.csv(here::here("PXD012162", "data", "sample_metadata.csv"))

############################################################
# Filtering proteins
############################################################

# Identify sample columns for each condition
nor_col = grep("MDA", colnames(log_lfq_matrix))
tum_col = grep("MCF10", colnames(log_lfq_matrix))

# Retain proteins detected in at least two samples
# within either biological condition
filtered_matrix = log_lfq_matrix[
  rowSums(!is.na(log_lfq_matrix[, nor_col])) >= 2 |
    rowSums(!is.na(log_lfq_matrix[, tum_col])) >= 2, ]

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

n5v0_tumor   = sum(tumor_count == 5 & normal_count == 0)
n5v0_normal  = sum(tumor_count == 0 & normal_count == 5)

n5v0_count   = n5v0_tumor + n5v0_normal
n_protein    = nrow(filtered_matrix)

prop_5v0 = n5v0_count / n_protein * 100
prop_5v0

# Approximately 10.91% of proteins exhibited complete 3v0
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

pattern_df$Pattern[normal_count == 5 & tumor_count == 0] =
  "Normal_only_5v0"

pattern_df$Pattern[tumor_count == 5 & normal_count == 0] =
  "Tumor_only_5v0"

############################################################
# Export results
############################################################

exclusive_df = subset(
  pattern_df,
  pattern_df$Pattern == "Normal_only_5v0" |
    pattern_df$Pattern == "Tumor_only_5v0"
)

write.csv(pattern_df,
          here::here("PXD012162", "data", "missing_pattern_summary.csv"),
          row.names = TRUE)

write.csv(exclusive_df,
          here::here("PXD012162", "data", "exclusive_3v0_proteins.csv"),
          row.names = TRUE)


