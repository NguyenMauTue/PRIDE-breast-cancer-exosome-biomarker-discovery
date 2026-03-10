############################################################
# 04 Missingness mechanism analysis
############################################################

library(tidyverse)

log_lfq_matrix = readRDS("../data/log_lfq_matrix.rds")
metadata = read.csv("../data/sample_metadata.csv")

############################################################
# 1 Convert matrix to long format
############################################################

lfq_long =
  as.data.frame(log_lfq_matrix) %>%
  rownames_to_column("Protein") %>%
  pivot_longer(
    -Protein,
    names_to = "Sample",
    values_to = "Intensity"
  )

lfq_long =
  lfq_long %>%
  left_join(metadata, by = c("Sample" = "Sample"))

############################################################
# 2 Detection variable
############################################################

lfq_long$Detected =
  ifelse(is.na(lfq_long$Intensity), 0, 1)

############################################################
# 3 Detection probability vs intensity
############################################################

det_plot =
  ggplot(lfq_long,
         aes(x = Intensity, y = Detected)) +
  geom_jitter(height = 0.05, alpha = 0.2) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              color = "red") +
  theme_minimal() +
  labs(
    title = "Detection probability vs intensity",
    x = "log2 LFQ intensity",
    y = "Detection probability"
  )

ggsave(
  "../results/detection_probability_curve.png",
  det_plot,
  width = 6,
  height = 5
)

############################################################
# 4 Logistic regression model
############################################################

det_model =
  glm(
    Detected ~ Intensity,
    data = lfq_long,
    family = binomial
  )

capture.output(
  summary(det_model),
  file = "../results/missingness_logistic_model.txt"
)

############################################################
# 5 Missingness per protein
############################################################

missing_summary =
  lfq_long %>%
  group_by(Protein) %>%
  summarise(
    detected = sum(Detected),
    missing = sum(Detected == 0),
    missing_fraction = mean(Detected == 0)
  )

write.csv(
  missing_summary,
  "../results/protein_missingness_summary.csv",
  row.names = FALSE
)

############################################################
# 6 Detect 3v0 pattern (complete separation)
############################################################

presence_matrix =
  !is.na(log_lfq_matrix)

group_vector =
  metadata$Group

protein_pattern =
  data.frame(
    Protein = rownames(log_lfq_matrix),
    Tumor_presence =
      rowSums(presence_matrix[, group_vector == "Tumor"]),
    Normal_presence =
      rowSums(presence_matrix[, group_vector == "Normal"])
  )

three_vs_zero =
  protein_pattern %>%
  filter(
    (Tumor_presence > 0 & Normal_presence == 0) |
    (Tumor_presence == 0 & Normal_presence > 0)
  )

write.csv(
  three_vs_zero,
  "../results/three_vs_zero_proteins.csv",
  row.names = FALSE
)

############################################################
# 7 Summary statistics
############################################################

summary_stats =
  data.frame(
    total_proteins = nrow(log_lfq_matrix),
    proteins_with_missing =
      sum(rowSums(is.na(log_lfq_matrix)) > 0),
    proteins_3v0 = nrow(three_vs_zero)
  )

write.csv(
  summary_stats,
  "../results/missingness_summary_stats.csv",
  row.names = FALSE
)
