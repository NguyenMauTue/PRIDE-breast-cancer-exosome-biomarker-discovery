############################################################
# 03 Technical replicate QC
############################################################

library(ggplot2)
library(ggrepel)

log_lfq_matrix = readRDS("../data/log_lfq_matrix.rds")
metadata = read.csv("../data/sample_metadata.csv")

# Correlation matrix
cor_matrix = cor(log_lfq_matrix, use = "pairwise.complete.obs")

write.csv(
  cor_matrix,
  "../results/technical_correlation_matrix.csv"
)

# PCA (only complete proteins)
complete_idx = rowSums(is.na(log_lfq_matrix)) == 0

pca_matrix = log_lfq_matrix[complete_idx, ]

pca = prcomp(t(pca_matrix), scale = TRUE)

pca_df = as.data.frame(pca$x)

pca_df$Group = metadata$Group

# Variance explained
percent_var = round(
  100 * (pca$sdev^2 / sum(pca$sdev^2)),
  2
)

# PCA plot
p = ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal() +
  geom_text_repel(
    aes(label = rownames(pca_df)),
    size = 3
  ) +
  labs(
    title = "PCA quality testing",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)")
  )

ggsave(
  "../results/pca_quality_control.png",
  plot = p,
  width = 6,
  height = 5
)
