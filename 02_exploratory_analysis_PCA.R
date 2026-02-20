# Script: 02_exploratory_analysis_PCA.R
# Purpose:
#   - Perform Principal Component Analysis (PCA)
#   - Visualize global proteomic variation
#   - Assess sample clustering by condition
#
# Input:
#   data/processed/proteomics_processed.RData
#
# Output:
#   plots/PCA_plot.png

#Load procesed dât
load("data/processed/proteomics_processed.RData")

#PCA
pca_res = prcomp(t(expMatrixNorm), scale = T)

#Convert to Data.frame
pca_df = as.data.frame(pca_res$x)
pca_df$Group = group

#Variance explained by Scree plot
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
Scree_df = data.frame(PC = colnames(pca_df[, 1:12]),
                      Variance = percentVar)
ggplot(Scree_df, aes(x = 1:nrow(Scree_df), y = Variance)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = 1:nrow(Scree_df),
                     labels = Scree_df$PC) +
  theme_minimal() +
  labs(
    title = "Scree Plot of global proteomic PCA",
    x = "Principal Component",
    y = "Variance Explained (%)"
)
ggsave("plots/Scree_plot.png", scree, width = 6, height = 5, dpi = 300)

#Ploting PCA
library(ggplot2)
plot = ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA reveals global proteomic differences between tumor and normal samples",
       x = paste0("PC1 (", percentVar[1], "% variance)"),
       y = paste0("PC2 (", percentVar[2], "% variance)")) 
ggsave("plots/PCA_plot.png", width = 6, height = 5, dpi = 300)

#Statistical test: whether PC1 differs between tumor and normal samples
t.test(PC1 ~ Group, data = pca_df)
