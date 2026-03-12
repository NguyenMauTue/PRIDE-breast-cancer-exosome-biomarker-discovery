############################################################
# 05 Missing value imputation
############################################################

library(ggplot2)
library(ggrepel)

log_lfq_matrix = readRDS("../data/log_lfq_matrix.rds")
metadata = read.csv("../data/sample_metadata.csv")


############################################################
# Split matrix by group
############################################################

nor = log_lfq_matrix[,1:3]
tum = log_lfq_matrix[,4:6]


############################################################
# Left-censored Gaussian imputation
############################################################

impute_left_gaussian = function(mat, shift = 1.8, width = 0.3){

  observed = mat[!is.na(mat)]

  mu = mean(observed)
  sigma = sd(observed)

  impute_values = rnorm(
    sum(is.na(mat)),
    mean = mu - shift * sigma,
    sd = width * sigma
  )

  mat[is.na(mat)] = impute_values

  return(mat)
}


############################################################
# Perform imputation
############################################################

imputed_nor = impute_left_gaussian(nor)
imputed_tum = impute_left_gaussian(tum)

imputed_matrix = cbind(imputed_nor, imputed_tum)

saveRDS(imputed_matrix,
        "../data/imputed_matrix.rds")


############################################################
# Distribution check
############################################################

hist(as.vector(log_lfq_matrix),
     breaks = 50,
     col = rgb(0,0,1,0.5),
     xlim = range(c(as.vector(imputed_matrix),
                    as.vector(log_lfq_matrix)),
                  na.rm = TRUE),
     main = "Distribution of original vs imputed values",
     xlab = "Log2 LFQ intensity")

hist(as.vector(imputed_matrix),
     breaks = 50,
     col = rgb(1,0,0,0.5),
     add = TRUE)

legend("topright",
       legend = c("Original", "Imputed"),
       fill = c(rgb(0,0,1,0.5),
                rgb(1,0,0,0.5)))

dev.off()


############################################################
# PCA stability test
############################################################

complete_inx = rowSums(is.na(log_lfq_matrix)) == 0

pca_before = prcomp(
  t(log_lfq_matrix[complete_inx,]),
  scale = TRUE
)

pca_after = prcomp(
  t(imputed_matrix),
  scale = TRUE
)


############################################################
# PCA before imputation
############################################################

percent_var = round(
  100 * (pca_before$sdev^2 / sum(pca_before$sdev^2)),
  1
)

p1 = ggplot(as.data.frame(pca_before$x),
       aes(x = PC1, y = PC2)) +
  geom_point(size = 4) +
  theme_minimal() +
  geom_text_repel(
    aes(label = rownames(as.data.frame(pca_before$x))),
    size = 3
  ) +
  labs(
    title = "PCA before imputation",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  )

ggsave("../results/pca_before_imputation.png",
       p1,
       width = 6,
       height = 5)


############################################################
# PCA after imputation
############################################################

percent_var = round(
  100 * (pca_after$sdev^2 / sum(pca_after$sdev^2)),
  1
)

p2 = ggplot(as.data.frame(pca_after$x),
       aes(x = PC1, y = PC2)) +
  geom_point(size = 4) +
  theme_minimal() +
  geom_text_repel(
    aes(label = rownames(as.data.frame(pca_after$x))),
    size = 3
  ) +
  labs(
    title = "PCA after imputation",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  )

ggsave("../results/pca_after_imputation.png",
       p2,
       width = 6,
       height = 5)


############################################################
# PCA structure comparison
############################################################

pc1_correlation =
  cor(pca_before$rotation[,1],
      pca_after$rotation[complete_inx,1])

write.csv(pc1_correlation,
          "../results/imputation_pca_pc1_correlation.csv")
