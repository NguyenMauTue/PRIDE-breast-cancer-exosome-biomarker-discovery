############################################################
# 05 Missing value imputation
############################################################

library(ggplot2)
library(ggrepel)
library(patchwork)

set.seed(30032026)
collapsed_matrix = readRDS("../data/collapsed_matrix.rds")
log_lfq_matrix = readRDS("../data/log_lfq_matrix.rds")
metadata = read.csv("../data/sample_metadata.csv")


############################################################
# Split matrix by group
############################################################

nor = collapsed_matrix[, grep("^Normal", colnames(collapsed_matrix))]
tum = collapsed_matrix[, grep("^Tumor", colnames(collapsed_matrix))]

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

complete_inx = rowSums(is.na(collapsed_matrix)) == 0

pca_before = prcomp(
  t(collapsed_matrix[complete_inx,]),
  scale = TRUE
)

pca_after = prcomp(
  t(imputed_matrix),
  scale = TRUE
)


# Helper function để tránh lặp code
make_pca_plot <- function(pca_obj, title) {
  df <- as.data.frame(pca_obj$x)
  percent_var <- round(100 * pca_obj$sdev^2 / sum(pca_obj$sdev^2), 1)
  
  df$sample <- rownames(df)
  df$group  <- ifelse(grepl("Normal", df$sample), "MCF10A (Normal)", "MDA-MB-231 (Tumor)")
  
  ggplot(df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("MCF10A (Normal)" = "#1A9E8F",
                                  "MDA-MB-231 (Tumor)" = "#7B2D8B")) +
    labs(
      title = title,
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Cell line"
    ) +
    theme_bw() +
    theme(
      plot.title   = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
}

p_before <- make_pca_plot(pca_before, "Before imputation")
p_after  <- make_pca_plot(pca_after,  "After imputation")

fig_supp <- p_before + p_after +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("../results/Fig_S1_PCA_before_after_imputation.png",
       fig_supp, width = 10, height = 5, dpi = 300)


############################################################
# PCA structure comparison
############################################################

pc1_correlation =
  cor(pca_before$rotation[,1],
      pca_after$rotation[complete_inx,1], method = "spearman")

write.csv(pc1_correlation,
          "../results/imputation_pca_pc1_correlation.csv")

