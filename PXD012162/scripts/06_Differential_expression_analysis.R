############################################################
# 06 Differential expression analysis (limma)
############################################################

library(limma)
library(dplyr)

imputed_matrix = readRDS(here::here("PXD012162", "data", "imputed_matrix.rds"))
log_lfq_matrix = readRDS(here::here("PXD012162", "data", "log_lfq_matrix.rds"))
metadata = readRDS(here::here("PXD012162", "data", "metadata.rds"))

############################################################
# Design matrix
############################################################

group = factor(
  ifelse(grepl("MCF10", colnames(imputed_matrix)), "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)

design = model.matrix(~0 + group)
colnames(design) = levels(group)


############################################################
# Limma model
############################################################

fit = lmFit(imputed_matrix, design)

contrast_matrix =
  makeContrasts(
    Tumor_VS_Normal = Tumor - Normal,
    levels = design
  )

fit2 = contrasts.fit(fit, contrast_matrix)
fit2 = eBayes(fit2)


############################################################
# Extract DE results
############################################################

imputed_result =
  topTable(
    fit2,
    coef = "Tumor_VS_Normal",
    number = Inf,
    sort.by = "P"
  )

imputed_result$X = row.names(imputed_result)
imputed_result <- imputed_result |>
  mutate(UNIPROT = sapply(strsplit(X, ";"), `[`, 1))

write.csv(
  imputed_result,
  here::here("PXD012162", "results", "tables", "differential_expression_imputed.csv")
)


############################################################
# P-value distribution check
############################################################

hist(imputed_result$P.Value,
     breaks = 50,
     col = rgb(0,0,1,0.5),
     main = "Comparing P.Value and FDR",
     xlab = "")

hist(imputed_result$adj.P.Val,
     breaks = 50,
     col = rgb(1,0,0,0.5),
     add = TRUE)

legend("topright",
       legend = c("P-value","FDR"),
       col = c(rgb(0,0,1,0.5),
               rgb(1,0,0,0.5)),
       pch = 16)

dev.off()


############################################################
# Imputation validation (complete case comparison)
############################################################

complete_inx =
  rowSums(is.na(log_lfq_matrix)) == 0

group_cc = factor(
  metadata$Group,
  levels = c("Normal", "Tumor")
)

design_cc = model.matrix(~0 + group_cc)
colnames(design_cc) = levels(group_cc)

fit_cc =
  lmFit(log_lfq_matrix[complete_inx,], design_cc)

contrast_matrix_cc = makeContrasts(
  Tumor_VS_Normal = Tumor - Normal,
  levels = design_cc
)


fit_cc =
  contrasts.fit(fit_cc, contrast_matrix_cc)

fit_cc =
  eBayes(fit_cc)

result_cc =
  topTable(fit_cc, number = Inf)

write.csv(
  result_cc,
  here::here("PXD012162", "results", "tables", "complete_case_DEA.csv")
)

############################################################
# Compare statistics
############################################################

common_genes =
  intersect(
    rownames(imputed_result),
    rownames(result_cc)
  )

length(common_genes)


############################################################
# Correlation of moderated t statistics
############################################################

t_cor =
  cor(
    imputed_result[common_genes,"t"],
    result_cc[common_genes,"t"]
  )

write.csv(
  t_cor,
  here::here("PXD012162", "results", "tables", "imputation_tstat_correlation.csv")
)

