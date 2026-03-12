############################################################
# 06 Differential expression analysis (limma)
############################################################

library(limma)

imputed_matrix = readRDS("../data/imputed_matrix.rds")
collapsed_matrix = readRDS("../data/collapsed_matrix.rds")


############################################################
# Design matrix
############################################################

group = factor(
  c(rep("Normal",3), rep("Tumor",3)),
  levels = c("Normal","Tumor")
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

write.csv(
  imputed_result,
  "../results/differential_expression_imputed.csv"
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
  rowSums(is.na(collapsed_matrix)) == 0

fit_cc =
  lmFit(collapsed_matrix[complete_inx,], design)

fit_cc =
  contrasts.fit(fit_cc, contrast_matrix)

fit_cc =
  eBayes(fit_cc)

result_cc =
  topTable(fit_cc, number = Inf)

write.csv(
  result_cc,
  "../results/complete_case_DEA.csv"
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
  "../results/imputation_tstat_correlation.csv"
)
