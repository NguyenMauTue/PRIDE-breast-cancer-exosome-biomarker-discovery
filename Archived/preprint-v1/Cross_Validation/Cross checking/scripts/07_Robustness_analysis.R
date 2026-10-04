############################################################
# 07 Robustness analysis (imputation sensitivity)
############################################################

imputed_result = read.csv(
  "../results/differential_expression_imputed.csv",
  row.names = 1
)

result_cc = read.csv(
  "../results/complete_case_DEA.csv",
  row.names = 1
)


############################################################
# Overlap
############################################################

common_genes =
  intersect(
    rownames(imputed_result),
    rownames(result_cc)
  )

sig_imp = rownames(imputed_result)
sig_cc  = rownames(result_cc)


############################################################
# Robustness metrics
############################################################

robustness_summary = list(
  
  cor_t =
    cor(
      imputed_result[common_genes,"t"],
      result_cc[common_genes,"t"]
    ),
  
  cor_FDR =
    cor(
      -log10(imputed_result[common_genes,"adj.P.Val"]),
      -log10(result_cc[common_genes,"adj.P.Val"])
    ),
  
  n_sig_imp = length(sig_imp),
  n_sig_cc  = length(sig_cc),
  
  overlap =
    length(intersect(sig_imp, sig_cc)),
  
  jaccard =
    length(intersect(sig_imp, sig_cc)) /
    length(union(sig_imp, sig_cc))
  
)


############################################################
# Save robustness metrics
############################################################

saveRDS(
  robustness_summary,
  "../results/limma_robustness_metrics.rds"
)

write.csv(
  data.frame(metric = names(robustness_summary),
             value  = unlist(robustness_summary)),
  "../results/limma_robustness_metrics.csv",
  row.names = FALSE
)
  
