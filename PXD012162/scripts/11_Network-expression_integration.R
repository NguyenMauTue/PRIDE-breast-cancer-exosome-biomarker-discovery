############################################################
# 11 Network-expression integration
############################################################

library(biomaRt)
library(dplyr)
library(stringr)

network_summary =
  read.csv(here::here("PXD056161", "results", "tables", "network_summary.csv"))

imputed_result =
  read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"))

#Extracting majority accession (using imputed datasets) 
limma_network_df = imputed_result |>
  inner_join(network_summary, by = c("UNIPROT" =  "uniprotswissprot")) |>
  dplyr::select(UNIPROT, Symbol, logFC, adj.P.Val, degree, betweenness)


write.csv(
  limma_network_df,
  here::here("PXD056161", "results", "tables", "limma_network_table.csv"),
  row.names=FALSE
)
