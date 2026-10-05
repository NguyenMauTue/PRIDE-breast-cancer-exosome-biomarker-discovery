############################################################
# 10 Protein-protein interaction network (STRING)
############################################################
library(igraph); library(here)
source(here::here("R", "Helper", "string_api_helper.R"))

# Gene symbols exported from pathway enrichment
annotated_gene = read.csv(here::here("PXD056161", "results", "tables", "annotated_gene_pool.csv"))
gene_symbols = annotated_gene$hgnc_symbol

network_summary = build_string_network(gene_symbols)
network_summary = network_summary |>
  left_join(annotated_gene, by = c("Symbol" = "hgnc_symbol")) |>
  dplyr::select(Symbol, uniprotswissprot, degree, betweenness)
############################################################
# Save results
############################################################
write.csv(network_summary,
          here::here("PXD056161", "results", "tables", "network_summary.csv"),
          row.names = FALSE)

