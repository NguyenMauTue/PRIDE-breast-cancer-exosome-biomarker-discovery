############################################################
# 09 Pathway thematic grouping and network preparation
############################################################

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

gsea_reactome_res =
  read.csv("../results/reactome_gsea_full.csv")


gene_lists =
  strsplit(
    gsea_reactome_res$core_enrichment,
    "/"
  )

genes =
  unique(
    unlist(gene_lists)
  )

############################################################
# Convert ENTREZ -> SYMBOL
############################################################

symbols =
  bitr(
    genes,
    fromType = "ENTREZID",
    toType   = "SYMBOL",
    OrgDb    = org.Hs.eg.db
  )

gene_symbols =
  symbols$SYMBOL

############################################################
# Export gene list for STRING network (unbiased)
############################################################

write.table(
  gene_symbols,
  "../results/genes_for_STRING.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
