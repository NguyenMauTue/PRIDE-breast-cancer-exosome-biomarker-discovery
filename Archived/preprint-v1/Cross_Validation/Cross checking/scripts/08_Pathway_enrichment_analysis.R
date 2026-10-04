############################################################
# 08 Pathway enrichment analysis (Reactome GSEA)
############################################################

library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

imputed_result =
  read.csv("../results/differential_expression_imputed.csv",
           row.names = 1)


############################################################
# Extract UNIPROT accession
############################################################

protein_ids = rownames(imputed_result)

uniprot_major <-
  str_split(protein_ids, ";") |>
  sapply(`[`, 1) |>
  str_trim()

imputed_result$UNIPROT = uniprot_major


############################################################
# Map UNIPROT -> ENTREZID
############################################################

gene_map =
  bitr(
    unique(imputed_result$UNIPROT),
    fromType = "UNIPROT",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )

rank_df =
  merge(imputed_result,
        gene_map,
        by = "UNIPROT")


############################################################
# Remove duplicated gene mapping
############################################################

rank_df =
  rank_df[order(abs(rank_df$t),
                decreasing = TRUE), ]

rank_df =
  rank_df[!duplicated(rank_df$ENTREZID), ]


############################################################
# Ranked gene list
############################################################

gene_list = rank_df$t
names(gene_list) = as.character(rank_df$ENTREZID)

gene_list =
  sort(gene_list,
       decreasing = TRUE)

gene_list =
  gene_list[!is.na(gene_list)]


############################################################
# Reactome GSEA
############################################################

gsea_reactome =
  gsePathway(
    geneList     = gene_list,
    organism     = "human",
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )

gsea_reactome_res =
  as.data.frame(gsea_reactome)


############################################################
# Pathway filtering
############################################################

filtered_reactome_res =
  gsea_reactome_res[
    abs(gsea_reactome_res$NES) > 1.5 &
    gsea_reactome_res$p.adjust < 0.05,
  ]


############################################################
# Save results
############################################################

write.csv(
  gsea_reactome_res,
  "../results/reactome_gsea_full.csv",
  row.names = FALSE
)

write.csv(
  filtered_reactome_res,
  "../results/reactome_gsea_filtered.csv",
  row.names = FALSE
)

saveRDS(
  gsea_reactome,
  "../results/reactome_gsea_object.rds"
)
