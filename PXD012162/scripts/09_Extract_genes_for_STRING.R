############################################################
# 09 Extract genes for STRING
############################################################

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)

gsea_reactome_res =
  read.csv(here::here("PXD056161", "results", "tables", "reactome_gsea_filtered.csv"))
gene_lists =
  strsplit(
    gsea_reactome_res$core_enrichment,
    "/"
  )
genes =
  unique(
    unlist(gene_lists)
  )
mart = NULL
mirrors = c("asia", "useast", "uswest", "www")
for (m in mirrors) {
  tryCatch({
    mart = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = m)
    message("Connected via mirror: ", m)
    break
  }, error = function(e) message("Mirror failed: ", m))
}
if (is.null(mart)) stop("All Ensembl mirrors failed.")
results <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "uniprotswissprot"),
  filters = "entrezgene_id",
  values = genes,
  mart = mart
)
results <- results %>%
  filter(uniprotswissprot != "" & !is.na(uniprotswissprot)) %>%
  distinct()
#Query annotation
annotations_all = read.csv(here::here("PXD056161", "data", "annotation_raw.csv"))
annotations = annotations_all %>%
  dplyr::filter(uniprotswissprot %in% results$uniprotswissprot)

#Group GO terms
annotations_grouped = annotations %>%
  group_by(uniprotswissprot) %>%
  summarise(
    Symbol_biomat = paste(unique(external_gene_name), collapse="; "),
    Protein_Description = paste(unique(description), collapse="; "),
    GO_BP = paste(name_1006[namespace_1003=="biological_process"], collapse="; "),
    GO_CC = paste(name_1006[namespace_1003=="cellular_component"], collapse="; "),
    GO_MF = paste(name_1006[namespace_1003=="molecular_function"], collapse="; "),
    # add GO ID columns
    GO_BP_IDs = paste(go_id[namespace_1003=="biological_process"], collapse="; "),
    GO_CC_IDs = paste(go_id[namespace_1003=="cellular_component"], collapse="; "),
    GO_MF_IDs = paste(go_id[namespace_1003=="molecular_function"], collapse="; ")
  ) %>%
  ungroup()

#Merge annotation
Symbol_annotated =
  merge(
    results,
    annotations_grouped,
    by.x="uniprotswissprot",
    by.y="uniprotswissprot",
    all.x=TRUE
  )

# Contaminant filter
contaminant_pattern <- "histone|keratin|actin|tubulin"

# Manual blacklist
manual_blacklist <- c("FMNL1", "H2AX", "H4C6", "H3C1", "H3-3B", "H2AZ2")

Symbol_annotated <- Symbol_annotated |>
  filter(!grepl(contaminant_pattern, Protein_Description, ignore.case = TRUE)) |>
  filter(!grepl(contaminant_pattern, as.character(Symbol_biomat), ignore.case = TRUE)) |>
  filter(!Symbol_biomat %in% manual_blacklist) 
Symbol_annotated <- na.omit(Symbol_annotated)

############################################################
# Export gene list for STRING network 
############################################################
write.csv(
  Symbol_annotated,
  here::here("PXD056161", "results", "tables", "annotated_gene_pool.csv"),
  row.names = FALSE,
  quote = TRUE
)