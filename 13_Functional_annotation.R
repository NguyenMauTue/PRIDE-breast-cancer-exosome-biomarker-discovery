############################################################
# 13 Functional annotation
############################################################

library(biomaRt)
library(dplyr)

#Load driver scores
Driver_table = read.csv("../results/CDS_candidates_annotated.csv")

#Chọn mart cho human
mart = useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

#List UniProt
uniprot_list = Driver_table$UNIPROT

#Query annotation
annotations = getBM(
  attributes = c(
    "uniprotswissprot",
    "external_gene_name",
    "description",
    "go_id",
    "name_1006",
    "namespace_1003"
  ),
  filters = "uniprotswissprot",
  values = uniprot_list,
  mart = mart
)

#Group GO terms
annotations_grouped = annotations %>%
  group_by(uniprotswissprot) %>%
  summarise(
    Symbol_biomat = paste(unique(external_gene_name), collapse="; "),
    Protein_Description = paste(unique(description), collapse="; "),
    GO_BP = paste(name_1006[namespace_1003=="biological_process"], collapse="; "),
    GO_CC = paste(name_1006[namespace_1003=="cellular_component"], collapse="; "),
    GO_MF = paste(name_1006[namespace_1003=="molecular_function"], collapse="; ")
  ) %>%
  ungroup()

#Merge annotation
BiomarkerCandidates_annotated =
  merge(
    Driver_table,
    annotations_grouped,
    by.x="UNIPROT",
    by.y="uniprotswissprot",
    all.x=TRUE
  )

#Export
write.csv(
  BiomarkerCandidates_annotated,
  "../results/BiomarkerCandidates_with_annotations.csv",
  row.names = FALSE
)
