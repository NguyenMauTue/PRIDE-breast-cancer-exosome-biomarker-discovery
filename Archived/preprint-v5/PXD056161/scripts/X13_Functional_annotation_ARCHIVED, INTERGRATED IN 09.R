############################################################
# 13 Functional annotation
############################################################

library(biomaRt)
library(dplyr)

#Load driver scores
Driver_table = read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_annotated.csv"))

#List UniProt
uniprot_list = Driver_table$UNIPROT

#Query annotation
annotations_all = read.csv(here::here("PXD056161", "data", "annotation_raw.csv"))
annotations = annotations_all %>%
  filter(uniprotswissprot %in% uniprot_list)

#Checking annotation
missing = setdiff(uniprot_list, annotations$uniprotswissprot)
if (length(missing) > 0) {
  warning("Các UniProt ID sau chưa có trong annotation_raw.csv, cần chạy lại 00_fetch_annotation_table.R: ",
          paste(missing, collapse = ", "))
}
#Group GO terms
annotations_grouped = annotations %>%
  group_by(uniprotswissprot) %>%
  summarise(
    Symbol_biomat = paste(unique(external_gene_name), collapse="; "),
    Protein_Description = paste(unique(description), collapse="; "),
    GO_BP = paste(name_1006[namespace_1003=="biological_process"], collapse="; "),
    GO_CC = paste(name_1006[namespace_1003=="cellular_component"], collapse="; "),
    GO_MF = paste(name_1006[namespace_1003=="molecular_function"], collapse="; "),
    # Thêm GO ID columns
    GO_BP_IDs = paste(go_id[namespace_1003=="biological_process"], collapse="; "),
    GO_CC_IDs = paste(go_id[namespace_1003=="cellular_component"], collapse="; "),
    GO_MF_IDs = paste(go_id[namespace_1003=="molecular_function"], collapse="; ")
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
  here::here("PXD056161", "results", "tables", "BiomarkerCandidates_with_annotations.csv"),
  row.names = FALSE
)

