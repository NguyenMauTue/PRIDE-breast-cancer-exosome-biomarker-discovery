############################################################
# 11 Network-expression integration
############################################################

library(biomaRt)
library(dplyr)
library(stringr)

network_summary =
  read.csv("../results/network_summary.csv")

imputed_result =
  read.csv("../results/differential_expression_imputed.csv")

#Extracting majority accession (using imputed datasets) 
protein_ids = imputed_result$X

uniprot_major <- str_split(protein_ids, ";") |>
  sapply(`[`, 1) |>
  str_trim()
imputed_result$UNIPROT = uniprot_major

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

mapping =
  getBM(
    attributes=c(
      "hgnc_symbol",
      "uniprotswissprot"
    ),
    filters="hgnc_symbol",
    values=network_summary$Symbol,
    mart=mart
  )

colnames(mapping)=
  c("SYMBOL","UNIPROT")

mapping =
  mapping[mapping$UNIPROT!="", ]

mapping =
  mapping[!duplicated(mapping$SYMBOL), ]


network_summary2 =
  merge(
    network_summary,
    mapping,
    by.x="Symbol",
    by.y="SYMBOL"
  )


limmanetwork_df =
  merge(
    network_summary2,
    imputed_result,
    by="UNIPROT",
  )


write.csv(
  limmanetwork_df,
  "../results/limma_network_table.csv",
  row.names=FALSE
)
