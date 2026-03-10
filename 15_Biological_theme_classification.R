############################################################
# 15 Biological theme classification
############################################################

library(dplyr)
library(rentrez)

BiomarkerCandidates_annotated =
  read.csv("../results/BiomarkerCandidates_with_annotations.csv")

#Theme classification
BiomarkerCandidates_themed =
  BiomarkerCandidates_annotated %>%
  mutate(

    Theme_Vesicle_Trafficking =
      grepl("vesicle|endosomal|ESCRT|exosome",
            paste(GO_BP,GO_CC),
            ignore.case=TRUE),

    Theme_Motility_Signaling =
      grepl("RHO|motility|migration|cytoskeleton",
            GO_BP,
            ignore.case=TRUE),

    Theme_Cell_Adhesion =
      grepl("integrin|adhesion|junction|ECM",
            paste(GO_BP,GO_CC),
            ignore.case=TRUE),

    Theme_Immune_Evasion =
      grepl("interferon|cytokine|antigen|MHC|PD|CTLA",
            GO_BP,
            ignore.case=TRUE)
  )

#Theme count
BiomarkerCandidates_themed$Num_Themes_Match =
  rowSums(
    BiomarkerCandidates_themed[
      ,grep("Theme_",colnames(BiomarkerCandidates_themed))
    ],
    na.rm=TRUE
  )

#Localization
BiomarkerCandidates_themed$Localization = case_when(

  grepl("extracellular",
        BiomarkerCandidates_themed$GO_CC,
        ignore.case=TRUE) ~ "extracellular",

  grepl("plasma membrane",
        BiomarkerCandidates_themed$GO_CC,
        ignore.case=TRUE) ~ "membrane",

  grepl("vesicle|endosome|multivesicular|exosome",
        BiomarkerCandidates_themed$GO_CC,
        ignore.case=TRUE) ~ "vesicle/exosome",

  grepl("cytosol|cytoplasm|nucleus",
        BiomarkerCandidates_themed$GO_CC,
        ignore.case=TRUE) ~ "intracellular",

  TRUE ~ "other"
)

#PubMed mining
check_pubmed = function(gene){

  query = paste(gene,"cancer")

  res = entrez_search(
    db="pubmed",
    term=query
  )

  return(res$count)
}

BiomarkerCandidates_themed$PubMed_hits =
  sapply(
    BiomarkerCandidates_themed$Symbol,
    check_pubmed
  )

#Save
write.csv(
  BiomarkerCandidates_themed,
  "../results/BiomarkerCandidates_themed.csv",
  row.names = FALSE
)
