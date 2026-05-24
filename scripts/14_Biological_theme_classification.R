############################################################
# 14 Biological theme classification
############################################################
library(dplyr)
library(rentrez)
library(GO.db)
library(AnnotationDbi)

BiomarkerCandidates_annotated <-
  read.csv("../results/BiomarkerCandidates_with_annotations.csv")

# Define root GO terms and its offspring
get_offspring_bp <- function(go_id) {
  offspring <- tryCatch(
    as.character(GOBPOFFSPRING[[go_id]]),
    error = function(e) character(0)
  )
  return(c(go_id, offspring))
}

get_offspring_cc <- function(go_id) {
  offspring <- tryCatch(
    as.character(GOCCOFFSPRING[[go_id]]),
    error = function(e) character(0)
  )
  return(c(go_id, offspring))
}

# ECM root terms
ecm_cc_terms <- get_offspring_cc("GO:0031012")  # extracellular matrix
ecm_bp_terms <- c(
  get_offspring_bp("GO:0030198"),  # ECM organization
  get_offspring_bp("GO:0007160"),  # cell-matrix adhesion
  get_offspring_bp("GO:0043062")   # extracellular structure organization
)

# Motility/Signaling root terms
ms_bp_terms <- c(
  get_offspring_bp("GO:0016477"),  # cell migration
  get_offspring_bp("GO:0000165"),  # MAPK cascade
  get_offspring_bp("GO:0007265"),  # Ras protein signal transduction
  get_offspring_bp("GO:0035023")   # regulation of Rho protein signal transduction
)
ms_cc_terms <- c(
  get_offspring_cc("GO:0015629"),  # actin cytoskeleton
  get_offspring_cc("GO:0001725")   # stress fiber
)

# Vesicle Trafficking root terms
vt_bp_terms <- c(
  get_offspring_bp("GO:0016192"),  # vesicle-mediated transport
  get_offspring_bp("GO:0036258")   # multivesicular body assembly
)
vt_cc_terms <- c(
  get_offspring_cc("GO:0005765"),  # lysosomal membrane
  get_offspring_cc("GO:0070971")   # extracellular vesicle biogenesis
)

# Helper function: check if any GO IDs match term list
has_term <- function(id_string, term_list) {
  if (is.na(id_string) || id_string == "") return(FALSE)
  ids <- trimws(strsplit(id_string, ";")[[1]])
  any(ids %in% term_list)
}

# Theme classification
BiomarkerCandidates_themed <-
  BiomarkerCandidates_annotated %>%
  rowwise() %>%
  mutate(
    Theme_Cell_Adhesion =
      has_term(GO_CC_IDs, ecm_cc_terms) |
      has_term(GO_BP_IDs, ecm_bp_terms),
    
    Theme_Motility_Signaling =
      has_term(GO_BP_IDs, ms_bp_terms) |
      has_term(GO_CC_IDs, ms_cc_terms),
    
    Theme_Vesicle_Trafficking =
      has_term(GO_BP_IDs, vt_bp_terms) |
      has_term(GO_CC_IDs, vt_cc_terms)
  ) %>%
  ungroup()

# Theme count
BiomarkerCandidates_themed$Num_Themes_Match <-
  rowSums(
    BiomarkerCandidates_themed[
      , grep("Theme_", colnames(BiomarkerCandidates_themed))
    ],
    na.rm = TRUE
  )

# Biological Module — multi-tag, comma-separated
BiomarkerCandidates_themed$Biological_Module <-
  apply(
    BiomarkerCandidates_themed[, grep("Theme_", colnames(BiomarkerCandidates_themed))],
    1,
    function(x) {
      theme_map <- c(
        Theme_Cell_Adhesion      = "ECM",
        Theme_Motility_Signaling = "MS",
        Theme_Vesicle_Trafficking = "VT"
      )
      tags <- theme_map[as.logical(x)]
      if (length(tags) == 0) "Other" else paste(tags, collapse = ",")
    }
  )

# PubMed mining
check_pubmed <- function(gene) {
  query <- paste(gene, "exosome")
  res <- entrez_search(db = "pubmed", term = query)
  return(res$count)
}

BiomarkerCandidates_themed$PubMed_hits <-
  sapply(
    BiomarkerCandidates_themed$Symbol,
    check_pubmed
  )

# Contaminant filter
contaminant_pattern <- "histone|keratin|actin|tubulin"

BiomarkerCandidates_themed[
  grepl(contaminant_pattern, BiomarkerCandidates_themed$Protein_Description, ignore.case = TRUE) |
    grepl(contaminant_pattern, as.character(BiomarkerCandidates_themed$Symbol), ignore.case = TRUE),
  "Symbol"
]

# Manual blacklist
manual_blacklist <- c("FMNL1", "H2AX", "H4C6", "H3C1", "H3-3B", "H2AZ2")

BiomarkerCandidates_themed <- BiomarkerCandidates_themed |>
  filter(!grepl(contaminant_pattern, Protein_Description, ignore.case = TRUE)) |>
  filter(!grepl(contaminant_pattern, as.character(Symbol), ignore.case = TRUE)) |>
  filter(!Symbol %in% manual_blacklist)

BiomarkerCandidates_themed = BiomarkerCandidates_themed |>
  arrange(desc(CDS))

# Save
write.csv(
  BiomarkerCandidates_themed,
  "../results/BiomarkerCandidates_themed.csv",
  row.names = FALSE
)
