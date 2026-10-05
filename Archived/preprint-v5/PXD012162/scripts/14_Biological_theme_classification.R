############################################################
# 14 Biological theme classification
############################################################
library(dplyr)
library(GO.db)
library(AnnotationDbi)

CDS_candidates = read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_annotated.csv"))
annotated_gene = read.csv(here::here("PXD056161", "results", "tables", "annotated_gene_pool.csv"))
CDS_candidates = CDS_candidates |>
  inner_join(annotated_gene, by = c("UNIPROT" = "uniprotswissprot"))

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
CDS_candidates_themed <-
  CDS_candidates %>%
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

# Biological Module — multi-tag, comma-separated
CDS_candidates_themed$Biological_Module <-
  apply(
    CDS_candidates_themed[, grep("Theme_", colnames(CDS_candidates_themed))],
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

CDS_candidates_themed = CDS_candidates_themed |>
  arrange(desc(CDS))

# Save
write.csv(
  CDS_candidates_themed,
  here::here("PXD056161", "results", "tables", "CDS_candidates_themed.csv"),
  row.names = FALSE
)
