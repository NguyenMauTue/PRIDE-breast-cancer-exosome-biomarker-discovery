suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
})

ns <- read.csv(here::here("PXD056161", "results", "tables", "network_summary.csv"))
genes <- ns$Symbol

# same root terms as script 14 (Biological_theme_classification)
get_offspring_bp <- function(go_id) {
  tryCatch(c(go_id, as.character(GOBPOFFSPRING[[go_id]])), error = function(e) go_id)
}
get_offspring_cc <- function(go_id) {
  tryCatch(c(go_id, as.character(GOCCOFFSPRING[[go_id]])), error = function(e) go_id)
}

ecm_cc_terms <- get_offspring_cc("GO:0031012")
ecm_bp_terms <- c(get_offspring_bp("GO:0030198"), get_offspring_bp("GO:0007160"), get_offspring_bp("GO:0043062"))
ms_bp_terms  <- c(get_offspring_bp("GO:0016477"), get_offspring_bp("GO:0000165"), get_offspring_bp("GO:0007265"), get_offspring_bp("GO:0035023"))
ms_cc_terms  <- c(get_offspring_cc("GO:0015629"), get_offspring_cc("GO:0001725"))
vt_bp_terms  <- c(get_offspring_bp("GO:0016192"), get_offspring_bp("GO:0036258"))
vt_cc_terms  <- c(get_offspring_cc("GO:0005765"), get_offspring_cc("GO:0070971"))

all_relevant_go <- unique(c(ecm_cc_terms, ecm_bp_terms, ms_bp_terms, ms_cc_terms, vt_bp_terms, vt_cc_terms))

# map SYMBOL -> GO (local SQLite annotation DB, no internet needed)
go_map <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("GO"))
go_map <- go_map[!is.na(go_map$GO), ]

theme_hits <- go_map |>
  filter(GO %in% all_relevant_go) |>
  distinct(SYMBOL) |>
  pull(SYMBOL)

cat("Genes in network with ECM/Motility-Signaling/Vesicle-Trafficking GO annotation:\n")
cat(length(theme_hits), "of", length(genes), "network genes\n\n")
print(sort(theme_hits))

writeLines(sort(theme_hits), here::here("Ablation_Table", "results", "tables", "go_filtered_seed_candidates.txt"))
