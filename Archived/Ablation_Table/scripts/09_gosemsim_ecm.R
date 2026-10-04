############################################################
# 09 GOSemSim-ECM benchmark
#
# Ground truth definition (confirmed with user): GOSemSim measured
# against the ECM reference term set specifically (not generic
# GOSemSim), since the AHP-CDS framework's core claim is recovery of
# a coherent ECM invasion module. Reuses the exact ECM GO root terms
# from script 14 (Biological_theme_classification.R) for consistency
# with the rest of the pipeline.
############################################################

suppressMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(dplyr)
})

TOP_N <- 20   # candidates evaluated per condition; RWR uses its own N (see below)

## ---- ECM reference term set (identical to script 14) ----
get_offspring_bp <- function(go_id) tryCatch(c(go_id, as.character(GOBPOFFSPRING[[go_id]])), error = function(e) go_id)
get_offspring_cc <- function(go_id) tryCatch(c(go_id, as.character(GOCCOFFSPRING[[go_id]])), error = function(e) go_id)

ecm_bp_terms <- unique(c(
  get_offspring_bp("GO:0030198"),  # ECM organization
  get_offspring_bp("GO:0007160"),  # cell-matrix adhesion
  get_offspring_bp("GO:0043062")   # extracellular structure organization
))
cat("N ECM reference BP terms:", length(ecm_bp_terms), "\n")

## ---- Build semantic data (local, no internet) ----
semData <- godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE)

## ---- Helper: mean Wang similarity of a gene list vs the ECM term set ----
# For each candidate gene, compare its own GO-BP annotation set to the ECM
# term set using Wang measure with BMA (best-match average) combination,
# then average across the top-N candidates of that condition.
score_condition <- function(entrez_ids, label) {
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids)])
  sims <- sapply(entrez_ids, function(g) {
    go_terms <- tryCatch(semData@geneAnno$GO[semData@geneAnno$ENTREZID == g], error = function(e) NA)
    go_terms <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(g), keytype = "ENTREZID", columns = "GO")
    go_terms <- unique(go_terms$GO[go_terms$ONTOLOGY == "BP" & !is.na(go_terms$GO)])
    if (length(go_terms) == 0) return(NA)
    tryCatch(
      mgoSim(go_terms, ecm_bp_terms, semData = semData, measure = "Wang", combine = "BMA"),
      error = function(e) NA
    )
  })
  cat(sprintf("  %-30s mean=%.4f  (n_scored=%d of %d)\n", label, mean(sims, na.rm = TRUE), sum(!is.na(sims)), length(entrez_ids)))
  data.frame(Condition = label, GOSemSim_ECM = mean(sims, na.rm = TRUE), n_scored = sum(!is.na(sims)))
}

## ---- Load rankings and map to ENTREZID ----
ablation_ranks <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))
rwr_df <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"))

sym_to_entrez <- function(symbols) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENTREZID")
  m[!duplicated(m$SYMBOL), ]
}

map_top <- function(df, rank_col, n = TOP_N) {
  top <- df |> arrange(.data[[rank_col]]) |> head(n)
  em <- sym_to_entrez(top$Symbol)
  merge(top, em, by.x = "Symbol", by.y = "SYMBOL", all.x = TRUE)$ENTREZID
}

cat("\nScoring GOSemSim-ECM for each condition (top", TOP_N, "candidates):\n")

results <- bind_rows(
  score_condition(map_top(ablation_ranks, "rank_AHPCDS"),            "AHP-CDS (full)"),
  score_condition(map_top(ablation_ranks, "rank_FConly"),            "FC-only"),
  score_condition(map_top(ablation_ranks, "rank_Centralityonly"),    "Centrality-only"),
  score_condition(map_top(ablation_ranks, "rank_EqualMCDA"),         "Equal-weight MCDA"),
  score_condition(map_top(rwr_df, "rank_RWR"),                        "Random walk")
)

## ---- Null model: average GOSemSim-ECM across the 9 exact permutations ----
uniprot_to_entrez <- function(uniprot_ids) {
  uniprot_ids <- as.character(uniprot_ids)
  uniprot_ids <- unique(uniprot_ids[!is.na(uniprot_ids) & uniprot_ids != ""])
  
  if (length(uniprot_ids) == 0) return(data.frame(UNIPROT = character(), ENTREZID = character()))
  
  m <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = uniprot_ids, keytype = "UNIPROT", columns = "ENTREZID"),
    error = function(e) data.frame(UNIPROT = character(), ENTREZID = character())
  )
  m <- m[!is.na(m$ENTREZID), ]
  m[!duplicated(m$UNIPROT), ]
}
null_ranks <- readRDS(here::here("Ablation_Table", "data", "null_permutation_ranks_CORRECTED.rds"))
null_scores <- sapply(seq_along(null_ranks), function(k) {
  nr <- null_ranks[[k]] |> arrange(rank_null) |> head(TOP_N)
  em <- uniprot_to_entrez(nr$UNIPROT)
  entrez_ids <- merge(nr, em, by.x = "UNIPROT", by.y = "UNIPROT", all.x = TRUE)$ENTREZID
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids)])
  sims <- sapply(entrez_ids, function(g) {
    go_terms <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(g), keytype = "ENTREZID", columns = "GO")
    go_terms <- unique(go_terms$GO[go_terms$ONTOLOGY == "BP" & !is.na(go_terms$GO)])
    if (length(go_terms) == 0) return(NA)
    tryCatch(mgoSim(go_terms, ecm_bp_terms, semData = semData, measure = "Wang", combine = "BMA"), error = function(e) NA)
  })
  mean(sims, na.rm = TRUE)
})
cat(sprintf("  %-30s mean=%.4f  (range across 9 perms: %.4f-%.4f)\n",
            "Null model (mean of 9 perms)", mean(null_scores), min(null_scores), max(null_scores)))

results <- bind_rows(results, data.frame(
  Condition = "Null model (mean of 9 perms)",
  GOSemSim_ECM = mean(null_scores),
  n_scored = TOP_N
))

print(results)
write.csv(results, here::here("Ablation_Table", "results", "tables", "table3_gosemsim_ecm.csv"), row.names = FALSE)
saveRDS(null_scores, here::here("Ablation_Table", "data", "null_gosemsim_ecm_9perms.rds"))
cat("\nSaved: results/tables/table3_gosemsim_ecm.csv\n")

