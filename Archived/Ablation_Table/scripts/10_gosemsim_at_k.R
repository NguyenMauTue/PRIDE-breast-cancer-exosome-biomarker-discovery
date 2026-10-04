suppressMessages({
  library(GOSemSim); library(org.Hs.eg.db); library(AnnotationDbi); library(GO.db); library(dplyr)
})

get_offspring_bp <- function(go_id) tryCatch(c(go_id, as.character(GOBPOFFSPRING[[go_id]])), error = function(e) go_id)
ecm_bp_terms <- unique(c(get_offspring_bp("GO:0030198"), get_offspring_bp("GO:0007160"), get_offspring_bp("GO:0043062")))
semData <- godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE)

per_gene_ecm_sim <- function(symbols) {
  em <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENTREZID")
  em <- em[!duplicated(em$SYMBOL) & !is.na(em$ENTREZID), ]
  sims <- sapply(em$ENTREZID, function(g) {
    go_terms <- AnnotationDbi::select(org.Hs.eg.db, keys = g, keytype = "ENTREZID", columns = "GO")
    go_terms <- unique(go_terms$GO[go_terms$ONTOLOGY == "BP" & !is.na(go_terms$GO)])
    if (length(go_terms) == 0) return(NA)
    tryCatch(mgoSim(go_terms, ecm_bp_terms, semData = semData, measure = "Wang", combine = "BMA"), error = function(e) NA)
  })
  setNames(sims, em$SYMBOL)
}

ablation_ranks <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))
rwr_df <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"))
null_ranks <- readRDS(here::here("Ablation_Table", "data", "null_permutation_ranks.rds"))
lnt <- read.csv(here::here("PXD056161", "results", "tables", "limma_network_table.csv")) |> distinct(UNIPROT, .keep_all = TRUE) |> dplyr::select(UNIPROT, Symbol)
lnt <- lnt[lnt$Symbol %in% ablation_ranks$Symbol, ]
# per-gene similarity, computed once for every distinct gene across all conditions
all_symbols <- unique(c(ablation_ranks$Symbol, rwr_df$Symbol))
sim_lookup <- per_gene_ecm_sim(all_symbols)

curve_for <- function(df, rank_col, label, max_k = 40) {
  ord <- df |> arrange(.data[[rank_col]]) |> pull(Symbol)
  ord <- ord[ord %in% names(sim_lookup)]
  s <- sim_lookup[ord]
  s[is.na(s)] <- NA
  cum_mean <- sapply(1:min(max_k, length(s)), function(k) mean(s[1:k], na.rm = TRUE))
  data.frame(k = 1:length(cum_mean), GOSemSim_at_k = cum_mean, Condition = label)
}

max_k <- 40
curves <- bind_rows(
  curve_for(ablation_ranks, "rank_AHPCDS", "AHP-CDS (full)", max_k),
  curve_for(ablation_ranks, "rank_FConly", "FC-only", max_k),
  curve_for(ablation_ranks, "rank_Centralityonly", "Centrality-only", max_k),
  curve_for(ablation_ranks, "rank_EqualMCDA", "Equal-weight MCDA", max_k),
  curve_for(rwr_df, "rank_RWR", "Random walk", max_k)
)

# null: average GOSemSim@k across the 9 permutations
null_curves <- bind_rows(lapply(seq_along(null_ranks), function(i) {
  nr <- null_ranks[[i]] |> left_join(lnt, by = "UNIPROT")
  curve_for(nr, "rank_null", paste0("null_perm", i), max_k)
}))
null_mean_curve <- null_curves |> group_by(k) |> summarise(GOSemSim_at_k = mean(GOSemSim_at_k, na.rm = TRUE)) |> mutate(Condition = "Null model (mean of 9)")

curves <- bind_rows(curves, null_mean_curve)
write.csv(curves, here::here("Ablation_Table", "results", "tables", "gosemsim_at_k_curves.csv"), row.names = FALSE)

cat("GOSemSim@k at k = 5, 10, 15, 20, 25, 30:\n")
checkpoints <- curves |> filter(k %in% c(5,10,15,20,25,30)) |>
  tidyr::pivot_wider(names_from = Condition, values_from = GOSemSim_at_k)
print(as.data.frame(checkpoints), digits = 3)
