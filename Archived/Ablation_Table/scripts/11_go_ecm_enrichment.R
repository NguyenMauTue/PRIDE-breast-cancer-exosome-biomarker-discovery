############################################################
# 11 GO-based ECM enrichment (Odds Ratio + Fisher's exact test)
#
# Complements GOSemSim (continuous semantic similarity, shown to be
# diluted by loosely-related "cell adhesion" genes like MCAM/ITGAV)
# with a hard membership-based test: is a top-N candidate list
# significantly ENRICHED for genes annotated to the ECM GO terms,
# relative to the full detected-protein background? This is the
# standard over-representation analysis (ORA) approach and does not
# have GOSemSim's averaging/dilution problem.
############################################################

suppressMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(GO.db); library(dplyr)
})

get_offspring_bp <- function(go_id) tryCatch(c(go_id, as.character(GOBPOFFSPRING[[go_id]])), error = function(e) go_id)
get_offspring_cc <- function(go_id) tryCatch(c(go_id, as.character(GOCCOFFSPRING[[go_id]])), error = function(e) go_id)

ecm_bp_terms <- unique(c(get_offspring_bp("GO:0030198"), get_offspring_bp("GO:0007160"), get_offspring_bp("GO:0043062")))
ecm_cc_terms <- get_offspring_cc("GO:0031012")  # extracellular matrix (CC)
all_ecm_terms <- unique(c(ecm_bp_terms, ecm_cc_terms))

## ---- Background: full detected proteome (866 proteins from script 03_de.R) ----
de <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
bg_symbols <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = unique(de$UNIPROT), keytype = "UNIPROT", columns = "SYMBOL")$SYMBOL)
bg_symbols <- bg_symbols[!is.na(bg_symbols)]

## ---- Flag ECM membership for the whole background once ----
go_map_bg <- AnnotationDbi::select(org.Hs.eg.db, keys = bg_symbols, keytype = "SYMBOL", columns = "GO")
ecm_genes_bg <- go_map_bg |> filter(GO %in% all_ecm_terms) |> distinct(SYMBOL) |> pull(SYMBOL)

N_bg <- length(bg_symbols)
K_bg <- length(ecm_genes_bg)  # total ECM-annotated genes in background
cat("Background: N =", N_bg, " proteins, K =", K_bg, "ECM-annotated (", round(100*K_bg/N_bg,1), "%)\n\n")

fisher_enrich <- function(top_symbols, label) {
  top_symbols <- unique(top_symbols)
  n <- length(top_symbols)
  k <- sum(top_symbols %in% ecm_genes_bg)  # ECM genes in top-N
  # 2x2 table: [in_top & ECM, in_top & not-ECM; not_top & ECM, not_top & not-ECM]
  tab <- matrix(c(k, n-k, K_bg-k, N_bg-n-(K_bg-k)), nrow=2)
  ft <- fisher.test(tab)
  data.frame(Condition = label, n_top = n, k_ecm_hits = k, pct_ecm = round(100*k/n,1),
             odds_ratio = round(unname(ft$estimate),2), p_value = signif(ft$p.value,3))
}

## ---- Run for all 6 conditions, top-20 candidates ----
ablation_ranks <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))
rwr_df <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"))
null_ranks <- readRDS(here::here("Ablation_Table", "data", "null_permutation_ranks.rds"))
lnt <- read.csv(here::here("PXD056161", "results", "tables", "limma_network_table.csv")) |> distinct(UNIPROT, .keep_all = TRUE) |> dplyr::select(UNIPROT, Symbol)
lnt <- lnt[lnt$Symbol %in% ablation_ranks$Symbol, ]

TOP_N <- 20
top_syms <- function(df, rank_col, n = TOP_N) df |> arrange(.data[[rank_col]]) |> head(n) |> pull(Symbol)

results <- bind_rows(
  fisher_enrich(top_syms(ablation_ranks, "rank_AHPCDS"), "AHP-CDS (full)"),
  fisher_enrich(top_syms(ablation_ranks, "rank_FConly"), "FC-only"),
  fisher_enrich(top_syms(ablation_ranks, "rank_Centralityonly"), "Centrality-only"),
  fisher_enrich(top_syms(ablation_ranks, "rank_EqualMCDA"), "Equal-weight MCDA"),
  fisher_enrich(top_syms(rwr_df, "rank_RWR"), "Random walk")
)

null_results <- bind_rows(lapply(seq_along(null_ranks), function(i) {
  nr <- null_ranks[[i]] |> left_join(lnt, by = "UNIPROT")
  fisher_enrich(top_syms(nr, "rank_null"), paste0("null_perm", i))
}))
null_summary <- data.frame(
  Condition = "Null model (mean of 9)",
  n_top = TOP_N,
  k_ecm_hits = round(mean(null_results$k_ecm_hits), 1),
  pct_ecm = round(mean(null_results$pct_ecm), 1),
  odds_ratio = round(mean(null_results$odds_ratio), 2),
  p_value = signif(mean(null_results$p_value), 3)
)

results <- bind_rows(results, null_summary) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

cat("=== GO ECM enrichment, top-", TOP_N, " candidates per condition ===\n", sep="")
print(results, digits = 3)

write.csv(results, here::here("Ablation_Table", "results", "tables", "table3_go_ecm_enrichment.csv"), row.names = FALSE)
write.csv(null_results, here::here("Ablation_Table", "results", "tables", "null_go_ecm_enrichment_9perms.csv"), row.names = FALSE)
cat("\nSaved: results/tables/table3_go_ecm_enrichment.csv\n")
