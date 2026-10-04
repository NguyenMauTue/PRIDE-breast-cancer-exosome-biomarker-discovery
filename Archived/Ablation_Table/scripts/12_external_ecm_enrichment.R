############################################################
# 12 External ECM reference set enrichment (MSigDB)
#
# Tests top-N candidates of each condition for enrichment against
# THREE independent, externally-curated ECM gene sets (not derived
# from this study's own pipeline), fetched directly from MSigDB:
#   - NABA_CORE_MATRISOME (275 genes, Naba et al. 2012)
#   - NABA_MATRISOME_ASSOCIATED (751 genes, Naba et al. 2012)
#   - REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION (321 genes, R-HSA-1474244)
#
# This resolves the earlier GOSemSim dilution concern (MCAM/ITGAV-type
# false positives averaging into the score) by using hard set
# membership + Fisher's exact test instead of continuous semantic
# similarity.
############################################################

suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi); library(dplyr) })

core_matrisome <- readLines(here::here("Ablation_Table", "data", "genesets", "naba_core_matrisome.txt"), warn = FALSE)
matrisome_assoc <- readLines(here::here("Ablation_Table", "data", "genesets", "naba_matrisome_associated.txt"), warn = FALSE)
reactome_ecm <- readLines(here::here("Ablation_Table", "data", "genesets", "reactome_ecm_organization.txt"), warn = FALSE)

## ---- Background: full detected proteome ----
de <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
bg_symbols <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = unique(de$UNIPROT), keytype = "UNIPROT", columns = "SYMBOL")$SYMBOL)
bg_symbols <- bg_symbols[!is.na(bg_symbols)]
N_bg <- length(bg_symbols)

fisher_enrich <- function(top_symbols, ref_set, label, ref_label) {
  top_symbols <- unique(top_symbols)
  n <- length(top_symbols)
  K_bg <- sum(bg_symbols %in% ref_set)
  k <- sum(top_symbols %in% ref_set)
  tab <- matrix(c(k, n - k, K_bg - k, N_bg - n - (K_bg - k)), nrow = 2)
  ft <- fisher.test(tab)
  data.frame(Condition = label, Reference = ref_label, n_top = n, k_hits = k,
             pct_hit = round(100 * k / n, 1), odds_ratio = round(unname(ft$estimate), 2),
             p_value = signif(ft$p.value, 3))
}

ablation_ranks <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))
rwr_df <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"))
null_ranks <- readRDS(here::here("Ablation_Table", "data", "null_permutation_ranks.rds"))
lnt <- read.csv(here::here("PXD056161", "results", "tables", "limma_network_table.csv")) |> distinct(UNIPROT, .keep_all = TRUE) |> dplyr::select(UNIPROT, Symbol)
lnt <- lnt[lnt$Symbol %in% ablation_ranks$Symbol, ]

TOP_N <- 20
top_syms <- function(df, rank_col, n = TOP_N) df |> arrange(.data[[rank_col]]) |> head(n) |> pull(Symbol)

refs <- list(NABA_Core = core_matrisome, NABA_Associated = matrisome_assoc, Reactome_ECM = reactome_ecm)

conditions <- list(
  `AHP-CDS (full)`     = top_syms(ablation_ranks, "rank_AHPCDS"),
  `FC-only`            = top_syms(ablation_ranks, "rank_FConly"),
  `Centrality-only`    = top_syms(ablation_ranks, "rank_Centralityonly"),
  `Equal-weight MCDA`  = top_syms(ablation_ranks, "rank_EqualMCDA"),
  `Random walk`        = top_syms(rwr_df, "rank_RWR")
)

all_results <- list()
for (cond_name in names(conditions)) {
  for (ref_name in names(refs)) {
    all_results[[paste(cond_name, ref_name)]] <- fisher_enrich(conditions[[cond_name]], refs[[ref_name]], cond_name, ref_name)
  }
}

# Null model: average across 9 permutations, per reference set
for (ref_name in names(refs)) {
  null_rows <- bind_rows(lapply(seq_along(null_ranks), function(i) {
    nr <- null_ranks[[i]] |> left_join(lnt, by = "UNIPROT")
    fisher_enrich(top_syms(nr, "rank_null"), refs[[ref_name]], paste0("null_perm", i), ref_name)
  }))
  all_results[[paste("Null model", ref_name)]] <- data.frame(
    Condition = "Null model (mean of 9)", Reference = ref_name, n_top = TOP_N,
    k_hits = round(mean(null_rows$k_hits), 1), pct_hit = round(mean(null_rows$pct_hit), 1),
    odds_ratio = round(mean(null_rows$odds_ratio), 2), p_value = signif(mean(null_rows$p_value), 3)
  )
}

results <- bind_rows(all_results) |> mutate(FDR = p.adjust(p_value, method = "BH"))

cat("=== ECM reference enrichment, top-", TOP_N, " candidates, 3 external MSigDB sources ===\n\n", sep = "")
for (ref_name in names(refs)) {
  cat("---", ref_name, "---\n")
  sub <- results |> filter(Reference == ref_name) |> arrange(desc(odds_ratio))
  print(sub[, c("Condition","pct_hit","odds_ratio","p_value","FDR")], row.names = FALSE)
  cat("\n")
}

write.csv(results, here::here("Ablation_Table", "results", "tables", "table3_external_ecm_enrichment.csv"), row.names = FALSE)
cat("Saved: results/tables/table3_external_ecm_enrichment.csv\n")

