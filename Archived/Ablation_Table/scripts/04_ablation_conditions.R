library(dplyr)

data056 <- read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_themed.csv"))
#Same normalization method as AHP-CDS
norm01 <- function(x) (x - min(x)) / (max(x) - min(x))


data056 <- data056 |>
  mutate(
    n_FC  = norm01(abs(logFC)),
    n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)),
    n_Deg = norm01(log1p(degree))
  )
## ---- Condition: FC-only ----
# rank purely by |logFC| magnitude (direction-agnostic prioritization, matches
# how AHP-CDS treats FC via n_FC = norm01(abs(logFC)))
fc_only <- data056 |>
  mutate(score_FConly = n_FC) |>
  arrange(desc(score_FConly)) |>
  mutate(rank_FConly = row_number()) |>
  dplyr::select(UNIPROT, score_FConly, rank_FConly)

## ---- Condition: Centrality-only ----
# no expression signal at all
centrality_only <- data056 |>
  mutate(score_Centralityonly = n_Deg) |>
  arrange(desc(score_Centralityonly)) |>
  mutate(rank_Centralityonly = row_number()) |>
  dplyr::select(UNIPROT, score_Centralityonly, rank_Centralityonly)

## ---- Condition: Equal-weight MCDA ----
# same 4 criteria as AHP-CDS (FC, FDR, Bet, Deg) but weights = 1/4 each
# instead of AHP-derived weights -> isolates the value AHP weighting adds
# on top of "just use all 4 criteria"
equal_mcda <- data056 |>
  mutate(score_EqualMCDA = 0.25*n_FC + 0.25*n_FDR + 0.25*n_Bet + 0.25*n_Deg) |>
  arrange(desc(score_EqualMCDA)) |>
  mutate(rank_EqualMCDA = row_number()) |>
  dplyr::select(UNIPROT, score_EqualMCDA, rank_EqualMCDA)

## ---- Reference: AHP-CDS (full), already computed in CDS_candidates_annotated.csv ----
ahp_cds <- data056 |>
  distinct(UNIPROT, .keep_all = TRUE) |>
  arrange(desc(CDS)) |>
  mutate(rank_AHPCDS = row_number()) |>
  dplyr::select(UNIPROT, Symbol, CDS, rank_AHPCDS)

## ---- Merge all into one comparison table ----
ablation_ranks <- ahp_cds |>
  full_join(fc_only, by = "UNIPROT") |>
  full_join(centrality_only, by = "UNIPROT") |>
  full_join(equal_mcda, by = "UNIPROT") |>
  arrange(rank_AHPCDS)

write.csv(ablation_ranks, here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"), row.names = FALSE)

cat("N candidates in pool:", nrow(data056), "\n")
cat("N in AHP-CDS table:", nrow(ahp_cds), "\n\n")
cat("=== Top 10 by AHP-CDS (full) ===\n")
print(head(ablation_ranks[,c("Symbol","rank_AHPCDS","rank_FConly","rank_Centralityonly","rank_EqualMCDA")], 10))

cat("\n=== Rank agreement (Spearman) between AHP-CDS and each ablated condition ===\n")
cat("AHP-CDS vs FC-only:         ", cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_FConly, method="spearman", use="pairwise.complete.obs"), "\n")
cat("AHP-CDS vs Centrality-only: ", cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_Centralityonly, method="spearman", use="pairwise.complete.obs"), "\n")
cat("AHP-CDS vs Equal-weight MCDA:", cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_EqualMCDA, method="spearman", use="pairwise.complete.obs"), "\n")

