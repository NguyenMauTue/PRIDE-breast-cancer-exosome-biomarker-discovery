############################################################
# 12 Compounded Driver score via AHP
############################################################

library(dplyr)
limma_network_df = read.csv(here::here("PXD056161", "results", "tables", "limma_network_table.csv"))
source(here::here("R", "Helper", "ahp_weights.R"))

norm01 <- function(x) (x - min(x)) / (max(x) - min(x))


limma_network_df <- limma_network_df |>
  mutate(
    n_FC  = norm01(abs(logFC)),
    n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)),
    n_Deg = norm01(log1p(degree))
  )

limma_network_df <- limma_network_df |>
  dplyr::mutate(
    CDS = AHP_WEIGHTS["FC"]  * n_FC  +
      AHP_WEIGHTS["FDR"] * n_FDR +
      AHP_WEIGHTS["Bet"] * n_Bet +
      AHP_WEIGHTS["Deg"] * n_Deg
  ) |>
  dplyr::arrange(dplyr::desc(CDS))

limma_network_df <- limma_network_df |>
  dplyr::distinct(UNIPROT, .keep_all = TRUE)

write.csv(limma_network_df, here::here("PXD056161", "results", "tables", "CDS_candidates_final_result.csv"), row.names = FALSE)

