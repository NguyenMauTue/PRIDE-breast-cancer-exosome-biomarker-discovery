############################################################
# 13 Compounded Driver score via AHP
############################################################

library(dplyr)
limmanetwork_df = read.csv("results/limma_network_table.csv")


ahp_weights <- function(M) {
  n      <- nrow(M)
  norm   <- sweep(M, 2, colSums(M), "/")  # normalize từng cột
  w      <- rowMeans(norm)                 # row mean = priority vector
  lambda <- mean((M %*% w) / w)
  CI     <- (lambda - n) / (n - 1)
  RI     <- c(0, 0, 0.58, 0.90, 1.12, 1.24)[n]
  CR     <- CI / RI
  list(weights = w, lambda_max = lambda, CI = CI, CR = CR)
}

Pairwise.mat <- matrix(c(
  1,    4,   5,   2.5,
  1/4,  1,   2,   1,
  1/5,  1/2, 1,   1/5,
  1/2.5,1,   5,   1
), nrow = 4, byrow = TRUE,
dimnames = list(c("FC","FDR","Bet","Deg"),
                c("FC","FDR","Bet","Deg")))

ahp   <- ahp_weights(Pairwise.mat)
w     <- ahp$weights
cat("Weights:", round(w, 4), "\n")
cat("CR =", round(ahp$CR, 4), "| OK:", ahp$CR < 0.1, "\n")


norm01 <- function(x) (x - min(x)) / (max(x) - min(x))


limmanetwork_df <- limmanetwork_df |>
  mutate(
    n_FC  = norm01(abs(logFC)),
    n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)),
    n_Deg = norm01(log1p(degree))
  )

limmanetwork_df <- limmanetwork_df |>
  dplyr::mutate(
    CDS = w["FC"]  * n_FC  +
      w["FDR"] * n_FDR +
      w["Bet"] * n_Bet +
      w["Deg"] * n_Deg
  ) |>
  dplyr::arrange(dplyr::desc(CDS))

limmanetwork_df <- limmanetwork_df |>
  dplyr::distinct(UNIPROT, .keep_all = TRUE)


top20 <- limmanetwork_df |>
  dplyr::arrange(dplyr::desc(CDS)) |>
  dplyr::select(UNIPROT, Symbol, logFC, adj.P.Val, degree, betweenness, CDS) |>
  head(20)


# Base ranks
base_ranks <- limmanetwork_df |>
  dplyr::arrange(dplyr::desc(CDS)) |>
  dplyr::mutate(rank_base = dplyr::row_number()) |>
  dplyr::select(UNIPROT, Symbol, rank_base)

# Perturb từng weight, tính CDS mới, lấy rank
perturb_grid <- expand.grid(
  w_idx   = 1:4,
  p_factor = seq(0.8, 1.2, 0.1)
)

all_ranks <- purrr::pmap(perturb_grid, function(w_idx, p_factor) {
  w_p <- w
  w_p[w_idx] <- w_p[w_idx] * p_factor
  w_p <- w_p / sum(w_p)
  
  limmanetwork_df |>
    dplyr::mutate(
      CDS_p = w_p[1]*n_FC + w_p[2]*n_FDR + w_p[3]*n_Bet + w_p[4]*n_Deg
    ) |>
    dplyr::arrange(dplyr::desc(CDS_p)) |>
    dplyr::mutate(rank_p = dplyr::row_number()) |>
    dplyr::select(UNIPROT, Symbol, rank_p)
}) |>
  dplyr::bind_rows() |>
  dplyr::group_by(UNIPROT, Symbol) |>
  dplyr::summarise(
    rank_sd  = round(sd(rank_p), 2),
    rank_min = min(rank_p),
    rank_max = max(rank_p),
    .groups  = "drop"
  ) |>
  dplyr::left_join(base_ranks, by = c("UNIPROT", "Symbol")) |>
  dplyr::arrange(rank_base)

robust_cutoff <- 1.5

output_df <- limmanetwork_df |>
  dplyr::select(UNIPROT, Symbol, logFC, adj.P.Val, degree, betweenness, CDS) |>
  dplyr::left_join(
    all_ranks |> dplyr::select(UNIPROT, rank_base, rank_sd, rank_min, rank_max),
    by = "UNIPROT"
  ) |>
  dplyr::mutate(
    robustness_label = ifelse(rank_sd <= robust_cutoff,
                              "robust_candidate",
                              "weight_sensitive_candidate")
  ) |>
  dplyr::arrange(rank_base)

write.csv(output_df, "results/CDS_candidates_annotated.csv", row.names = FALSE)

