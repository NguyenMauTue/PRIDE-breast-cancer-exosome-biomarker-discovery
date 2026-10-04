suppressMessages({
  library(limma)
  library(dplyr)
})

imputed_matrix <- readRDS(here::here("PXD056161", "data", "imputed_matrix.rds"))   # 6 samples: Normal1-3, Tumor1-3
network_summary <- read.csv(here::here("PXD056161", "results", "tables", "network_summary.csv"))
candidates_list <- read.csv(here::here("PXD056161", "results", "tables", "BiomarkerCandidates_themed.csv"))
network_summary <- network_summary[network_summary$Symbol %in% candidates_list$Symbol, ]
network_summary <- network_summary %>%
  left_join(
    candidates_list %>% dplyr::select(Symbol, UNIPROT) %>% distinct(Symbol, .keep_all = TRUE), 
    by = "Symbol"
  )


samples <- colnames(imputed_matrix)
stopifnot(identical(sort(samples), sort(c("Normal1","Normal2","Normal3","Tumor1","Tumor2","Tumor3"))))
imputed_matrix <- imputed_matrix[, c("Normal1","Normal2","Normal3","Tumor1","Tumor2","Tumor3")]

# enumerate the 10 unique 3v3 splits (sample "Normal1" anchored to side A to avoid double counting)
idx <- 1:6
combos <- combn(6, 3)
splits <- list()
seen <- c()
for (i in 1:ncol(combos)) {
  A <- sort(combos[, i]); B <- sort(setdiff(idx, A))
  if (!(1 %in% A)) { tmp <- A; A <- B; B <- tmp }
  key <- paste(A, collapse = ",")
  if (!(key %in% seen)) { seen <- c(seen, key); splits[[length(splits) + 1]] <- list(A = A, B = B) }
}
stopifnot(length(splits) == 10)

# splits[[1]] is the TRUE split (Normal1-3 vs Tumor1-3) -- exclude it, keep the 9 null splits
null_splits <- splits[-1]
cat("Running", length(null_splits), "null permutations (exact enumeration, n=3v3)\n")

# Same UNIPROT dedup + Symbol canonicalization used for the other conditions
symbol_map <- network_summary |> distinct(Symbol, .keep_all = TRUE)

norm01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# Recompute the exact AHP weights from script 12's own pairwise comparison
# matrix, rather than hardcoding -- guarantees the null model uses the
# identical weighting scheme as the real AHP-CDS, only the labels differ.
ahp_weights <- function(M) {
  n    <- nrow(M)
  norm <- sweep(M, 2, colSums(M), "/")
  w    <- rowMeans(norm)
  list(weights = w)
}
Pairwise.mat <- matrix(c(
  1,    4,   5,   2.5,
  1/4,  1,   2,   1,
  1/5,  1/2, 1,   1/5,
  1/2.5,1,   5,   1
), nrow = 4, byrow = TRUE,
dimnames = list(c("FC","FDR","Bet","Deg"), c("FC","FDR","Bet","Deg")))
w <- ahp_weights(Pairwise.mat)$weights
cat("AHP weights used for null model (matches script 12):", round(w, 4), "\n\n")

run_de_for_split <- function(A_idx, B_idx) {
  group <- factor(rep(NA, 6), levels = c("GroupA", "GroupB"))
  group[A_idx] <- "GroupA"; group[B_idx] <- "GroupB"
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(imputed_matrix, design)
  cm <- makeContrasts(GroupB_VS_GroupA = GroupB - GroupA, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  res <- topTable(fit2, coef = "GroupB_VS_GroupA", number = Inf, sort.by = "none")
  res$X <- rownames(res)
  res <- res |> mutate(UNIPROT = sapply(strsplit(X, ";"), `[`, 1))
  res
}

null_cds_ranks <- list()

for (k in seq_along(null_splits)) {
  sp <- null_splits[[k]]
  de_null <- run_de_for_split(sp$A, sp$B)

  # restrict to the same 95-row candidate pool (STRING network / GSEA-derived),
  # matched by UNIPROT, then attach degree/betweenness which are expression-
  # independent (network topology doesn't change under label shuffle)
  lnt_null <- de_null |>
    inner_join(
      network_summary |>
        dplyr::select(UNIPROT, degree, betweenness) |>
        distinct(UNIPROT, .keep_all = TRUE),
      by = "UNIPROT"
    ) |>
    distinct(UNIPROT, .keep_all = TRUE)

  lnt_null <- lnt_null |>
    mutate(
      n_FC  = norm01(abs(logFC)),
      n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
      n_Bet = norm01(log1p(betweenness)),
      n_Deg = norm01(log1p(degree))
    )

  lnt_null <- lnt_null |>
    mutate(CDS_null = w["FC"]*n_FC + w["FDR"]*n_FDR + w["Bet"]*n_Bet + w["Deg"]*n_Deg) |>
    arrange(desc(CDS_null)) |>
    mutate(rank_null = row_number())

  null_cds_ranks[[k]] <- lnt_null |> dplyr::select(UNIPROT, CDS_null, rank_null)
  cat("split", k, "of 9 done -- N candidates:", nrow(lnt_null), "\n")
}

saveRDS(null_cds_ranks, here::here("Ablation_Table", "data", "null_permutation_ranks.rds"))
cat("\nSaved 9 null permutation rankings to data/null_permutation_ranks.rds\n")
