suppressMessages({ library(dplyr); library(openxlsx) })

## ---- Reference: 026 dataset fractional rank ----
module_tables_pxd012162 <- read.xlsx(here::here("PXD012162", "results", "tables", "Module_tables.xlsx"), sheet = "All")
n_pxd012162 <- nrow(module_tables_pxd012162)
module_tables_pxd012162 <- module_tables_pxd012162 |> arrange(desc(CDS)) |> mutate(fracrank_pxd012162 = row_number() / n_pxd012162) |>
  dplyr::select(UNIPROT, fracrank_pxd012162)

## ---- True split (AHP-CDS) ----
ahp <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))
N_ahp <- nrow(ahp)
true_frac <- ahp |> mutate(frac_AHPCDS = rank_AHPCDS / N_ahp) |> dplyr::select(UNIPROT, frac_AHPCDS)
true_m <- inner_join(true_frac, module_tables_pxd012162, by = "UNIPROT")
true_r <- cor(true_m$frac_AHPCDS, true_m$fracrank_pxd012162, method = "pearson")
cat("True split cross-dataset r:", round(true_r, 4), " (n=", nrow(true_m), ")\n\n")

## ---- 9 null permutations, INDIVIDUAL values (not mean) ----
null_ranks <- readRDS(here::here("Ablation_Table", "data", "null_permutation_ranks_CORRECTED.rds"))
null_rs <- sapply(seq_along(null_ranks), function(k) {
  nr <- null_ranks[[k]]
  n_k <- nrow(nr)
  nf <- nr |> mutate(frac_null = rank_null / n_k) |> dplyr::select(UNIPROT, frac_null)
  m <- inner_join(nf, module_tables_pxd012162, by = "UNIPROT")
  cor(m$frac_null, m$fracrank_pxd012162, method = "pearson")
})

cat("=== Individual null cross-dataset r values (9 permutations) ===\n")
for (i in seq_along(null_rs)) cat(sprintf("  perm %d: r = %.4f\n", i, null_rs[i]))

cat("\n=== Distribution summary ===\n")
cat("Mean:   ", round(mean(null_rs), 4), "\n")
cat("SD:     ", round(sd(null_rs), 4), "\n")
cat("SE:     ", round(sd(null_rs) / sqrt(length(null_rs)), 4), "\n")
cat("Median: ", round(median(null_rs), 4), "\n")
cat("Min:    ", round(min(null_rs), 4), "\n")
cat("Max:    ", round(max(null_rs), 4), "\n")
cat("Range:  ", round(max(null_rs) - min(null_rs), 4), "\n")


# skewness (Fisher-Pearson, no external package needed)
skew <- function(x) {
  n <- length(x); m <- mean(x); s <- sd(x)
  (sum((x - m)^3) / n) / (s^3)
}
cat("Skewness:", round(skew(null_rs), 3), "(0 = symmetric; this is only descriptive with n=9)\n\n")

## ---- Rank-based exact permutation p-value (robust to skew, doesn't use mean) ----
all_10 <- c(true_r, null_rs)
rank_of_true <- sum(all_10 >= true_r)  # 1 = highest of all 10
cat("=== Exact permutation test (rank-based, distribution-free) ===\n")
cat("True split rank among all 10 (1 = highest):", rank_of_true, "\n")
cat("Exact permutation p-value:", rank_of_true / 10, "\n\n")

## ---- z-score of true value relative to null distribution (uses mean+SD; sensitive to n=9 noise, report alongside rank test) ----
z <- (true_r - mean(null_rs)) / sd(null_rs)
cat("z-score of true vs null distribution: ", round(z, 3), "\n")
cat("(For reference only -- with n=9, SD estimate itself has ~30% relative\n")
cat(" uncertainty, so treat this z-score as illustrative, not a formal test.\n")
cat(" The rank-based exact test above is the valid inference for this design.)\n\n")

## ---- Save for plotting ----
out <- data.frame(
  split = c("TRUE (Normal vs Tumor)", paste0("null_perm_", 1:9)),
  cross_dataset_r = c(true_r, null_rs),
  is_true = c(TRUE, rep(FALSE, 9))
)
write.csv(out, here::here("Ablation_Table", "results", "tables", "null_cross_dataset_r_distribution.csv"), row.names = FALSE)
cat("Saved: results/tables/null_cross_dataset_r_distribution.csv\n")
