suppressMessages({ library(limma); library(dplyr); library(openxlsx) })

imputed_matrix <- readRDS(here::here("PXD056161", "data", "imputed_matrix.rds"))
imputed_matrix <- imputed_matrix[, c("Normal1","Normal2","Normal3","Tumor1","Tumor2","Tumor3")]

candidates_list <- read.csv(here::here("PXD056161", "results", "tables", "BiomarkerCandidates_themed.csv"))
pool_uniprot <- unique(candidates_list$UNIPROT)   # same 61-protein candidate pool used everywhere else

## ---- Enumerate the same 10 splits as script 05 (split 1 = TRUE) ----
idx <- 1:6
combos <- combn(6, 3)
splits <- list(); seen <- c()
for (i in 1:ncol(combos)) {
  A <- sort(combos[, i]); B <- sort(setdiff(idx, A))
  if (!(1 %in% A)) { tmp <- A; A <- B; B <- tmp }
  key <- paste(A, collapse = ",")
  if (!(key %in% seen)) { seen <- c(seen, key); splits[[length(splits) + 1]] <- list(A = A, B = B) }
}
stopifnot(length(splits) == 10)

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

## ---- FC-only fractional rank for a given split, restricted to the 61-protein pool ----
fconly_fracrank_for_split <- function(A_idx, B_idx) {
  de <- run_de_for_split(A_idx, B_idx)
  de <- de |> filter(UNIPROT %in% pool_uniprot) |> distinct(UNIPROT, .keep_all = TRUE)
  n <- nrow(de)
  de |> arrange(desc(abs(logFC))) |> mutate(frac_FConly = row_number() / n) |>
    dplyr::select(UNIPROT, frac_FConly)
}

## ---- 026 reference ----
module_tables_pxd012162 <- read.xlsx(here::here("PXD012162", "results", "tables", "Module_tables.xlsx"), sheet = "All")
n_pxd012162 <- nrow(module_tables_pxd012162)
module_tables_pxd012162 <- module_tables_pxd012162 |> arrange(desc(CDS)) |> mutate(fracrank_pxd012162 = row_number() / n_pxd012162) |>
  dplyr::select(UNIPROT, fracrank_pxd012162)

cross_r_for_split <- function(A_idx, B_idx) {
  fc <- fconly_fracrank_for_split(A_idx, B_idx)
  m <- inner_join(fc, module_tables_pxd012162, by = "UNIPROT")
  cor(m$frac_FConly, m$fracrank_pxd012162, method = "pearson")
}

cat("Computing FC-only cross-dataset r for all 10 splits (split 1 = TRUE)...\n")
all_r <- sapply(seq_along(splits), function(i) {
  r <- cross_r_for_split(splits[[i]]$A, splits[[i]]$B)
  cat(sprintf("  split %d %s: r = %.4f\n", i, ifelse(i == 1, "(TRUE)", "(null)"), r))
  r
})

true_r_fconly <- all_r[1]
null_r_fconly <- all_r[-1]

cat("\n=== FC-only cross-dataset r: distribution ===\n")
cat("True:  ", round(true_r_fconly, 4), "\n")
cat("Null mean:", round(mean(null_r_fconly), 4), " SD:", round(sd(null_r_fconly), 4), "\n")

rank_of_true <- sum(all_r >= true_r_fconly)
cat("\n=== Exact permutation test (FC-only) ===\n")
cat("True split rank among all 10 (1 = highest):", rank_of_true, "\n")
cat("Exact permutation p-value:", rank_of_true / 10, "\n\n")

## ---- Compare directly to the full-CDS result already computed ----
full_cds <- read.csv(here::here("Ablation_Table", "results", "tables", "null_cross_dataset_r_distribution.csv"))
cat("=== Side-by-side: full-CDS vs FC-only ===\n")
cat("Full-CDS   true r =", round(full_cds$cross_dataset_r[full_cds$is_true], 4),
    " | rank among 10 = 3 (from earlier script)\n")
cat("FC-only    true r =", round(true_r_fconly, 4),
    " | rank among 10 =", rank_of_true, "\n")

out <- data.frame(
  split = c("TRUE (Normal vs Tumor)", paste0("null_perm_", 1:9)),
  cross_dataset_r_FConly = all_r,
  is_true = c(TRUE, rep(FALSE, 9))
)
write.csv(out, here::here("Ablation_Table", "results", "tables", "null_cross_dataset_r_FConly_distribution.csv"), row.names = FALSE)
cat("\nSaved: results/tables/null_cross_dataset_r_FConly_distribution.csv\n")
