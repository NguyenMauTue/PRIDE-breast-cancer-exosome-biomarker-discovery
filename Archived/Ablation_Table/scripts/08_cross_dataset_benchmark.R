############################################################
# 08 Cross-dataset benchmark (fractional rank Pearson r)
#
# Ground truth definition:
#   - "Cross-dataset r" = Pearson correlation of FRACTIONAL RANK
#     (rank / N_condition), NOT raw Spearman on overlap proteins.
#     This corrects the earlier survivorship-bias issue where only
#     Spearman on the overlap subset was used (inflates r because
#     it silently drops proteins that "disappeared" due to
#     missingness, which itself correlates with effect size).
#   - Reference = the independently-generated 026 dataset's own
#     CDS-based candidate ranking (Module_tables026.xlsx, "All"
#     sheet), matched by UNIPROT accession.
############################################################

suppressMessages({
  library(dplyr)
  library(openxlsx)
})

frac <- function(rank_vec, N) rank_vec / N

## ============================================================
## ---- Load all 6 condition rankings from PXD012162 ----
## ============================================================
ablation_cross_ranks <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_cross_ranks_partial.csv"))   # AHP-CDS, FC-only, Centrality-only, Equal-weight MCDA
rwr_cross_df         <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_cross_ranks_full.csv"))           # Random walk (needs UNIPROT via Symbol)
N_012 <- nrow(ablation_cross_ranks)

conditions_012 <- ablation_cross_ranks |>
  transmute(
    UNIPROT,
    Symbol,
    frac_AHPCDS_012      = frac(rank_AHPCDS, N_012),
    frac_FConly_012      = frac(rank_FConly, N_012),
    frac_Centrality_012  = frac(rank_Centralityonly, N_012),
    frac_EqualMCDA_012   = frac(rank_EqualMCDA, N_012)
  )

rwr_cross_frac <- rwr_cross_df |>
  mutate(frac_RWR_012 = rank_RWR / nrow(rwr_cross_df)) |>
  dplyr::select(Symbol, frac_RWR_012)
## ============================================================
## ---- Load all 6 condition rankings from PXD056161 ----
## ============================================================
ablation_ranks <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))   # AHP-CDS, FC-only, Centrality-only, Equal-weight MCDA
rwr_df         <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"))           # Random walk (needs UNIPROT via Symbol)
null_ranks     <- readRDS(here::here("Ablation_Table", "data", "null_permutation_ranks_CORRECTED.rds"))    # 9 exact permutations

N_056 <- nrow(ablation_ranks)

conditions_056 <- ablation_ranks |>
  transmute(
    UNIPROT,
    Symbol,
    frac_AHPCDS      = frac(rank_AHPCDS, N_056),
    frac_FConly      = frac(rank_FConly, N_056),
    frac_Centrality  = frac(rank_Centralityonly, N_056),
    frac_EqualMCDA   = frac(rank_EqualMCDA, N_056)
  )

rwr_frac <- rwr_df |>
  mutate(frac_RWR = rank_RWR / nrow(rwr_df)) |>
  dplyr::select(Symbol, frac_RWR)

# Null model: average fractional rank across the 9 exact permutations
null_frac_all <- bind_rows(lapply(seq_along(null_ranks), function(k) {
  nr <- null_ranks[[k]]
  n_k <- nrow(nr)
  data.frame(UNIPROT = nr$UNIPROT, frac_null = nr$rank_null / n_k, perm = k)
}))
null_frac <- null_frac_all |>
  group_by(UNIPROT) |>
  summarise(frac_Null = mean(frac_null), .groups = "drop")

without_RWR_056 <- conditions_056 |>
  left_join(null_frac, by = "UNIPROT")

## ============================================================
## ---- Merge 056 with 012, method-matched, and compute r ----
## ============================================================
merged <- without_RWR_056 |>
  inner_join(conditions_012, by = "UNIPROT") 

merged_rwr <- rwr_frac |>
  inner_join(
    rwr_cross_frac |>
      dplyr::select(Symbol, frac_RWR_012),
    by = "Symbol"
  )

cat("N protein overlap (056 RWR pool vs 012 RWR pool):", nrow(merged_rwr), "\n\n")
cat("N protein overlap (056 ablation pool vs 012 ablation pool):", nrow(merged), "\n\n")

compute_r <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  n <- sum(ok)
  if (n < 4) return(c(r = NA, n = n, p = NA, ci_lo = NA, ci_hi = NA))
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  ci <- ct$conf.int
  c(r = unname(ct$estimate), n = n, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2])
}


results <- tibble::tibble(
  Condition = c("AHP-CDS (full)", "FC-only", "Centrality-only", "Equal-weight MCDA", "Random walk", "Null model (mean of 9 perms)"),
  metric = list(
    compute_r(merged$frac_AHPCDS,     merged$frac_AHPCDS_012),      # method-matched: AHP-CDS vs AHP-CDS
    compute_r(merged$frac_FConly,     merged$frac_FConly_012),      # method-matched: FC-only vs FC-only
    compute_r(merged$frac_Centrality, merged$frac_Centrality_012),  # method-matched: Centrality vs Centrality
    compute_r(merged$frac_EqualMCDA,  merged$frac_EqualMCDA_012),   # method-matched: Equal-weight vs Equal-weight
    compute_r(merged_rwr$frac_RWR,        merged_rwr$frac_RWR_012),         # method-matched: RWR vs RWR
    compute_r(merged$frac_Null,       merged$frac_AHPCDS_012)       # deliberate exception: null vs TRUE AHP-CDS reference
  )
) |>
  mutate(
    cross_dataset_r = sapply(metric, function(m) round(m["r"], 3)),
    n_overlap        = sapply(metric, function(m) m["n"]),
    p_value          = sapply(metric, function(m) round(m["p"], 3)),
    ci_95            = sapply(metric, function(m) paste0("[", round(m["ci_lo"],3), ", ", round(m["ci_hi"],3), "]"))
  ) |>
  dplyr::select(-metric)

print(results)

write.csv(results, here::here("Ablation_Table", "results", "tables", "table3_cross_dataset_r.csv"), row.names = FALSE)
write.csv(merged, here::here("Ablation_Table", "results", "tables", "cross_dataset_merged_fracranks.csv"), row.names = FALSE)
cat("\nSaved: results/tables/table3_cross_dataset_r.csv\n")
