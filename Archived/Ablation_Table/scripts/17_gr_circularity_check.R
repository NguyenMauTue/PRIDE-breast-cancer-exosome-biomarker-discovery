# 17_gr_circularity_check
# ============================================================
# Table 3 robustness check: does the AHP-CDS <-> GR (Ren 2019)
# rank correlation reflect genuine convergence, or is it
# inflated by circularity (GR seed proteins are, by
# construction, ranked highest by GR -- and many are also
# top-ranked by AHP-CDS because both use DE significance as
# an input)?
#
# Approach: compute Spearman correlation between rank_AHPCDS
# and rank_RWR (GR) twice --
#   (a) on the FULL overlap between the AHP-CDS candidate pool
#       and the GR-ranked giant component
#   (b) on the SAME set EXCLUDING any protein used as a GR seed
# A large drop from (a) to (b) would indicate the full-set
# correlation is substantially circularity-driven.
#
# 95% CI via Fisher z-transform. NOTE: Fisher z is exact for
# Pearson r under bivariate normality; applied to Spearman rho
# it is a widely used but approximate large-sample method
# (adequate here for descriptive CI, not a formal test).
#
# Inputs (paths follow here::here() convention, PR1):
#   PXD012162/results/tables/ablation_ranks_partial.csv
#     -- must contain: Symbol, rank_AHPCDS
#   PXD012162/results/tables/rwr_ranks_full.csv
#     -- GR (Ren 2019) ranking; filename kept as "rwr_ranks_full.csv"
#        for downstream compatibility even though the algorithm
#        is Ren et al. (2019) GR, not classic RWR -- see
#        07_ren2019_gr.R header for full deviation disclosure.
#     -- must contain: Symbol, rank_RWR
#   PXD012162/results/tables/ren2019_seed_genes.csv
#     -- must contain: SYMBOL (the Bonferroni-DEP seed set used
#        to build the GR ranking)
#
# Output:
#   PXD012162/results/tables/table3_gr_circularity_check.csv
#     -- two rows: "Full set" and "Seed excluded", each with
#        rho, p_value, n, ci_lower, ci_upper
# ============================================================

library(dplyr)
library(here)

# ---- 1. Load ----

ahp  <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"),
                 stringsAsFactors = FALSE)
gr   <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"),
                 stringsAsFactors = FALSE)
seed <- read.csv(here::here("Ablation_Table", "results", "tables", "ren2019_seed_genes.csv"),
                 stringsAsFactors = FALSE)

stopifnot(all(c("Symbol", "rank_AHPCDS") %in% colnames(ahp)))
stopifnot(all(c("Symbol", "rank_RWR")    %in% colnames(gr)))
stopifnot("SYMBOL" %in% colnames(seed))

seed_set <- unique(seed$SYMBOL)
n_seed   <- length(seed_set)
cat(sprintf("GR seed set: %d proteins\n", n_seed))

# ---- 2. Merge AHP-CDS candidate pool with GR ranking ----

merged <- ahp |>
  inner_join(gr, by = "Symbol") |>
  mutate(is_seed = Symbol %in% seed_set)

n_pool     <- nrow(ahp)
n_overlap  <- nrow(merged)
n_missing  <- n_pool - n_overlap

cat(sprintf(
  "AHP-CDS candidate pool: %d | overlapping with GR giant component: %d | outside giant component: %d\n",
  n_pool, n_overlap, n_missing
))
if (n_missing > 0) {
  cat("Candidates outside GR giant component (excluded from this comparison):\n")
  print(ahp$Symbol[!ahp$Symbol %in% gr$Symbol])
}

# ---- 3. Spearman + Fisher z 95% CI, helper ----

spearman_with_ci <- function(x, y, label) {
  n   <- length(x)
  fit <- suppressWarnings(cor.test(x, y, method = "spearman"))
  rho <- unname(fit$estimate)
  pval <- fit$p.value
  
  # Fisher z-transform CI (approximate for Spearman; see header note)
  z  <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1, 1) * qnorm(0.975) * se)
  
  data.frame(
    Condition = label,
    rho       = rho,
    p_value   = pval,
    n         = n,
    ci_lower  = ci[1],
    ci_upper  = ci[2]
  )
}

# ---- 4. Full set vs seed-excluded ----

full_row <- spearman_with_ci(merged$rank_AHPCDS, merged$rank_RWR, "Full set")

non_seed <- merged |> filter(!is_seed)
n_seed_in_overlap <- sum(merged$is_seed)
cat(sprintf(
  "Of %d overlapping candidates, %d (%.0f%%) are GR seed proteins\n",
  n_overlap, n_seed_in_overlap, 100 * n_seed_in_overlap / n_overlap
))

if (nrow(non_seed) >= 4) {  # need n>=4 for Fisher z SE (n-3 > 0) to be meaningful
  excl_row <- spearman_with_ci(non_seed$rank_AHPCDS, non_seed$rank_RWR, "Seed excluded")
} else {
  warning("Too few non-seed overlapping proteins for a stable CI; reporting NA row.")
  excl_row <- data.frame(
    Condition = "Seed excluded", rho = NA, p_value = NA,
    n = nrow(non_seed), ci_lower = NA, ci_upper = NA
  )
}

result <- bind_rows(full_row, excl_row)
print(result)

# ---- 5. Save ----

out_path <- here::here("Ablation_Table", "results", "tables", "table3_gr_circularity_check.csv")
write.csv(result, out_path, row.names = FALSE)
cat(sprintf("\nSaved circularity check table to: %s\n", out_path))
