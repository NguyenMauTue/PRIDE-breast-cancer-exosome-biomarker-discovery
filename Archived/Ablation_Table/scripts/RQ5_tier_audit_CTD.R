# =============================================================================
# RQ4 Tier Audit — CTD version (no API, no rate limits)
# Replaces DisGeNET-based tier 1/2 assignment entirely. Requires ONE manual
# download step from you — see instructions below — then everything runs
# 100% locally.
#
# Tier 1: curated GDA record for TNBC (CTD DirectEvidence non-empty, disease
#         = Triple-Negative Breast Neoplasms)
# Tier 2: no Tier 1, but curated GDA record for Neoplasms broadly
# Tier 3: STRING-neighbor permutation test (unchanged from before, no API)
# Tier 4: none of the above
# =============================================================================

library(dplyr)
library(here)
library(igraph)

# ---- 0. MANUAL STEP (do this once, outside R) --------------------------------
# 1. Go to https://ctdbase.org/downloads/
# 2. Find "Gene–disease associations" and download CTD_curated_genes_diseases.csv
#    (do NOT unzip — R can read .gz directly)
# 3. Save it into: Ablation_Table/data/raw/CTD_curated_genes_diseases.csv
CTD_PATH <- here::here("Ablation_Table", "data", "CTD_curated_genes_diseases.csv")

if (!file.exists(CTD_PATH)) {
  stop(sprintf(
    "CTD file not found at %s.\nDownload CTD_curated_genes_diseases.csv from https://ctdbase.org/downloads/ and place it there first.",
    CTD_PATH
  ))
}

# ---- 1. Load, skipping comment lines -----------------------------------------
ctd <- read.csv(
  CTD_PATH, comment.char = "#", header = FALSE, stringsAsFactors = FALSE,
  col.names = c("GeneSymbol", "GeneID", "DiseaseName", "DiseaseID",
                "DirectEvidence", "OmimIDs", "PubMedIDs")
)

ctd_curated <- ctd %>% filter(DirectEvidence != "" & !is.na(DirectEvidence))

# ---- 2. Disease definition — keyword-based, CONFIRMED against the uploaded file
# REVISED DESIGN (locked after auditing the actual CTD file):
#   Tier 1 = curated evidence for BREAST cancer broadly (not TNBC-specific —
#            curators almost never tag the narrow "Triple Negative Breast
#            Neoplasms" term; 19 genes vs. 548 under the broad "Breast
#            Neoplasms" term. FN1, for example, has 0 records under the
#            narrow term but 3 under the broad one.)
#   Tier 2 = curated evidence for cancer broadly (any non-breast term
#            matching the cancer keyword list below).
#
# Keyword list was audited term-by-term against every matching MeSH
# DiseaseName in the file to rule out false positives. Excluded on purpose:
#   - "malignant" -> false-positives on "Malignant Hyperthermia" and
#     "Hypertension, Malignant" (neither is cancer)
#   - "tumor"/"tumour" -> false-positives on "Tumoral Calcinosis" (a calcium
#     deposition disorder) and "Tumor Lysis Syndrome" (a treatment
#     complication, not itself a cancer diagnosis)
# Both were dropped. Remaining keywords were spot-checked and found clean.
TIER1_KEYWORDS <- c("neoplasm", "carcinoma")
TIER2_KEYWORDS <- c("neoplasm", "carcinoma", "sarcoma", "melanoma", "leukemia",
                     "lymphoma", "glioma", "myeloma", "blastoma", "mesothelioma")

dname_lower <- tolower(ctd_curated$DiseaseName)
is_breast   <- grepl("breast", dname_lower, fixed = TRUE)
is_tier1_kw <- Reduce(`|`, lapply(TIER1_KEYWORDS, function(k) grepl(k, dname_lower, fixed = TRUE)))
is_tier2_kw <- Reduce(`|`, lapply(TIER2_KEYWORDS, function(k) grepl(k, dname_lower, fixed = TRUE)))

tnbc_genes <- unique(ctd_curated$GeneSymbol[is_breast & is_tier1_kw])
# Tier 2 pool excludes anything already captured by Tier 1's gene set.
neoplasm_genes <- unique(ctd_curated$GeneSymbol[is_tier2_kw & !(ctd_curated$GeneSymbol %in% tnbc_genes)])

cat(sprintf("Tier 1 gene pool (breast cancer, broad, genome-wide): %d (expect 569)\nTier 2 gene pool (any other cancer keyword, genome-wide): %d (expect 3102)\n",
            length(tnbc_genes), length(neoplasm_genes)))
# These counts are genome-wide (all genes in CTD), not restricted to your
# 143 candidates yet — the actual overlap with your candidate list happens
# in step 3 below and will be much smaller.

# ---- 3. Assign Tier 1 / Tier 2 to your 143 candidates ------------------------
candidates <- read.csv(
  here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"),
  stringsAsFactors = FALSE
)

candidates <- candidates %>%
  mutate(
    tier_12 = case_when(
      Symbol %in% tnbc_genes                                  ~ 1L,
      !(Symbol %in% tnbc_genes) & Symbol %in% neoplasm_genes  ~ 2L,
      TRUE                                                     ~ NA_integer_
    )
  )

tier1_symbols <- candidates$Symbol[candidates$tier_12 == 1 & !is.na(candidates$tier_12)]

cat(sprintf("Of your 143 candidates: Tier 1 = %d, Tier 2 = %d, unresolved (Tier 3/4 pending) = %d\n",
            sum(candidates$tier_12 == 1, na.rm = TRUE),
            sum(candidates$tier_12 == 2, na.rm = TRUE),
            sum(is.na(candidates$tier_12))))

# ---- 4. Build igraph object from STRING edges -----
species = 9606
required_score = 700
caller_identity = "AHP_CDS_null_permutation_test"
id_resp <- httr::POST(
  url = "https://string-db.org/api/tsv/get_string_ids",
  body = list(
    identifiers = paste(candidates$Symbol, collapse = "\r"),
    species = species,
    limit = 1,
    echo_query = 1,
    caller_identity = caller_identity
  ),
  encode = "form"
)
httr::stop_for_status(id_resp)
id_map <- read.delim(text = httr::content(id_resp, as = "text", encoding = "UTF-8"),
                     stringsAsFactors = FALSE)

if (nrow(id_map) == 0) {
  warning("No STRING IDs resolved for this gene set.")
  return(data.frame(Symbol = character(0), degree = numeric(0), betweenness = numeric(0)))
}

string_ids <- unique(id_map$stringId)


net_resp <- httr::POST(
  url = "https://string-db.org/api/tsv/network",
  body = list(
    identifiers = paste(string_ids, collapse = "\r"),
    species = species,
    required_score = required_score,
    caller_identity = caller_identity
  ),
  encode = "form"
)
httr::stop_for_status(net_resp)
edges <- read.delim(text = httr::content(net_resp, as = "text", encoding = "UTF-8"),
                    stringsAsFactors = FALSE)

if (nrow(edges) == 0) {
  warning("STRING returned zero edges above required_score for this gene set.")
  return(data.frame(Symbol = gene_symbols, degree = 0, betweenness = 0))
}

# sanity check: API-side filtering should already respect required_score,
# but re-filter defensively in case of API version drift (score is on
# the 0-1000 scale here, same convention as script 10's combined_score)
edges <- edges[edges$score >= required_score / 1000, ]

g <- igraph::graph_from_data_frame(
  edges[, c("preferredName_A", "preferredName_B")],
  directed = FALSE
)
g <- igraph::simplify(g)
# ---- 5. Tier 3: STRING-neighbor of Tier 1, permutation null (unchanged) -----
set.seed(21082026)
n_perm <- 1000
n <- nrow(candidates)
n_tier1 <- length(tier1_symbols)
all_symbols <- igraph::V(g)$name

count_tier1_neighbors <- function(sym, tier1_set) {
  if (!(sym %in% all_symbols)) return(0L)
  nbrs <- names(neighbors(g, sym))
  sum(nbrs %in% tier1_set)
}

candidates$n_tier1_neighbors_real <- sapply(candidates$Symbol, count_tier1_neighbors, tier1_set = tier1_symbols)

perm_matrix <- matrix(NA_integer_, nrow = n, ncol = n_perm)
for (p in seq_len(n_perm)) {
  perm_tier1 <- sample(all_symbols, n_tier1)
  perm_matrix[, p] <- sapply(candidates$Symbol, count_tier1_neighbors, tier1_set = perm_tier1)
}

candidates$tier3_pval <- sapply(seq_len(n), function(i) {
  mean(perm_matrix[i, ] >= candidates$n_tier1_neighbors_real[i])
})

# ---- 6. Final tier assignment ------------------------------------------------
candidates <- candidates %>%
  mutate(
    tier = case_when(
      tier_12 == 1 ~ 1L,
      tier_12 == 2 ~ 2L,
      is.na(tier_12) & n_tier1_neighbors_real > 0 & tier3_pval < 0.05 ~ 3L,
      TRUE ~ 4L
    )
  )

write.csv(candidates, here::here("Ablation_Table", "results", "tables", "candidates_tiered_CTD.csv"), row.names = FALSE)

# ---- 7a. Merge in RWR ranks (Ren et al. 2019 baseline) ----------------------
# RWR was run on its own seed network (291 Bonferroni-DEP genes, not just the
# 143 candidates) to respect the original protocol — see 07_ren2019_gr.R.
# Here we extract the subset relevant to the 143-candidate universe and
# re-rank locally within it, since RQ5's Tier framework and Monte Carlo null
# are both defined on this 143-candidate universe specifically.
rwr_raw <- read.csv(here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"),
                     stringsAsFactors = FALSE)

rwr_in_candidates <- rwr_raw[rwr_raw$Symbol %in% candidates$Symbol, ]
rwr_in_candidates$rank_RWR_143 <- rank(-rwr_in_candidates$RWR_score, ties.method = "min")

# 11 candidates have no STRING edge (score >= 700) within the RWR seed
# network and therefore receive no diffusion signal (confirmed via id_map /
# edges lookup — not a bug). Assign them the lowest rank tier (143) by
# convention, consistent with treatment of degree-zero nodes elsewhere.
isolated_in_candidates <- setdiff(candidates$Symbol, rwr_in_candidates$Symbol)
cat("Candidates with no RWR signal (network-isolated in seed set):",
    length(isolated_in_candidates), "\n")
print(isolated_in_candidates)

isolated_df <- data.frame(Symbol = isolated_in_candidates, rank_RWR_143 = 143)
rwr_combined <- rbind(rwr_in_candidates[, c("Symbol", "rank_RWR_143")], isolated_df)

candidates <- merge(candidates, rwr_combined, by = "Symbol", all.x = TRUE)
stopifnot(sum(is.na(candidates$rank_RWR_143)) == 0)  # every candidate must have a rank now

# ---- 7. G* sweep across k and FOUR rankings (AHP-CDS, FC-only, Equal-weight, RWR)
ks <- c(5, 10, 15, 20, 25, 30)
rank_cols <- c(AHP_CDS = "rank_AHPCDS", FC_only = "rank_FConly",
               Equal_weight = "rank_EqualMCDA", RWR = "rank_RWR_143")

g_star <- function(A, C, N) {
  if (A <= 0 || C <= 0 || N <= 0) return(0)
  3 * (A * C * N)^(1/3)
}

sweep_results <- list()
for (method_name in names(rank_cols)) {
  rc <- rank_cols[[method_name]]
  for (k in ks) {
    topk <- candidates[candidates[[rc]] <= k, ]
    tier_counts <- table(factor(topk$tier, levels = 1:4))
    A <- (tier_counts["1"]) / k
    C <- (tier_counts["2"] + tier_counts["3"]) / k
    N <- (tier_counts["4"]) / k
    sweep_results[[length(sweep_results) + 1]] <- data.frame(
      method = method_name, k = k,
      tier1 = tier_counts["1"], tier2 = tier_counts["2"],
      tier3 = tier_counts["3"], tier4 = tier_counts["4"],
      A = A, C = C, N = N, G_star = g_star(A, C, N)
    )
  }
}
sweep_df <- do.call(rbind, sweep_results)

# ---- 8. Monte Carlo null distribution of G* (10,000 draws per k) -----------
set.seed(21082026)
n_mc <- 10000
tier_vec <- candidates$tier
null_results <- list()

for (k in ks) {
  g_star_null <- numeric(n_mc)
  for (i in seq_len(n_mc)) {
    draw <- sample(tier_vec, k, replace = FALSE)
    tc <- table(factor(draw, levels = 1:4))
    A <- tc["1"] / k; C <- (tc["2"] + tc["3"]) / k; N <- tc["4"] / k
    g_star_null[i] <- g_star(A, C, N)
  }
  null_results[[length(null_results) + 1]] <- data.frame(
    k = k,
    null_mean = mean(g_star_null),
    null_p2_5 = quantile(g_star_null, 0.025),
    null_p97_5 = quantile(g_star_null, 0.975)
  )
}
null_df <- do.call(rbind, null_results)

write.csv(sweep_df, here::here("Ablation_Table", "results", "tables", "RQ5_Gstar_sweep_results.csv"), row.names = FALSE)
write.csv(null_df, here::here("Ablation_Table", "results", "tables", "RQ5_null_distribution.csv"), row.names = FALSE)
print(sweep_df)
print(null_df)
