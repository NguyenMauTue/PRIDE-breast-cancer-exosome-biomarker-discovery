suppressMessages({
  library(limma); library(dplyr)
  library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
  library(AnnotationDbi); library(GO.db)
  library(igraph); library(httr)
})

TOP_N <- 20
STRING_SPECIES <- 9606
STRING_SCORE   <- 700   # 0.7 on the 0-1000 scale, matches script 10
CALLER_ID      <- "tue1661@gmail.com"  # STRING fair-use requirement

## ============================================================
## Part 0 -- static, label-independent resources (build ONCE)
## ============================================================

imputed_matrix <- readRDS(here::here("PXD056161", "data", "imputed_matrix.rds"))
imputed_matrix <- imputed_matrix[, c("Normal1","Normal2","Normal3","Tumor1","Tumor2","Tumor3")]

# 0a. Fixed ECM term set + fixed 866-protein background 
get_offspring_bp <- function(go_id) tryCatch(c(go_id, as.character(GOBPOFFSPRING[[go_id]])), error = function(e) go_id)
get_offspring_cc <- function(go_id) tryCatch(c(go_id, as.character(GOCCOFFSPRING[[go_id]])), error = function(e) go_id)
ecm_bp_terms <- unique(c(get_offspring_bp("GO:0030198"), get_offspring_bp("GO:0007160"), get_offspring_bp("GO:0043062")))
ecm_cc_terms <- get_offspring_cc("GO:0031012")
all_ecm_terms <- unique(c(ecm_bp_terms, ecm_cc_terms))

de_bg <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
bg_symbols <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = unique(de_bg$UNIPROT), keytype = "UNIPROT", columns = "SYMBOL")$SYMBOL)
bg_symbols <- bg_symbols[!is.na(bg_symbols)]
go_map_bg <- AnnotationDbi::select(org.Hs.eg.db, keys = bg_symbols, keytype = "SYMBOL", columns = "GO")
ecm_genes_bg <- go_map_bg |> filter(GO %in% all_ecm_terms) |> distinct(SYMBOL) |> pull(SYMBOL)
N_bg <- length(bg_symbols); K_bg <- length(ecm_genes_bg)
cat("Fixed background: N =", N_bg, " K_ecm =", K_bg, "\n\n")

fisher_enrich <- function(top_symbols, label) {
  top_symbols <- unique(top_symbols)
  n <- length(top_symbols)
  k <- sum(top_symbols %in% ecm_genes_bg)
  tab <- matrix(c(k, n-k, K_bg-k, N_bg-n-(K_bg-k)), nrow=2)
  ft <- fisher.test(tab)
  data.frame(Condition = label, n_top = n, k_ecm_hits = k, pct_ecm = round(100*k/n,1),
             odds_ratio = round(unname(ft$estimate),2), p_value = signif(ft$p.value,3))
}

# 0b. Static full-proteome annotation table for the contaminant filter
#     (Protein_Description only -- GO columns not needed, script 14 doesn't use them
#     for filtering). Built ONCE over the WHOLE quantified universe, not just the
#     true-label 61-pool, so any permutation's own genes can be looked up.
full_uniprot_universe <- unique(sapply(strsplit(rownames(de_bg), ";"), `[`, 1))
annotations_all <- read.csv(here::here("PXD056161", "data", "annotation_raw.csv"))
static_annotation <- annotations_all |>
  filter(uniprotswissprot %in% full_uniprot_universe) |>
  distinct(external_gene_name, .keep_all = TRUE) |>
  transmute(Symbol = external_gene_name, Protein_Description = description)
cat("Static annotation table built for", nrow(static_annotation), "symbols (full proteome universe)\n\n")
missing_bg <- setdiff(full_uniprot_universe, annotations_all$uniprotswissprot)
if (length(missing_bg) > 0) {
  message(length(missing_bg), " UniProt ID does not have Ensembl translation model (often pseudogene, e.g. Q8IZP2/ST13P4): ",
          paste(missing_bg, collapse = ", "))
}

# 0c. AHP weights
source(here::here("R", "Helper", "ahp_weights.R"))
w <- AHP_WEIGHTS  
norm01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}
contaminant_pattern <- "histone|keratin|actin|tubulin"
manual_blacklist <- c("FMNL1", "H2AX", "H4C6", "H3C1", "H3-3B", "H2AZ2")

## ============================================================
## Part 1 -- STRING REST API network builder, sourced from the
## standalone/testable helper (single source of truth -- do not
## duplicate the function body here, edit string_api_helper.R instead)
## ============================================================
source(here::here("R", "Helper", "string_api_helper.R"))
build_string_network <- local({
  .impl <- build_string_network  # capture the sourced implementation
  function(gene_symbols) {
    .impl(gene_symbols, species = STRING_SPECIES, required_score = STRING_SCORE, caller_identity = CALLER_ID)
  }
})

## ============================================================
## Part 2 -- per-split DE
## ============================================================
source(here::here("R", "Helper", "split_helper.R"))

## ============================================================
## Part 3 -- per-split GSEA -> gene extraction (script 08 + 09 logic)
## ============================================================
genes_from_gsea <- function(de_res, seed) {
  gene_map <- bitr(unique(de_res$UNIPROT), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  rank_df <- merge(de_res, gene_map, by = "UNIPROT")
  rank_df <- rank_df[order(abs(rank_df$t), decreasing = TRUE), ]
  rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]
  gene_list <- sort(setNames(rank_df$t, as.character(rank_df$ENTREZID)), decreasing = TRUE)
  gene_list <- gene_list[!is.na(gene_list)]
  
  set.seed(seed)  # reproducible but split-specific (base seed + split index)
  gsea <- gsePathway(gene_list, organism = "human", minGSSize = 10, maxGSSize = 500,
                     pvalueCutoff = 0.05, verbose = FALSE)
  res <- as.data.frame(gsea)
  if (nrow(res) == 0) return(character(0))
  
  genes <- unique(unlist(strsplit(res$core_enrichment, "/")))
  symbols <- bitr(genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  symbols$SYMBOL
}

## ============================================================
## Part 4 -- build one split's independent candidate pool + score it
## ============================================================
build_pool_and_score <- function(A_idx, B_idx, seed) {
  # Dedup by UNIPROT right away -- protein-group notation ("P12345;Q9Y6K9")
  # collapses to the same first-UNIPROT across multiple original rows after
  # the strsplit(";")[1] step in run_de_for_split(); without this, every join
  # downstream (gene_map, pool) fans out and silently duplicates logFC/FDR
  # onto multiple Symbols. Matches the distinct(UNIPROT, .keep_all=TRUE)
  # pattern already used (twice) in the original 05_null_permutations.R.
  # FIX: run_de_for_split() (from split_helper.R) takes imputed_matrix as its
  # first argument -- calling it with just (A_idx, B_idx) silently shifted
  # every argument over by one and would error with "argument B_idx is missing".
  de <- run_de_for_split(imputed_matrix, A_idx, B_idx) |> distinct(UNIPROT, .keep_all = TRUE)
  gene_symbols <- genes_from_gsea(de, seed)
  net <- build_string_network(gene_symbols)
  if (nrow(net) == 0) return(NULL)
  
  pool <- net |>
    left_join(static_annotation, by = "Symbol") |>
    filter(!grepl(contaminant_pattern, Protein_Description, ignore.case = TRUE) | is.na(Protein_Description)) |>
    filter(!grepl(contaminant_pattern, Symbol, ignore.case = TRUE)) |>
    filter(!Symbol %in% manual_blacklist)
  
  # map DE stats onto the pool via UNIPROT<->Symbol (bitr, same as script 08/09)
  # dedup gene_map by UNIPROT *first* (arbitrary but consistent pick when a
  # UNIPROT maps to >1 symbol), THEN by SYMBOL -- doing SYMBOL-only dedup
  # first (as before) still let one UNIPROT fan out across several Symbol
  # rows and duplicate its logFC/FDR onto each; deduping UNIPROT first
  # guarantees a strict 1 UNIPROT : 1 SYMBOL mapping before the join.
  gene_map <- bitr(unique(de$UNIPROT), fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db) |>
    distinct(UNIPROT, .keep_all = TRUE) |>
    distinct(SYMBOL, .keep_all = TRUE)
  de_sym <- de |> left_join(gene_map, by = "UNIPROT") |> filter(!is.na(SYMBOL)) |> distinct(SYMBOL, .keep_all = TRUE)
  
  scored <- pool |>
    distinct(Symbol, .keep_all = TRUE) |>
    inner_join(de_sym, by = c("Symbol" = "SYMBOL")) |>
    mutate(n_FC = norm01(abs(logFC)), n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
           n_Bet = norm01(log1p(betweenness)), n_Deg = norm01(log1p(degree))) |>
    mutate(CDS_null   = w["FC"]*n_FC + w["FDR"]*n_FDR + w["Bet"]*n_Bet + w["Deg"]*n_Deg,
           FC_null     = n_FC,
           EqualMCDA_null = 0.25*n_FC + 0.25*n_FDR + 0.25*n_Bet + 0.25*n_Deg) |>
    # full ranks across the WHOLE independent pool, not just top-20 --
    # needed downstream for cross-dataset correlation analysis, which
    # compares full rankings, not just enrichment-test candidate sets
    mutate(rank_ahp   = rank(-CDS_null,      ties.method = "first"),
           rank_fc    = rank(-FC_null,       ties.method = "first"),
           rank_equal = rank(-EqualMCDA_null, ties.method = "first"))
  
  list(
    pool_size = nrow(scored),
    full_ranks = scored |>  # full pool, every column downstream scripts might need
      dplyr::select(UNIPROT, Symbol, logFC, adj.P.Val, degree, betweenness,
                    CDS_null, FC_null, EqualMCDA_null, rank_ahp, rank_fc, rank_equal),
    top20_ahp    = scored |> arrange(desc(CDS_null))      |> head(TOP_N) |> pull(Symbol),
    top20_fc     = scored |> arrange(desc(FC_null))       |> head(TOP_N) |> pull(Symbol),
    top20_equal  = scored |> arrange(desc(EqualMCDA_null))|> head(TOP_N) |> pull(Symbol)
  )
}

## ============================================================
## Part 5 -- enumerate splits, run all 9 nulls (split 1 = TRUE, excluded)
## ============================================================
sp <- enumerate_3v3_splits()
true_split  <- sp$true_split   # GroupA = c(1,2,3), i.e. columns Normal1-3, Tumor1-3
null_splits <- sp$null_splits

# Overlap diagnostic: how many samples does each null split's GroupA share
# with the TRUE GroupA? A null split that swaps only 1 sample keeps 2/3 of
# the true grouping intact -> its DE signal (and everything downstream:
# GSEA pathways, STRING network, pool, top-20, OR) will correlate with the
# true result much more than a split that swaps 2-3 samples. This is an
# inherent property of exhaustive 3-vs-3 permutation with n=6 (only 10
# unique partitions exist total), not a code bug -- worth stating explicitly
# as a Methods/Discussion caveat regardless of what the numbers below show.
split_overlap <- sapply(null_splits, function(sp) length(intersect(sp$A, true_split$A)))
cat("True GroupA:", paste(true_split$A, collapse=","), "\n")
cat("Null split overlap with true GroupA (out of 3):", paste(split_overlap, collapse=", "), "\n\n")

BASE_SEED <- 30032026  # same base as script 08, offset per split for reproducible-but-independent GSEA runs

results_ahp <- results_fc <- results_equal <- list()
pool_sizes <- integer(0)
top20_symbols_log <- list()  # keep the actual symbol lists for manual inspection of the 37.21/32.24 coincidences

for (k in seq_along(null_splits)) {
  sp_k <- null_splits[[k]]
  cat("=== Null split", k, "of 9 ===\n")
  Sys.sleep(1)  # be polite to STRING's API across 9 back-to-back calls
  out <- tryCatch(build_pool_and_score(sp_k$A, sp_k$B, seed = BASE_SEED + k),
                  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(out)) next
  cat("  independent pool size:", out$pool_size,
      " | n_top actually used -- AHP:", length(out$top20_ahp),
      " FC:", length(out$top20_fc),
      " Equal:", length(out$top20_equal), "\n")
  pool_sizes[k] <- out$pool_size
  top20_symbols_log[[k]] <- out  # AHP/FC/Equal symbol lists, for later diffing against the true top-20
  results_ahp[[k]]   <- fisher_enrich(out$top20_ahp,   paste0("null_perm", k))
  results_fc[[k]]    <- fisher_enrich(out$top20_fc,    paste0("null_perm", k))
  results_equal[[k]] <- fisher_enrich(out$top20_equal, paste0("null_perm", k))
}

# Combine OR with pool_size / n_top side by side, per condition.
diagnostic_table <- function(res_list, label) {
  bind_rows(res_list) |>
    mutate(pool_size = pool_sizes[seq_along(res_list)],
           overlap_with_true = split_overlap[seq_along(res_list)],
           Method = label)
}
diag_ahp   <- diagnostic_table(results_ahp,   "AHP-CDS")
diag_fc    <- diagnostic_table(results_fc,    "FC-only")
diag_equal <- diagnostic_table(results_equal, "Equal-weight")
cat("\n=== Diagnostic: OR vs n_top vs pool_size vs overlap-with-true, per split ===\n")
print(bind_rows(diag_ahp, diag_fc, diag_equal) |>
        dplyr::select(Method, Condition, overlap_with_true, pool_size, n_top, k_ecm_hits, odds_ratio),
      row.names = FALSE)

## ============================================================
## Part 5b -- full rankings export, for downstream cross-dataset analysis
## (not just the top-20 used by the enrichment test above)
## ============================================================

# Comprehensive: full pool + all 3 methods' ranks, per split
null_full_rankings <- lapply(top20_symbols_log, function(x) if (is.null(x)) NULL else x$full_ranks)
saveRDS(null_full_rankings, here::here("Ablation_Table", "data", "null_permutation_full_rankings_CORRECTED.rds"))
cat("\nSaved full per-split rankings (all 3 methods, whole independent pool) to\n",
    "  data/null_permutation_full_rankings_CORRECTED.rds\n")

# Legacy-compatible: mirrors the original null_permutation_ranks.rds structure
# (list of 9 data.frames with UNIPROT, CDS_null, rank_null) for any old
# downstream script that specifically expects that shape -- AHP-CDS only,
# since that's what the original file scored
null_cds_ranks_legacy_format <- lapply(null_full_rankings, function(df) {
  if (is.null(df)) return(NULL)
  df |> dplyr::select(UNIPROT, CDS_null, rank_null = rank_ahp)
})
saveRDS(null_cds_ranks_legacy_format, here::here("Ablation_Table", "data", "null_permutation_ranks_CORRECTED.rds"))
cat("Saved legacy-format (UNIPROT/CDS_null/rank_null) version to\n",
    "  data/null_permutation_ranks_CORRECTED.rds\n",
    "  (NOTE: any script reading the old null_permutation_ranks.rds should be\n",
    "   repointed at *_CORRECTED.rds -- the original file is still circular)\n\n")

null_ahp_df   <- bind_rows(results_ahp)
null_fc_df    <- bind_rows(results_fc)
null_equal_df <- bind_rows(results_equal)

cat("\n=== Corrected null distributions (9 independently-curated pools) ===\n")
cat("AHP-CDS:       "); print(null_ahp_df$odds_ratio)
cat("FC-only:       "); print(null_fc_df$odds_ratio)
cat("Equal-weight:  "); print(null_equal_df$odds_ratio)

write.csv(null_ahp_df,   here::here("Ablation_Table","results","tables","null_go_ecm_enrichment_9perms_CORRECTED.csv"), row.names = FALSE)
write.csv(null_fc_df,    here::here("Ablation_Table","results","tables","null_go_ecm_enrichment_FConly_9perms_CORRECTED.csv"), row.names = FALSE)
write.csv(null_equal_df, here::here("Ablation_Table","results","tables","null_go_ecm_enrichment_EqualMCDA_9perms_CORRECTED.csv"), row.names = FALSE)
cat("\nSaved 3 corrected null distribution files.\n")

## ---- Compare true-label result (already computed, pool unaffected by this fix)
## against the CORRECTED null. True numbers come from table3_go_ecm_enrichment.csv.
true_row <- read.csv(here::here("Ablation_Table","results","tables","table3_go_ecm_enrichment.csv"))
true_ahp_or <- true_row$odds_ratio[true_row$Condition == "AHP-CDS (full)"]
rank_ahp <- sum(c(true_ahp_or, null_ahp_df$odds_ratio) >= true_ahp_or)
cat("\nAHP-CDS true OR =", true_ahp_or, " | rank among 10 (corrected null) =", rank_ahp,
    " | exact p =", rank_ahp/10, "\n")
