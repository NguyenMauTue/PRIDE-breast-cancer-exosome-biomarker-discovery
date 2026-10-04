############################################################
# 14b_fconly_go_ecm_null.R
#
# FC-only equivalent of the GO-ECM enrichment null test that
# produced null_go_ecm_enrichment_9perms.csv (AHP-CDS, from script
# 11_go_ecm_enrichment.R). Reuses fconly_top20_for_split() logic from
# 14_fconly_null_test.R exactly (same 10-split enumeration, same DE
# pipeline, same 61-protein candidate pool restriction), and reuses
# fisher_enrich() logic verbatim from script 11 (same ECM term set,
# same N=866 full-detected-proteome background) so the two null
# tests are directly comparable.
############################################################

suppressMessages({
  library(limma); library(dplyr); library(openxlsx)
  library(org.Hs.eg.db); library(AnnotationDbi); library(GO.db)
})

TOP_N <- 20   # matches AHP-CDS top-N (script 11)

## ============================================================
## Part A — ECM term set + background (verbatim from script 11)
## ============================================================
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

N_bg <- length(bg_symbols)
K_bg <- length(ecm_genes_bg)
cat("Background: N =", N_bg, " proteins, K =", K_bg, "ECM-annotated (", round(100*K_bg/N_bg,1), "%)\n\n")

fisher_enrich <- function(top_symbols, label) {
  top_symbols <- unique(top_symbols)
  n <- length(top_symbols)
  k <- sum(top_symbols %in% ecm_genes_bg)
  tab <- matrix(c(k, n-k, K_bg-k, N_bg-n-(K_bg-k)), nrow=2)
  ft <- fisher.test(tab)
  data.frame(Condition = label, n_top = n, k_ecm_hits = k, pct_ecm = round(100*k/n,1),
             odds_ratio = round(unname(ft$estimate),2), p_value = signif(ft$p.value,3))
}

## ============================================================
## Part B — FC-only splits (verbatim from script 14)
## ============================================================
imputed_matrix <- readRDS(here::here("PXD056161", "data", "imputed_matrix.rds"))
imputed_matrix <- imputed_matrix[, c("Normal1","Normal2","Normal3","Tumor1","Tumor2","Tumor3")]

candidates_list <- read.csv(here::here("PXD056161", "results", "tables", "BiomarkerCandidates_themed.csv"))
pool_uniprot <- unique(candidates_list$UNIPROT)   # same 61-protein candidate pool used everywhere else

lnt <- read.csv(here::here("PXD056161", "results", "tables", "limma_network_table.csv")) |>
  distinct(UNIPROT, .keep_all = TRUE) |> dplyr::select(UNIPROT, Symbol)

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

## FC-only top-20 Symbols for a given split, restricted to the 61-protein pool
fconly_top20_symbols_for_split <- function(A_idx, B_idx, n = TOP_N) {
  de <- run_de_for_split(A_idx, B_idx)
  de <- de |> filter(UNIPROT %in% pool_uniprot) |> distinct(UNIPROT, .keep_all = TRUE)
  top <- de |> arrange(desc(abs(logFC))) |> head(n)
  top <- top |> left_join(lnt, by = "UNIPROT")
  top$Symbol
}

## ============================================================
## Part C — run all 10 splits, score with fisher_enrich
## ============================================================
cat("Computing FC-only GO-ECM enrichment OR for all 10 splits (split 1 = TRUE)...\n")
null_results <- bind_rows(lapply(seq_along(splits), function(i) {
  top_symbols <- fconly_top20_symbols_for_split(splits[[i]]$A, splits[[i]]$B)
  label <- if (i == 1) "TRUE (Normal vs Tumor)" else paste0("null_perm", i - 1)
  res <- fisher_enrich(top_symbols, label)
  cat(sprintf("  split %d %s: OR = %.2f, k_ecm_hits = %d\n", i, ifelse(i == 1, "(TRUE)", "(null)"), res$odds_ratio, res$k_ecm_hits))
  res
}))

true_or <- null_results$odds_ratio[null_results$Condition == "TRUE (Normal vs Tumor)"]
all_or <- null_results$odds_ratio
rank_of_true <- sum(all_or >= true_or)

cat("\n=== Exact permutation test (FC-only, GO-ECM OR) ===\n")
cat("True split OR:", true_or, "\n")
cat("True split rank among all 10 (1 = highest):", rank_of_true, "\n")
cat("Exact permutation p-value:", rank_of_true / 10, "\n\n")

## ---- Side-by-side with AHP-CDS's already-computed null (script 11) ----
ahp_null <- read.csv(here::here("Ablation_Table", "results", "tables", "null_go_ecm_enrichment_9perms.csv"))
ahp_true_row <- read.csv(here::here("Ablation_Table", "results", "tables", "table3_go_ecm_enrichment.csv"))
ahp_true_or <- ahp_true_row$odds_ratio[ahp_true_row$Condition == "AHP-CDS (full)"]
ahp_all_or <- c(ahp_true_or, ahp_null$odds_ratio)
ahp_rank <- sum(ahp_all_or >= ahp_true_or)

cat("=== Side-by-side: AHP-CDS full vs FC-only (GO-ECM OR null) ===\n")
cat("AHP-CDS full   true OR =", ahp_true_or, " | rank among 10 =", ahp_rank, "\n")
cat("FC-only        true OR =", true_or,     " | rank among 10 =", rank_of_true, "\n")

write.csv(null_results, here::here("Ablation_Table", "results", "tables", "null_go_ecm_enrichment_FConly_9perms.csv"), row.names = FALSE)
cat("\nSaved: results/tables/null_go_ecm_enrichment_FConly_9perms.csv\n")