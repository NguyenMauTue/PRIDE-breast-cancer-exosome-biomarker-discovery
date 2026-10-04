# 07_ren2019_gr.R
# ============================================================
# Global Ranking (GR) baseline for Table 3 ablation —
# faithful re-implementation of Ren et al. (2019), "Ranking
# Cancer Proteins by Integrating PPI Network and Protein
# Expression Profiles" (BioMed Research International, 3907195)
#
# Algorithm reproduced directly from their published MATLAB
# source (github.com/Mathlida/protein-ranking, GR_DEPs.m):
#
#   L = I - D^{-1}A          random-walk normalized Laplacian
#   M = I - (alpha/N) * L    one diffusion step
#   K = M^N                  discrete approx. of exp(-alpha*L)
#   p = p0 %*% K             p0 = binary seed indicator
#
# THREE DOCUMENTED DEVIATIONS from the original protocol
# (disclose all three in Methods; do not silently vary):
#
#   (1) NETWORK: our STRING subnetwork (combined_score > 0.7),
#       BINARIZED to match the unweighted property of Ren's
#       HINT interactome. Source and topology still differ from
#       HINT itself -- only "unweighted" is matched, not the
#       actual edge set.
#
#   (2) DEP DEFINITION: Bonferroni correction applied to our
#       limma-derived p-values (de$P.Value), rather than
#       re-deriving DEPs via Welch's t-test as in the original
#       paper. This keeps a single DE pipeline across the whole
#       manuscript while matching Ren's correction stringency.
#       Seed set is the FULL Bonferroni-significant DEP list
#       (no top-N cutoff), same as Ren.
#
#   (3) N = 3 is a coarse truncation of exp(-alpha*L). Kept
#       exactly as reported (not tuned/optimized) to stay
#       faithful to the cited protocol; sensitivity to larger N
#       was not explored.
#
# Requires (already present in environment via run_all.R chain):
#   g   — igraph object, giant component, V(g)$name = SYMBOL
#   de  — full limma DE table (866 obs), incl. de$SYMBOL, de$P.Value
# ============================================================

library(Matrix)
library(igraph)
library(dplyr)
library(here)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)

species = 9606
required_score = 700
caller_identity = "AHP_CDS_NMT"
chunk_size <- 20


de <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
sig <- de[de$adj.P.Val < 0.05, ]
sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys = sig$UNIPROT, keytype = "UNIPROT", columns = "SYMBOL")
sym_map <- sym_map[!duplicated(sym_map$UNIPROT), ]
de_map <- AnnotationDbi::select(org.Hs.eg.db, keys = de$UNIPROT, keytype = "UNIPROT", columns = "SYMBOL")
de_map <- de_map[!duplicated(de_map$UNIPROT), ]
de$SYMBOL = de_map[match(de$UNIPROT, de_map$UNIPROT), "SYMBOL"]

chunks <- split(sym_map$SYMBOL, ceiling(seq_along(sym_map$SYMBOL) / chunk_size))
id_map_list <- lapply(chunks, function(chunk) {
  resp <- httr::POST(
    url = "https://string-db.org/api/tsv/get_string_ids",
    body = list(
      identifiers = paste(chunk, collapse = "\r"),
      species = species, limit = 1, echo_query = 1,
      caller_identity = caller_identity
    ),
    encode = "form"
  )
  httr::stop_for_status(resp)
  df <- read.delim(text = httr::content(resp, as = "text", encoding = "UTF-8"),
                   stringsAsFactors = FALSE)
  Sys.sleep(1)
  df
})
id_map <- do.call(rbind, id_map_list)

if (nrow(id_map) == 0) {
  warning("No STRING IDs resolved for this gene set.")
  return(data.frame(Symbol = character(0), degree = numeric(0), betweenness = numeric(0)))
}

string_ids <- unique(id_map$stringId)

# --- Step 2: pull the high-confidence network for those IDs ---
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

stopifnot(exists("g"), exists("de"))
stopifnot("SYMBOL" %in% colnames(de), "P.Value" %in% colnames(de))

# ---- 1. Seed set: full Bonferroni-significant DEP set, restricted to giant component ----

sig_bonf <- de |>
  mutate(p_bonf = p.adjust(P.Value, method = "bonferroni")) |>
  filter(p_bonf < 0.05)

sig_bonf_graph <- sig_bonf |> filter(SYMBOL %in% V(g)$name)

seed_genes <- unique(sig_bonf_graph$SYMBOL)
n_seed <- length(seed_genes)

cat(sprintf(
  "Seed set (Bonferroni DEP, in giant component): %d / %d proteins in graph\n",
  n_seed, vcount(g)
))

stopifnot(n_seed > 0)  # guard against a degenerate (empty) seed set

# ---- 2. Binarized adjacency + random-walk normalized Laplacian ----

A <- as_adjacency_matrix(g, sparse = TRUE)
A[A > 0] <- 1  # DEVIATION (1): binarize to match Ren's unweighted HINT network
stopifnot(all(A@x == 1))

deg <- Matrix::rowSums(A)
stopifnot(all(deg > 0))  # no isolated nodes should exist in a giant component

Dinv <- Diagonal(x = 1 / deg)
Wm   <- Dinv %*% A
Id   <- Diagonal(nrow(A))
L    <- Id - Wm

# ---- 3. Heat kernel diffusion, alpha = 0.5, N = 3 (verbatim from GR_DEPs.m) ----

alpha <- 0.5
N     <- 3

M <- Id - (alpha / N) * L
K <- M
for (i in seq_len(N - 1)) K <- K %*% M

# ---- 4. Seed vector: binary indicator, NOT normalized to sum = 1 ----

p0 <- setNames(rep(0, vcount(g)), V(g)$name)
p0[seed_genes] <- 1
stopifnot(sum(p0) == n_seed)  # confirms binary, un-normalized seed vector

p0_mat <- Matrix(p0, nrow = 1, sparse = TRUE)
p <- as.vector(p0_mat %*% K)
names(p) <- V(g)$name

# ---- 5. Rank output ----

ren2019_gr_df <- data.frame(
  Symbol    = names(p),
  RWR_score  = p,
  row.names = NULL
) |>
  arrange(desc(RWR_score)) |>
  mutate(rank_RWR = row_number())

# ---- 6. Save ----

out_path <- here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(ren2019_gr_df, out_path, row.names = FALSE)

cat(sprintf(
  "Saved GR ranking (%d proteins) to: %s\n",
  nrow(ren2019_gr_df), out_path
))

write.csv(sig_bonf_graph[, "SYMBOL", drop = FALSE], 
          here::here("Ablation_Table", "results", "tables", "ren2019_seed_genes.csv"), 
          row.names = FALSE)

