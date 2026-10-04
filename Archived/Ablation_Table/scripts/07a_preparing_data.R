############################################################
#
# 07a_Preparing data for 08_cross_dataset_benchmark.R
#
############################################################
library(dplyr)
library(Matrix)
library(igraph)
library(dplyr)
library(here)
library(org.Hs.eg.db)

species = 9606
required_score = 700
caller_identity = "AHP_CDS_null_permutation_test"

data012 <- read.csv(here::here("PXD012162", "results", "tables", "CDS_candidates_annotated.csv"))
norm01 <- function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))

data012 <- data012 |>
  mutate(
    n_FC  = norm01(abs(logFC)),
    n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)),
    n_Deg = norm01(log1p(degree))
  )
## ---- Condition: FC-only ----
# rank purely by |logFC| magnitude (direction-agnostic prioritization, matches
# how AHP-CDS treats FC via n_FC = norm01(abs(logFC)))
fc_only <- data012 |>
  mutate(score_FConly = n_FC) |>
  arrange(desc(score_FConly)) |>
  mutate(rank_FConly = row_number()) |>
  dplyr::select(UNIPROT, score_FConly, rank_FConly)

## ---- Condition: Centrality-only ----
centrality_only <- data012 |>
  mutate(score_Centralityonly = n_Deg) |>
  arrange(desc(score_Centralityonly)) |>
  mutate(rank_Centralityonly = row_number()) |>
  dplyr::select(UNIPROT, score_Centralityonly, rank_Centralityonly)

## ---- Condition: Equal-weight MCDA ----
# same 4 criteria as AHP-CDS (FC, FDR, Bet, Deg) but weights = 1/4 each
# instead of AHP-derived weights -> isolates the value AHP weighting adds
# on top of "just use all 4 criteria"
equal_mcda <- data012 |>
  mutate(score_EqualMCDA = 0.25*n_FC + 0.25*n_FDR + 0.25*n_Bet + 0.25*n_Deg) |>
  arrange(desc(score_EqualMCDA)) |>
  mutate(rank_EqualMCDA = row_number()) |>
  dplyr::select(UNIPROT, score_EqualMCDA, rank_EqualMCDA)

## ---- Condition: AHP-CDS ----
ahp_cds <- data012 |>
  distinct(UNIPROT, .keep_all = TRUE) |>
  arrange(desc(CDS)) |>
  mutate(rank_AHPCDS = row_number()) |>
  dplyr::select(UNIPROT, Symbol, CDS, rank_AHPCDS)

## ---- Condition: RWR (Ren et al. 2019) ----
de <- read.csv(here::here("PXD012162", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
de_map <- AnnotationDbi::select(org.Hs.eg.db, keys = de$UNIPROT, keytype = "UNIPROT", columns = "SYMBOL")
de_map <- de_map[!duplicated(de_map$UNIPROT), ]
de$SYMBOL = de_map[match(de$UNIPROT, de_map$UNIPROT), "SYMBOL"]
sig <- de[de$adj.P.Val < 0.05, ]
sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys = sig$UNIPROT, keytype = "UNIPROT", columns = "SYMBOL")
sym_map <- sym_map[!duplicated(sym_map$UNIPROT), ]

id_resp <- httr::POST(
  url = "https://string-db.org/api/tsv/get_string_ids",
  body = list(
    identifiers = paste(sym_map$SYMBOL, collapse = "\r"),
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
# ---- 1. Seed set: full Bonferroni-significant DEP set, restricted to giant component ----
sig_bonf <- de |>
  mutate(p_bonf = p.adjust(P.Value, method = "bonferroni")) |>
  filter(p_bonf < 0.05)

sig_bonf_graph <- sig_bonf |> filter(SYMBOL %in% igraph::V(g)$name)

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

out_path <- here::here("Ablation_Table", "results", "tables", "rwr_cross_ranks_full.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(ren2019_gr_df, out_path, row.names = FALSE)

cat(sprintf(
  "Saved GR ranking (%d proteins) to: %s\n",
  nrow(ren2019_gr_df), out_path
))

write.csv(sig_bonf_graph[, "SYMBOL", drop = FALSE], 
          here::here("Ablation_Table", "results", "tables", "ren2019_cross_seed_genes.csv"), 
          row.names = FALSE)


## ---- Merge all into one comparison table ----
ablation_ranks <- ahp_cds |>
  full_join(fc_only, by = "UNIPROT") |>
  full_join(centrality_only, by = "UNIPROT") |>
  full_join(equal_mcda, by = "UNIPROT") |>
  arrange(rank_AHPCDS)

## ---- Save ----
write.csv(ablation_ranks, here::here("Ablation_Table", "results", "tables", "ablation_cross_ranks_partial.csv"), row.names = FALSE)





