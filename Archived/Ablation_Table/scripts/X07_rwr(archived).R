suppressMessages({
  library(igraph)
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

edges <- read.delim(here::here("Ablation_Table", "data", "string_interactions366.tsv"), check.names = FALSE)
colnames(edges)[1] <- "node1"
edges_hi <- edges[edges$combined_score > 0.7, ]
g <- graph_from_data_frame(edges_hi[, c("node1", "node2", "combined_score")], directed = FALSE)

cat("Network: vcount =", vcount(g), " ecount =", ecount(g), "\n")

# restrict to the giant connected component -- RWR can't propagate across
# disconnected components, and the disconnected fragments here are tiny (2-3
# nodes each), not part of the main biological signal
comp <- components(g)
giant_nodes <- names(comp$membership[comp$membership == which.max(comp$csize)])
g <- induced_subgraph(g, giant_nodes)
cat("Giant component: vcount =", vcount(g), " ecount =", ecount(g), "\n")

# --- seed set: top 20 most FDR-significant limma genes present in this graph ---
de <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
sig <- de[de$adj.P.Val < 0.05, ]
sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys = sig$UNIPROT, keytype = "UNIPROT", columns = "SYMBOL")
sym_map <- sym_map[!duplicated(sym_map$UNIPROT), ]
sig <- merge(sig, sym_map, by = "UNIPROT")

sig_in_graph <- sig |> filter(SYMBOL %in% V(g)$name)
seed_genes <- sig_in_graph |> arrange(adj.P.Val) |> head(20) |> pull(SYMBOL)
cat("\nSeed set (top 20 most significant, present in giant component):\n")
print(seed_genes)
cat(
  "Seed genes retained:",
  length(seed_genes),
  "/",
  nrow(sig),
  "\n"
)
# --- RWR via igraph personalized PageRank ---
# p_(t+1) = (1-r) W p_t + r p0   <=>   PageRank with damping d = 1-r, personalized = p0
r <- 0.7
p0 <- setNames(rep(0, vcount(g)), V(g)$name)
cat(paste("sum p0 =", sum(p0)))
p0[seed_genes] <- 1 / length(seed_genes)
stopifnot(abs(sum(p0)-1)<1e-8)
rwr <- page_rank(g, damping = 1 - r, personalized = p0, weights = E(g)$combined_score)
rwr_scores <- rwr$vector

rwr_df <- data.frame(Symbol = names(rwr_scores), RWR_score = rwr_scores) |>
  arrange(desc(RWR_score)) |>
  mutate(rank_RWR = row_number())

write.csv(rwr_df, here::here("Ablation_Table", "results", "tables", "rwr_ranks_full.csv"), row.names = FALSE)

cat("\nTop 15 by RWR stationary probability:\n")
print(head(rwr_df, 15))

# --- compare against AHP-CDS on the overlap subset (evaluation only, not method design) ---
ahp <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"))
merged <- merge(ahp[, c("Symbol","rank_AHPCDS")], rwr_df[, c("Symbol","rank_RWR")], by = "Symbol")
cat("\nOverlap between AHP-CDS candidate pool and RWR-ranked giant component:", nrow(merged), "of 73\n")
cat("Spearman(AHP-CDS rank, RWR rank) on overlap:", cor(merged$rank_AHPCDS, merged$rank_RWR, method = "spearman"), "\n")

