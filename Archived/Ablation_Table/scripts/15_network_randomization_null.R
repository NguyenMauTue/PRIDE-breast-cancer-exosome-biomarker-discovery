suppressMessages({ library(igraph); library(dplyr); library(openxlsx) })

set.seed(20260806)

## ---- Build the real network (same one behind network_summary.csv / BiomarkerCandidates_themed.csv) ----
edges <- read.delim(here::here("PXD056161", "data", "string_interactions.tsv"), check.names = FALSE)
colnames(edges)[1] <- "node1"
edges_hi <- edges[edges$combined_score > 0.7, ]
g_real <- graph_from_data_frame(edges_hi[, c("node1", "node2")], directed = FALSE)
g_real <- simplify(g_real)
cat("Real network: vcount =", vcount(g_real), " ecount =", ecount(g_real), "\n")

## ---- Real DE + candidate pool (TRUE labels, unchanged) ----
de <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1)
candidates_list <- read.csv(here::here("PXD056161", "results", "tables", "BiomarkerCandidates_themed.csv"))
pool <- candidates_list |> distinct(UNIPROT, .keep_all = TRUE) |>
  dplyr::select(UNIPROT, Symbol, logFC, adj.P.Val)

## ---- AHP weights (identical to script 12 / null_permutations.R) ----
ahp_weights <- function(M) { norm <- sweep(M, 2, colSums(M), "/"); rowMeans(norm) }
Pairwise.mat <- matrix(c(
  1,    4,   5,   2.5,
  1/4,  1,   2,   1,
  1/5,  1/2, 1,   1/5,
  1/2.5,1,   5,   1
), nrow = 4, byrow = TRUE, dimnames = list(c("FC","FDR","Bet","Deg"), c("FC","FDR","Bet","Deg")))
w <- ahp_weights(Pairwise.mat)
cat("AHP weights:", round(w, 4), "\n\n")

norm01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

## ---- 026 reference for cross-dataset r ----
module_tables_pxd012162 <- read.xlsx(here::here("PXD012162", "results", "tables", "Module_tables.xlsx"), sheet = "All")
n_pxd012162 <- nrow(module_tables_pxd012162)
module_tables_pxd012162 <- module_tables_pxd012162 |> arrange(desc(CDS)) |> mutate(fracrank_pxd012162 = row_number() / n_pxd012162) |>
  dplyr::select(UNIPROT, fracrank_pxd012162)

## ---- Given a graph (real or rewired), compute CDS + cross-dataset r ----
cds_and_crossr_for_graph <- function(g) {
  deg <- igraph::degree(g)
  bet <- igraph::betweenness(g)
  net_df <- data.frame(Symbol = names(deg), degree = deg, betweenness = bet[names(deg)])

  d <- pool |> inner_join(net_df, by = "Symbol")
  n <- nrow(d)
  d <- d |> mutate(
    n_FC  = norm01(abs(logFC)),
    n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)),
    n_Deg = norm01(log1p(degree)),
    CDS   = w["FC"]*n_FC + w["FDR"]*n_FDR + w["Bet"]*n_Bet + w["Deg"]*n_Deg
  ) |>
  arrange(desc(CDS)) |> mutate(frac_rank = row_number() / n)

  m <- inner_join(d[, c("UNIPROT","frac_rank")], module_tables_pxd012162, by = "UNIPROT")
  list(n_pool = n, r = cor(m$frac_rank, m$fracrank_pxd012162, method = "pearson"))
}

## ---- Real network result ----
real_result <- cds_and_crossr_for_graph(g_real)
cat("REAL network: n_pool =", real_result$n_pool, " cross-dataset r =", round(real_result$r, 4), "\n\n")

## ---- Degree-preserving rewiring null (Monte Carlo, since not exactly enumerable) ----
N_REWIRE <- 199
cat("Running", N_REWIRE, "degree-preserving network rewirings...\n")
null_r <- numeric(N_REWIRE)
for (i in 1:N_REWIRE) {
  g_null <- rewire(g_real, with = keeping_degseq(loops = FALSE, niter = ecount(g_real) * 10))
  res <- tryCatch(cds_and_crossr_for_graph(g_null), error = function(e) list(r = NA))
  null_r[i] <- res$r
  if (i %% 50 == 0) cat("  ", i, "/", N_REWIRE, "done\n")
}
null_r <- null_r[!is.na(null_r)]

cat("\n=== Network-randomization null distribution (n=", length(null_r), " rewirings) ===\n", sep="")
cat("Mean:  ", round(mean(null_r), 4), "\n")
cat("SD:    ", round(sd(null_r), 4), "\n")
cat("Median:", round(median(null_r), 4), "\n")
cat("Min:   ", round(min(null_r), 4), "\n")
cat("Max:   ", round(max(null_r), 4), "\n\n")

emp_p <- (sum(null_r >= real_result$r) + 1) / (length(null_r) + 1)
cat("=== Empirical p-value (real network vs", length(null_r), "degree-preserving rewirings) ===\n")
cat("Real network cross-dataset r:", round(real_result$r, 4), "\n")
cat("Rank among (null + real):", sum(c(null_r, real_result$r) >= real_result$r), "of", length(null_r)+1, "\n")
cat("Empirical p-value:", round(emp_p, 4), "\n")

out <- data.frame(rewire_id = 1:length(null_r), cross_dataset_r = null_r)
out <- rbind(out, data.frame(rewire_id = 0, cross_dataset_r = real_result$r))
write.csv(out, here::here("Ablation_Table", "results", "tables", "network_randomization_null.csv"), row.names = FALSE)
cat("\nSaved: results/tables/network_randomization_null.csv (rewire_id=0 is the REAL network)\n")
