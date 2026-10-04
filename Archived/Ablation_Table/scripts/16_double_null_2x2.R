suppressMessages({ library(limma); library(igraph); library(dplyr); library(openxlsx) })

set.seed(20260806)

## ---- Real network ----
edges <- read.delim(here::here("PXD056161", "data", "string_interactions.tsv"), check.names = FALSE)
colnames(edges)[1] <- "node1"
edges_hi <- edges[edges$combined_score > 0.7, ]
g_real <- simplify(graph_from_data_frame(edges_hi[, c("node1", "node2")], directed = FALSE))

## ---- Candidate pool + imputed matrix ----
imputed_matrix <- readRDS(here::here("PXD056161", "data", "imputed_matrix.rds"))
imputed_matrix <- imputed_matrix[, c("Normal1","Normal2","Normal3","Tumor1","Tumor2","Tumor3")]
candidates_list <- read.csv(here::here("PXD056161", "results", "tables", "BiomarkerCandidates_themed.csv"))
pool_symbols <- candidates_list |> distinct(UNIPROT, .keep_all = TRUE) |> dplyr::select(UNIPROT, Symbol)

## ---- AHP weights ----
ahp_weights <- function(M) { norm <- sweep(M, 2, colSums(M), "/"); rowMeans(norm) }
Pairwise.mat <- matrix(c(
  1,    4,   5,   2.5,
  1/4,  1,   2,   1,
  1/5,  1/2, 1,   1/5,
  1/2.5,1,   5,   1
), nrow = 4, byrow = TRUE, dimnames = list(c("FC","FDR","Bet","Deg"), c("FC","FDR","Bet","Deg")))
w <- ahp_weights(Pairwise.mat)
norm01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

## ---- 026 reference ----
module_tables_pxd012162 <- read.xlsx(here::here("PXD012162", "results", "tables", "Module_tables.xlsx"), sheet = "All")
n_pxd012162 <- nrow(module_tables_pxd012162)
module_tables_pxd012162 <- module_tables_pxd012162 |> arrange(desc(CDS)) |> mutate(fracrank_pxd012162 = row_number() / n_pxd012162) |>
  dplyr::select(UNIPROT, fracrank_pxd012162)

## ---- 9 null label splits (identical enumeration to script 05) ----
idx <- 1:6
combos <- combn(6, 3)
splits <- list(); seen <- c()
for (i in 1:ncol(combos)) {
  A <- sort(combos[, i]); B <- sort(setdiff(idx, A))
  if (!(1 %in% A)) { tmp <- A; A <- B; B <- tmp }
  key <- paste(A, collapse = ",")
  if (!(key %in% seen)) { seen <- c(seen, key); splits[[length(splits) + 1]] <- list(A = A, B = B) }
}
null_splits <- splits[-1]   # drop TRUE split, keep 9 null label splits

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
  res |> mutate(UNIPROT = sapply(strsplit(X, ";"), `[`, 1))
}

## ---- Pre-compute DE for each of the 9 null label splits, restricted to the pool ----
cat("Pre-computing DE for 9 null label splits...\n")
null_de_list <- lapply(null_splits, function(sp) {
  de <- run_de_for_split(sp$A, sp$B)
  de |> inner_join(pool_symbols, by = "UNIPROT") |> distinct(UNIPROT, .keep_all = TRUE)
})

## ---- Given DE data + a graph, compute cross-dataset r ----
cross_r_given <- function(de_df, g) {
  deg <- igraph::degree(g); bet <- igraph::betweenness(g)
  net_df <- data.frame(Symbol = names(deg), degree = deg, betweenness = bet[names(deg)])
  d <- de_df |> inner_join(net_df, by = "Symbol")
  n <- nrow(d)
  d <- d |> mutate(
    n_FC = norm01(abs(logFC)), n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)), n_Deg = norm01(log1p(degree)),
    CDS = w["FC"]*n_FC + w["FDR"]*n_FDR + w["Bet"]*n_Bet + w["Deg"]*n_Deg
  ) |> arrange(desc(CDS)) |> mutate(frac_rank = row_number()/n)
  m <- inner_join(d[, c("UNIPROT","frac_rank")], module_tables_pxd012162, by = "UNIPROT")
  cor(m$frac_rank, m$fracrank_pxd012162, method = "pearson")
}

## ---- Double null: 9 label splits x 20 network rewirings each = 180 draws ----
N_REWIRE_PER_SPLIT <- 20
cat("Running double null:", length(null_de_list), "label splits x", N_REWIRE_PER_SPLIT, "rewirings =",
    length(null_de_list)*N_REWIRE_PER_SPLIT, "draws...\n")

double_null_r <- c()
for (s in seq_along(null_de_list)) {
  for (k in 1:N_REWIRE_PER_SPLIT) {
    g_null <- rewire(g_real, with = keeping_degseq(loops = FALSE, niter = ecount(g_real) * 10))
    r <- tryCatch(
      cross_r_given(null_de_list[[s]], g_null),
      error = function(e) { cat("ERROR at s=", s, ":", conditionMessage(e), "\n"); NA }
    )
    double_null_r <- c(double_null_r, r)
  }
  cat("  label split", s, "/9 done\n")
}
double_null_r <- double_null_r[!is.na(double_null_r)]

cat("\n=== Double null (label random x network random), n=", length(double_null_r), " ===\n", sep="")
cat("Mean:  ", round(mean(double_null_r), 4), "\n")
cat("SD:    ", round(sd(double_null_r), 4), "\n")
cat("Median:", round(median(double_null_r), 4), "\n")
cat("Min:   ", round(min(double_null_r), 4), " Max:", round(max(double_null_r), 4), "\n\n")

## ---- Real (label true, network true) result for comparison ----
real_de <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"), row.names = 1) |>
  inner_join(pool_symbols, by = "UNIPROT") |> distinct(UNIPROT, .keep_all = TRUE)
real_r <- cross_r_given(real_de, g_real)
cat("REAL (label true x network true): cross-dataset r =", round(real_r, 4), "\n\n")

emp_p <- (sum(double_null_r >= real_r) + 1) / (length(double_null_r) + 1)
cat("=== Empirical p-value: real vs double-null (label random x network random) ===\n")
cat("Rank among (double_null + real):", sum(c(double_null_r, real_r) >= real_r), "of", length(double_null_r)+1, "\n")
cat("Empirical p-value:", round(emp_p, 4), "\n\n")

## ---- Full 2x2 summary table ----
label_true_net_true   <- real_r
label_true_net_null   <- read.csv(here::here("Ablation_Table", "results", "tables", "network_randomization_null.csv"))
label_true_net_null_mean <- mean(label_true_net_null$cross_dataset_r[label_true_net_null$rewire_id != 0])
label_null_net_true   <- read.csv(here::here("Ablation_Table", "results", "tables", "null_cross_dataset_r_distribution.csv"))
label_null_net_true_mean <- mean(label_null_net_true$cross_dataset_r[!label_null_net_true$is_true])

summary_2x2 <- data.frame(
  Network = c("TRUE", "TRUE", "RANDOM", "RANDOM"),
  Label   = c("TRUE", "RANDOM (mean of 9)", "TRUE (mean of 199)", "RANDOM (mean of 180)"),
  cross_dataset_r = c(label_true_net_true, label_null_net_true_mean, label_true_net_null_mean, mean(double_null_r))
)
cat("=== Full 2x2 factorial summary ===\n")
print(summary_2x2, row.names = FALSE)

write.csv(data.frame(draw_id = 1:length(double_null_r), cross_dataset_r = double_null_r),
          here::here("Ablation_Table", "results", "tables", "double_null_label_network.csv"), row.names = FALSE)
write.csv(summary_2x2, here::here("Ablation_Table", "results", "tables", "factorial_2x2_summary.csv"), row.names = FALSE)
cat("\nSaved: results/tables/double_null_label_network.csv, results/tables/factorial_2x2_summary.csv\n")
