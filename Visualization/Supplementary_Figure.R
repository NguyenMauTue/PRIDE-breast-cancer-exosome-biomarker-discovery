###################################################################
###################################################################
##### 
#####Visualization for supplementary
#####
#####
###################################################################


# ------------------- Technical replication heatmap ------------------------------
library(here)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
source(here::here("R", "Helper", "theme_paper.R"))
library(reshape2)
library(here)
library(scales)
library(ComplexHeatmap)

mat <- as.matrix(read.csv(
  here::here("PXD056161", "results", "tables", "technical_correlation_matrix.csv"),
  row.names = 1, check.names = FALSE
))

df <- melt(mat, varnames = c("Sample1", "Sample2"), value.name = "r")
r_range <- range(df$r)   

ggplot(df, aes(Sample2, Sample1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.3, color = "grey20") +
  scale_fill_gradient(
    low = "white", high = "#0072B2",
    limits = r_range,
    breaks = pretty_breaks(n = 5)(r_range),
    name = "Pearson r"
  ) +
  coord_fixed() +
  scale_y_discrete(limits = rev(levels(factor(df$Sample1)))) +
  theme_paper(base_size = 16) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.6, "cm"),
    legend.key.width  = unit(0.4, "cm")
  ) +
  guides(fill = guide_colorbar(title.position = "top"))


# ------------------- MNRA sigmoid curve ------------------------------
# Convert LFQ matrix to long format
log_lfq_matrix <- readRDS(here::here("PXD056161", "data", "log_lfq_matrix.rds"))
long_df <- melt(log_lfq_matrix)
colnames(long_df) <- c("Protein", "Sample", "Intensity")

# Detection status
long_df$Detected <- !is.na(long_df$Intensity)

# Mean abundance per protein
protein_mean <- rowMeans(log_lfq_matrix, na.rm = TRUE)
long_df$MeanIntensity <- protein_mean[long_df$Protein]

# Visualization of detection bias
# Convert logical detection to numeric
long_df$Detected_Number = 0
long_df$Detected_Number[long_df$Detected == TRUE] = 1
fit <- glm(Detected_Number ~ MeanIntensity, data = long_df, family = "binomial")
label_x <- 24
label_y <- predict(fit, newdata = data.frame(MeanIntensity = label_x), type = "response")

pMissing <- ggplot(long_df, aes(x = MeanIntensity, y = Detected_Number)) +
  geom_jitter(height = 0.05, width = 0.1, alpha = 0.15,
              color = "grey30", size = 0.8) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE,
              color = pal_sig["Up"],   
              linewidth = 1) +
  annotate("text", x = label_x, y = label_y - 0.01, label = "fitted P(detected)",
           color = pal_sig["Up"], size = 3, hjust = 0, family = "Arial") +
  labs(x = "Mean intensity", y = "Detected (0/1)") +
  theme_paper()

pMissing

# ------------------- Missingness boxplot -------------------------------
# ── Helper: matrix -> long dataframe ────────────────────────────────────
PXD056 <- readRDS(here::here("PXD056161", "data", "collapsed_matrix.rds"))

mat_to_long <- function(mat, dataset_label) {
  as.data.frame(mat) %>%
    pivot_longer(everything(), names_to = "Sample", values_to = "Intensity") %>%
    filter(!is.na(Intensity)) %>%
    mutate(Dataset = dataset_label)
}

# ── PXD056161: assign group ─────────────────────────────────────────────
df056 <- mat_to_long(PXD056, "PXD056161") %>%
  mutate(Group = case_when(
    grepl("Normal", Sample, ignore.case = TRUE) ~ "MCF10A (Normal)",
    grepl("Tumor",  Sample, ignore.case = TRUE) ~ "MDA-MB-231 (Tumor)",
    TRUE ~ Sample
  ))

# ── CVs (two different variables now — the original overwrote cv_056
# on the second computation, silently losing the mean/SD version) ──────
sample_means_056   <- colMeans(PXD056, na.rm = TRUE)
cv_mean_056        <- sd(sample_means_056) / mean(sample_means_056) * 100

sample_medians_056 <- apply(PXD056, 2, median, na.rm = TRUE)
cv_mad_056         <- mad(sample_medians_056) / median(sample_medians_056) * 100

cat("CV PXD056 (mean/SD):", round(cv_mean_056, 2), "%\n")
cat("CV PXD056 (median/MAD):", round(cv_mad_056, 2), "%\n")

# ── Plot ─────────────────────────────────────────────────────────────────
# theme_paper() already sets legend.position="bottom", so no need for
# the separate legend_plot + cowplot::get_legend() dance from before.
missingBox <- ggplot(df056, aes(x = Sample, y = Intensity, fill = Group)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = pal_condition, name = NULL) +
  scale_x_discrete(limits = c("Normal1", "Normal2", "Normal3",
                              "Tumor1", "Tumor2", "Tumor3")) +
  labs(
    x = NULL,
    y = expression(log[2] ~ "intensity"),
    caption = sprintf("CV (mean/SD) = %.2f%%   |   CV (median/MAD) = %.2f%%",
                      cv_mean_056, cv_mad_056)
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(size = 7, color = "grey30", hjust = 1)
  )
missingBox



# ------------------- Missingness heatmap (protein x sample) -------------
# ── Build status_mat from your raw (pre-filter, pre-imputation) matrix ──
# raw_mat: proteins x samples, NA = not detected. Adjust group regex to
# match your actual sample naming.
build_status_matrix <- function(raw_mat, group_regex_A = "Normal", group_regex_B = "Tumor") {
  grpA_cols <- grep(group_regex_A, colnames(raw_mat), ignore.case = TRUE)
  grpB_cols <- grep(group_regex_B, colnames(raw_mat), ignore.case = TRUE)
  
  frac_missing_A <- rowMeans(is.na(raw_mat[, grpA_cols, drop = FALSE]))
  frac_missing_B <- rowMeans(is.na(raw_mat[, grpB_cols, drop = FALSE]))
  
  status_mat <- matrix("Present", nrow = nrow(raw_mat), ncol = ncol(raw_mat),
                       dimnames = dimnames(raw_mat))
  status_mat[is.na(raw_mat)] <- "Missing"   
  
  status_mat
}
#     "Present"   — protein detected in this sample
#     "Missing"   — missing  


raw_mat <- readRDS(here::here("PXD056161", "data", "log_lfq_matrix.rds"))

status_mat <- build_status_matrix(raw_mat)

# ── Sample -> condition annotation (reuses pal_condition, same mapping
# as the boxplot panel, so the two panels read as one coherent story) ───
sample_condition <- case_when(
  grepl("MCF10", colnames(status_mat), ignore.case = TRUE) ~ "MCF10A (Normal)",
  grepl("MDA",   colnames(status_mat), ignore.case = TRUE) ~ "MDA-MB-231 (Tumor)",
  TRUE ~ NA_character_
)

col_anno <- HeatmapAnnotation(
  Condition = sample_condition,
  col = list(Condition = pal_condition),
  annotation_name_side = "left",
  simple_anno_size = unit(2, "mm")
)

binary_mat <- ifelse(status_mat == "Present", 0, 1)  # Imputed & Excluded both = "missing" for clustering
row_dend <- hclust(dist(binary_mat, method = "binary"), method = "ward.D2")
col_dend <- hclust(dist(t(binary_mat), method = "binary"), method = "ward.D2")
# ── Status colors — reuse the vocabulary already established elsewhere
# rather than inventing a 4th ad hoc palette: Present = neutral grey
# (baseline, nothing to flag), Imputed = pal_sig["Down"] (recoverable,
# not lost).
status_colors <- c(
  "Present"  = unname(pal_sig["NS"]),
  "Missing"  = "#3D3D3D"
)

# ── Plot ──────────────────────────────────────────────────────────────────
ht <- Heatmap(
  status_mat,
  name = "Status",
  col = status_colors,
  top_annotation = col_anno,
  show_row_names = FALSE,          
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  cluster_rows = row_dend,             
  cluster_columns = col_dend,          
  row_dend_width = unit(4, "mm"),      
  column_dend_height = unit(4, "mm"),  
  heatmap_legend_param = list(         
    title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
  )
)

draw(ht, merge_legend = TRUE)



# -------------------   PCA imputation stability test ------------------
collapsed_matrix = readRDS(
  here::here("PXD056161", "data", "collapsed_matrix.rds")
)

imputed_matrix = readRDS(
  here::here("PXD056161", "data", "imputed_matrix.rds")
)

complete_inx = rowSums(is.na(collapsed_matrix)) == 0

pca_before = prcomp(
  t(collapsed_matrix[complete_inx,]),
  scale = TRUE
)

pca_after = prcomp(
  t(imputed_matrix),
  scale = TRUE
)


make_pca_plot <- function(pca_obj, title) {
  df <- as.data.frame(pca_obj$x)
  percent_var <- round(100 * pca_obj$sdev^2 / sum(pca_obj$sdev^2), 1)
  
  df$sample <- rownames(df)
  df$group <- factor(
    ifelse(grepl("Normal", df$sample), "MCF10A (Normal)", "MDA-MB-231 (Tumor)"),
    levels = c("MCF10A (Normal)", "MDA-MB-231 (Tumor)")
  )
  
  ggplot(df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
    scale_color_manual(
      values = c("MCF10A (Normal)"    = pal_sig[["Down"]],
                 "MDA-MB-231 (Tumor)" = pal_sig[["Up"]]),
      drop = FALSE
    ) +
    labs(
      title = title,
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Cell line"
    ) +
    theme_paper()
}

p_before <- make_pca_plot(pca_before, "Before imputation")
p_after  <- make_pca_plot(pca_after,  "After imputation")

fig_supp <- (p_before + p_after) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig_supp


# ============================================================
# contaminant_filter_breakdown.R
# Standalone script — does NOT modify 09_extract_genes_for_STRING.R
# (frozen, manuscript already references it). Reads the same raw
# inputs to reconstruct the pre-filter checkpoint for visualization
# purposes only — not intended to reproduce the official reported
# numbers unless the sanity check below passes.
# ============================================================

library(dplyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggVennDiagram)

# ---- Step 1: recover the original ~193 core-enrichment genes (same as script 09) ----
gsea_reactome_res <- read.csv(
  here::here("PXD056161", "results", "tables", "reactome_gsea_filtered.csv")
)
genes <- unique(unlist(strsplit(gsea_reactome_res$core_enrichment, "/")))

# ---- Step 2: Entrez -> UniProt mapping, using org.Hs.eg.db instead of biomaRt ----
# NOTE: the original script used biomaRt::getBM() via an Ensembl mirror.
# That mirror is no longer usable, so this step falls back to org.Hs.eg.db
# (offline lookup). The two ID-mapping sources can disagree on a handful
# of edge-case genes (multi-mapping IDs, withdrawn symbols, etc.), so
# anything downstream here is APPROXIMATE and for illustrative Venn
# purposes only — not a substitute for the published 41/152 figures
# unless verified against the official output.
results <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(genes),
  keytype = "ENTREZID",
  columns = c("SYMBOL", "UNIPROT")
) %>%
  dplyr::rename(entrezgene_id = ENTREZID,
                hgnc_symbol   = SYMBOL,
                uniprotswissprot = UNIPROT) %>%
  dplyr::filter(!is.na(uniprotswissprot) & uniprotswissprot != "") %>%
  dplyr::distinct()

# ---- Step 3: reuse annotation_raw.csv unchanged ----
annotations_all <- read.csv(here::here("PXD056161", "data", "annotation_raw.csv"))
annotations <- annotations_all %>%
  filter(uniprotswissprot %in% results$uniprotswissprot)

annotations_grouped <- annotations %>%
  group_by(uniprotswissprot) %>%
  summarise(
    Symbol_biomat = paste(unique(external_gene_name), collapse = "; "),
    Protein_Description = paste(unique(description), collapse = "; "),
    .groups = "drop"
  )

Symbol_annotated_raw <- merge(
  results, annotations_grouped,
  by = "uniprotswissprot", all.x = TRUE
)

# ---- Step 4: re-apply the exact filter logic from script 09, but keep the removed rows ----
contaminant_pattern <- "histone|keratin|actin|tubulin"
manual_blacklist <- c("FMNL1", "H2AX", "H4C6", "H3C1", "H3-3B", "H2AZ2")

regex_flag <- with(Symbol_annotated_raw,
                   grepl(contaminant_pattern, Protein_Description, ignore.case = TRUE) |
                     grepl(contaminant_pattern, as.character(Symbol_biomat), ignore.case = TRUE)
)
blacklist_flag <- Symbol_annotated_raw$Symbol_biomat %in% manual_blacklist

after_contam <- Symbol_annotated_raw[!(regex_flag | blacklist_flag), ]
na_flag <- !complete.cases(after_contam)

# ---- Step 5: MANDATORY sanity check before trusting any number below ----
# Compare against the file already produced by the original (biomaRt-based) script
official <- read.csv(here::here("PXD056161", "results", "tables", "annotated_gene_pool.csv"))
n_reconstructed <- nrow(after_contam[!na_flag, ])
n_official <- nrow(official)

cat("Reconstructed row count: ", n_reconstructed, "\n")
cat("Official row count (annotated_gene_pool.csv): ", n_official, "\n")
if (n_reconstructed != n_official) {
  warning("MISMATCH between org.Hs.eg.db mapping and the original biomaRt mapping — ",
          "the breakdown/Venn below is illustrative only and should NOT be ",
          "reported as the official manuscript figure.")
}

# ---- Step 6: split out the removed groups, build the Venn ----
removed_regex     <- Symbol_annotated_raw[regex_flag, ]
removed_blacklist <- Symbol_annotated_raw[blacklist_flag & !regex_flag, ]
removed_naomit    <- after_contam[na_flag, ]
removed_contam    <- bind_rows(removed_regex, removed_blacklist)

sets <- list(
  Histone = removed_contam %>%
    filter(str_detect(str_to_lower(Protein_Description), "histone") |
             str_detect(str_to_lower(Symbol_biomat), "histone") |
             Symbol_biomat %in% c("H2AX","H4C6","H3C1","H3-3B","H2AZ2")) %>%
    pull(Symbol_biomat),
  Keratin = removed_contam %>%
    filter(str_detect(str_to_lower(Protein_Description), "keratin") |
             str_detect(str_to_lower(Symbol_biomat), "keratin")) %>%
    pull(Symbol_biomat),
  Actin = removed_contam %>%
    filter(str_detect(str_to_lower(Protein_Description), "actin") |
             str_detect(str_to_lower(Symbol_biomat), "actin") |
             Symbol_biomat == "FMNL1") %>%
    pull(Symbol_biomat),
  Tubulin = removed_contam %>%
    filter(str_detect(str_to_lower(Protein_Description), "tubulin") |
             str_detect(str_to_lower(Symbol_biomat), "tubulin")) %>%
    pull(Symbol_biomat)
)

ggVennDiagram(sets, label = "count") +
  scale_fill_gradient(low = "white", high = "#0072B2") +
  theme_paper(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.grid = element_blank()
  )

# ---------------------------- NES Distribution hist gsea_unfiltered ---------------------------
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
# ---- Steps 1-4: identical to the original script, unchanged ----
set.seed(30032026)
imputed_result <- read.csv(
  here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"),
  row.names = 1
)

protein_ids <- rownames(imputed_result)
uniprot_major <- str_split(protein_ids, ";") |> sapply(`[`, 1) |> str_trim()
imputed_result$UNIPROT <- uniprot_major

gene_map <- bitr(
  unique(imputed_result$UNIPROT),
  fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db
)
rank_df <- merge(imputed_result, gene_map, by = "UNIPROT")
rank_df <- rank_df[order(abs(rank_df$t), decreasing = TRUE), ]
rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]

gene_list <- rank_df$t
names(gene_list) <- as.character(rank_df$ENTREZID)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!is.na(gene_list)]

# ---- Step 5: SAME call, only pvalueCutoff changed 0.05 -> 1 ----
# pvalueCutoff only controls which rows the function RETURNS after
# the permutation test — it does not change the permutation itself,
# so NES/p-value for pathways that were already significant should
# come out identical to the original run (checked below).
gsea_reactome_unfiltered <- gsePathway(
  geneList     = gene_list,
  organism     = "human",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 1,       # <-- only change: keep every tested pathway
  verbose      = FALSE,
  seed         = TRUE
)

gsea_full_res <- as.data.frame(gsea_reactome_unfiltered)

# ---- Sanity check: re-applying the original filter should reproduce
# the official reactome_gsea_filtered.csv exactly ----
official_filtered <- read.csv(
  here::here("PXD056161", "results", "tables", "reactome_gsea_filtered.csv")
)
reconstructed_filtered <- gsea_full_res[
  abs(gsea_full_res$NES) > 1.5 & gsea_full_res$p.adjust < 0.05, 
]

cat("Official filtered pathway count:      ", nrow(official_filtered), "\n")
cat("Reconstructed filtered pathway count: ", nrow(reconstructed_filtered), "\n")
if (nrow(official_filtered) != nrow(reconstructed_filtered)) {
  warning("MISMATCH — the unfiltered rerun does not reproduce the official ",
          "filtered result exactly. Do not use this histogram for the ",
          "manuscript until resolved.")
}

# ---- Histogram: NES distribution, significant vs. not ----
gsea_full_res <- gsea_full_res %>%
  mutate(sig_status = ifelse(p.adjust < 0.05, "Significant", "Not significant"))

ggplot(gsea_full_res, aes(x = NES, fill = sig_status)) +
  geom_histogram(binwidth = 0.25, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", 
             color = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", 
             color = "grey60", linewidth = 0.4) +
  scale_fill_manual(values = c("Significant" = pal_sig[["Up"]], "Not significant" = "grey70")) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "Number of pathways", 
    fill = "FDR < 0.05",
    caption = "Dashed lines mark |NES| = 1.5, the effect-size cutoff used\nalongside FDR < 0.05 for the final candidate pathway list."
  ) +
  theme_paper(base_size = 16) +
  theme(legend.position = "bottom")
