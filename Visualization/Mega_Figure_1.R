library(ggplot2)
library(here)
library(tidyverse)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)

volcano_icon <- function(de_table,
                         logfc_col = "logFC",
                         padj_col  = "adj.P.Val",
                         fc_thresh  = 1,
                         fdr_thresh = 0.05,
                         point_size = 3) {
  
  de_table$category <- ifelse(
    de_table[[padj_col]] < fdr_thresh & de_table[[logfc_col]] > fc_thresh, "Upregulation",
    ifelse(
      de_table[[padj_col]] < fdr_thresh & de_table[[logfc_col]] < -fc_thresh, "Downregulation",
      "Non-significant"
    )
  )
  
  ggplot(de_table, aes(x = .data[[logfc_col]], 
                       y = -log10(.data[[padj_col]]), 
                       color = category)) +
    geom_point(size = point_size, alpha = 0.75, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    scale_color_manual(values = c(
      "Upregulation"    = "#0000FF80",
      "Downregulation"  = "#FF000080",
      "Non-significant" = "grey75"
    )) +
    theme_void(base_family = "Arial") +
    theme(legend.position = "none",
          plot.margin   = margin(1, 1, 1, 1))
}

de_table <- read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"))
volcano_plot_icon <- volcano_icon(de_table) 

volcano_plot_icon


#######

gsea_icon <- function(pathway, stats,
                      curve_color = "#0000FF",
                      tick_color  = "grey50",
                      smooth_n    = 200,
                      linewidth   = 1.2) {
  
  p <- plotEnrichment(pathway, stats)
  
  geom_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  line_idx  <- which(geom_classes == "GeomLine")
  seg_idx   <- which(geom_classes == "GeomSegment")
  hline_idx <- which(geom_classes == "GeomHline")
  
  built     <- ggplot_build(p)
  line_data <- built$data[[line_idx]]
  
  smoothed  <- as.data.frame(spline(line_data$x, line_data$y, n = smooth_n))
  names(smoothed) <- c("x", "y")
  
  p$layers[[line_idx]]$data    <- smoothed
  p$layers[[line_idx]]$mapping <- aes(x = x, y = y)
  
  p$layers[[line_idx]]$aes_params$colour    <- curve_color
  p$layers[[line_idx]]$aes_params$linewidth <- linewidth
  p$layers[[seg_idx]]$aes_params$colour     <- tick_color
  p$layers[hline_idx] <- NULL
  
  p +
    theme_void(base_family = "Arial") +
    theme(legend.position = "none",
          plot.margin = margin(1, 1, 1, 1))
}

gsea_df <- readRDS(here::here("PXD056161", "data", "reactome_gsea_object.rds"))
uniprot_to_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                           keys = de_table$UNIPROT,   
                                           keytype = "UNIPROT",
                                           columns = "ENTREZID")
uniprot_to_entrez <- uniprot_to_entrez[!duplicated(uniprot_to_entrez$UNIPROT), ]
de_table_entrez <- merge(de_table, uniprot_to_entrez,
                         by.x = "UNIPROT", by.y = "UNIPROT")


gene_sets <- gsea_df@geneSets

ranked_stats <- de_table$t
names(ranked_stats) <- de_table_entrez$ENTREZID
ranked_stats <- sort(ranked_stats, decreasing = TRUE)

pathway_genes <- gene_sets[["R-HSA-1474244"]]
  
  
gsea_plot_icon <- gsea_icon(pathway_genes, ranked_stats)
gsea_plot_icon


#________________________________STRING ___________________________________ #
library(igraph)
library(ggraph)
library(ggplot2)

network_icon <- function(n_nodes  = 26,
                         n_blocks = 4,
                         node_color = "#0000FF",
                         edge_color = "grey20",
                         seed = 30082026) {
  set.seed(seed)
  
  block_sizes <- as.vector(table(sample(seq_len(n_blocks), n_nodes, replace = TRUE)))
  block_sizes[block_sizes == 0] <- 1
  
  within_p  <- runif(1, 0.4, 0.7)
  between_p <- runif(1, 0.02, 0.1)
  
  pref_matrix <- matrix(between_p, n_blocks, n_blocks)
  diag(pref_matrix) <- within_p
  
  g <- sample_sbm(n_nodes, pref.matrix = pref_matrix, block.sizes = block_sizes)
  
  ggraph(g, layout = "stress") +
    geom_edge_link(color = edge_color, width = 0.4, alpha = 0.7) +
    geom_node_point(color = node_color, size = 4.2) +
    theme_void(base_family = "Arial") +
    theme(legend.position = "none",
          plot.margin = margin(1, 1, 1, 1))
}

string_plot_icon <- network_icon()
string_plot_icon

#________________________________AHP Overview ___________________________________ #
driver_landscape_overview <- function(candidates,
                                      logfc_col  = "logFC",
                                      degree_col = "degree",
                                      cds_col    = "CDS",
                                      gene_col   = "Symbol",
                                      top_n      = 20) {
  
  candidates$log_degree <- log1p(candidates[[degree_col]])
  candidates$top20 <- rank(-candidates[[cds_col]]) <= top_n
  
  ggplot(candidates, aes(x = .data[[logfc_col]],
                         y = log_degree,
                         color = .data[[cds_col]])) +
    geom_point(data = subset(candidates, !top20),
               size = 1.5, alpha = 0.5) +
    geom_point(data = subset(candidates, top20),
               size = 2.5, alpha = 0.9) +
    geom_text_repel(data = subset(candidates, top20),
                    aes(label = .data[[gene_col]]),
                    size = 2.8, family = "Arial",
                    max.overlaps = 20,
                    segment.size = 0.3, segment.color = "grey50") +
    scale_color_gradient(low = "#FFD9D9", high = "#8B0000", name = "CDS",
                         breaks = c(0.25, 0.50, 0.75),
                         guide = guide_colorbar(barwidth = 8, barheight = 0.6)) +
    labs(x = "log2 Fold-Change", y = "log1p(Degree)") +
    theme_paper()
}

scatter_plot <- driver_landscape_overview(CDS_df)
scatter_plot
