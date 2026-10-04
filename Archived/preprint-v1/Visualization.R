library(openxlsx)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(grid)

library(tidygraph)
library(ggraph)
library(graphlayouts)

library(forcats)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggvenn)
library(UpSetR)
library(dplyr)

PXD012162 = read.csv("Cross checking/results/BiomarkerCandidates_themed.csv")
PXD056161 = read.csv("results/BiomarkerCandidates_themed.csv")
string_df = read.delim("data/string_interactions.tsv")
cat_colors = c("ECM / adhesion"="#E8789A",
"Signaling / immune"="#F9A825",
"Cytoskeleton"="#2E7D32",
"EV / trafficking"="#6A1B9A",
"RNA-binding / nuclear"="#C62828")


visualize_df = read.xlsx("cross_validation_workbook.xlsx", sheet = 5, startRow = 2)

p1 = ggplot(visualize_df, aes(x=logFC_056, y=logFC_012, color=`Concordant?`)) +
  geom_point(size=3, alpha=0.8) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_vline(xintercept=0, color="grey80") +
  geom_hline(yintercept=0, color="grey80") +
  geom_text_repel(aes(label=Symbol), size=3, max.overlaps=20) +
  scale_color_manual(values=c("Yes"="#2196F3","No"="#F44336")) +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5,
           label=paste0("Spearman r = ", round(cor(visualize_df$logFC_056, visualize_df$logFC_012, method="spearman"),2)),
           size=3.5, fontface="italic") +
  labs(x="logFC (PXD056161)", y="logFC (PXD012162)", color="Concordant") +
  theme_classic() +
  theme(text=element_text(size=12, ))

p2 = ggplot(visualize_df, aes(x=Rank_056, y=Rank_012)) +
  geom_point(aes(color=Category), size=3, alpha=0.8) +
  geom_text_repel(
    data=subset(visualize_df, Rank_056 <= 20),  # chỉ label top
    aes(label=Symbol), size=3
  ) +
  scale_color_manual(values= cat_colors) +
  labs(x="Rank in PXD056161 (of 61)", y="Rank in PXD012162 (of 426)") +
  theme_classic() +
  theme(text=element_text(size=12, ))


p12 = p1 + p2 + plot_annotation(tag_levels = 'A')
p12

ggsave("cross_val_AB.png", plot = p12, 
       width = 13, height = 8, dpi = 600, units = "in")

library(ComplexHeatmap)
library(circlize)

#Heatmap

# 1. Sắp xếp lại dataframe theo Category
# Bạn có thể sắp xếp thêm theo logFC bên trong mỗi Category để nó mượt hơn nữa
visualize_df <- visualize_df[order(visualize_df$Category, -visualize_df$logFC_056), ]

# 2. Tạo lại ma trận và vector Symbol sau khi đã sắp xếp
mat <- as.matrix(visualize_df[, c("logFC_056", "logFC_012")])
rownames(mat) <- visualize_df$Symbol

# 3. Vẽ Heatmap (với show_row_dend = FALSE và cluster_rows = FALSE)
hp_sorted <- Heatmap(mat, 
                     name = "logFC", 
                     col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
                     column_labels = c("PXD056161", "PXD012162"),
                     cluster_rows = FALSE, 
                     show_row_dend = FALSE,
                     
                     right_annotation = rowAnnotation(
                       Category = visualize_df$Category,
                       RankDiff = anno_barplot(
                         as.numeric(gsub("\\+", "", visualize_df$Rank.shift)), 
                         gp = gpar(fill = "#455A64", col = "#455A64"), 
                         border = FALSE,
                         width = unit(1.5, "cm")
                       ),
                       col = list(Category = cat_colors),
                       show_annotation_name = TRUE
                     ),
                     
                     row_names_gp = gpar(fontsize = 9, fontface = "italic"),
                     column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                     cluster_columns = FALSE,
                     rect_gp = gpar(col = "white", lwd = 1),
                     column_title = "Supplementary Heatmap (Grouped by Category)"
)

# Xuất kết quả
draw(hp_sorted, heatmap_legend_side = "left", annotation_legend_side = "right")
hp_grob <- grid.grabExpr(draw(hp_sorted, 
                              heatmap_legend_side = "left", 
                              annotation_legend_side = "right"))
png("cross_val_C.png", width = 13, height = 6, 
    units = "in", res = 300)
draw(hp_sorted, heatmap_legend_side = "left", 
     annotation_legend_side = "right")
dev.off()



graph_data <- as_tbl_graph(string_df, directed = FALSE) %>%
  activate(nodes) %>%
  left_join(visualize_df, by = c("name" = "Symbol")) %>%
  mutate(degree = centrality_degree()) %>%
  filter(!is.na(Category)) %>%
  filter(!node_is_isolated())

p3 = ggraph(graph_data, layout = "stress") + 
  # Vẽ cạnh mờ thôi để làm nền
  geom_edge_link(alpha = 0.45, color = "grey70") +
  
  # Node to nhỏ theo CDS và màu theo Category
  geom_node_point(aes(size = CDS_056, color = Category), alpha = 0.8) +
  scale_color_manual(values = cat_colors) +
  scale_size_continuous(range = c(3, 12)) + 
  
  # Hiện tên nhiều hơn một chút để đỡ trống trải
  geom_node_text(aes(label = ifelse(CDS_056 > 0.4, name, "")), 
                 repel = TRUE, 
                 size = 3.5, 
                 fontface = "bold.italic",
                 box.padding = 0.6,
                 ) +
  
  theme_graph() +
  theme(legend.position = "right") +
  labs(
    color = "Functional Category", 
    size = "CDS Score",
    caption = "Edges represent STRING interactions",
    title = "PPI Network of Candidate Proteins",
    x = NULL,
    y = NULL) +
  theme(
    legend.title = element_text(face = "bold", size = 10, ),
    legend.text = element_text(size = 9, ),
    plot.margin = margin(10, 10, 10, 10),
    text = element_text(),
    plot.title = element_text(, face = "bold", size = 16, hjust = 0.5)
    )
  

p3

ggsave("PPI Network plot.png", plot = p3, 
       width = 7, height = 7, dpi = 300, units = "in")


genes_to_test <- bitr(PXD056161$Symbol, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)$ENTREZID

ego_BP <- enrichGO(gene          = genes_to_test,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

#CC (Cellular Component)
ego_CC <- enrichGO(gene          = genes_to_test,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

#MF (Molecular Function)
ego_MF <- enrichGO(gene          = genes_to_test,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)



df_BP_top <- ego_BP %>% as.data.frame() %>% top_n(10, wt = -p.adjust) %>% mutate(Ontology = "Biological Process")
df_CC_top <- ego_CC %>% as.data.frame() %>% top_n(10, wt = -p.adjust) %>% mutate(Ontology = "Cellular Component")
df_MF_top <- ego_MF %>% as.data.frame() %>% top_n(10, wt = -p.adjust) %>% mutate(Ontology = "Molecular Function")

# Combine 
df_final_go <- bind_rows(df_BP_top, df_CC_top, df_MF_top)
names(df_final_go)

# Vẽ biểu đồ thanh 3 tầng
p4 = ggplot(df_final_go, aes(x = -log10(p.adjust), 
                        y = reorder(Description, -p.adjust) 
                        )) +
  geom_bar(stat = "identity", width = 0.8, aes(fill = Ontology)) +
  facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
  geom_text(aes(label = Count), hjust = -0.2, size = 3) +
  scale_fill_manual(values = c("Biological Process" = "#d62728", 
                               "Cellular Component" = "#2ca02c", 
                               "Molecular Function" = "#1f77b4")) +
  xlim(0, max(-log10(df_final_go$p.adjust)) + 2) + # Chừa chỗ cho nhãn Count
  labs(x = "-log10(Adjusted P-value)",
       y = NULL,
       fill = "Category",
       title = "Gene Ontology Enrichment Analysis") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        text = element_text(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = element_text(face = "bold")) 
  
p4

ggsave("results/Gene Ontology plot.png", plot = p4, 
       width = 10, height = 7, dpi = 600, units = "in")

df_sensitivity = PXD056161[ , c("Symbol", "rank_base", "rank_sd", "rank_min", "rank_max", "robustness_label")]
df_sensitivity = df_sensitivity %>%
  arrange(rank_base) %>%
  slice_head(n = 25) %>%  # top 25 theo rank gốc
  mutate(Symbol = fct_reorder(Symbol, -rank_base))


  
p5 = ggplot(df_sensitivity, aes(y = Symbol)) +
  # error bar = range
  geom_segment(aes(x = rank_min, xend = rank_max,
                   yend = Symbol,
                   color = robustness_label), 
               linewidth = 1.2, alpha = 0.6) +
  # dot = rank gốc
  geom_point(aes(x = rank_base, color = robustness_label), 
             size = 3) +
  scale_color_manual(
    values = c("robust_candidate" = "#1565C0",
               "weight_sensitive_candidate" = "#E65100"),
    labels = c("robust_candidate" = "Robust",
               "weight_sensitive_candidate" = "Weight-sensitive")
  ) +
  scale_x_continuous(name = "CDS Rank (lower = better)") +
  labs(y = NULL, color = "Robustness",
       title = "Rank stability under AHP weight perturbation (±20%)") +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text())

p5

ggsave("results/AHP Robustness plot.png", plot = p5, 
       width = 7, height = 8, dpi = 300, units = "in")

