############################################################
# 15 Final tables and visualization
############################################################

library(ggplot2)
library(ggrepel)
library(openxlsx)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(patchwork)
library(circlize)


filtered_themed =
  read.csv("../results/BiomarkerCandidates_themed.csv")
imputed_matrix =
  readRDS("../data/imputed_matrix.rds")
imputed_result =
  read.csv("../results/differential_expression_imputed.csv")


#Summary table

summarized_df =
  filtered_themed[
    ,
    c(
      "UNIPROT",
      "Symbol",
      "logFC",
      "adj.P.Val",
      "degree",
      "betweenness",
      "CDS",
      "robustness_label",
      "Localization",
      "PubMed_hits"
    )
  ]


#Split by localization
listLocal = unique(summarized_df$Localization)

get_Subset_Local = function(df,loc){
  
  clean_name = gsub("[^A-Za-z0-9]","_",loc)
  
  assign(
    paste0(clean_name,"_df"),
    subset(df,df$Localization==loc),
    envir=.GlobalEnv
  )
}

for(loc in listLocal){
  get_Subset_Local(summarized_df,loc)
}

#Excel export
wb = createWorkbook()

addWorksheet(wb,"All")
writeData(wb,"All",summarized_df)

addWorksheet(wb,"Membrane")
writeData(wb,"Membrane",membrane_df)

addWorksheet(wb,"Extracellular")
writeData(wb,"Extracellular",extracellular_df)

addWorksheet(wb,"Vesicle_exosome")
writeData(wb,"Vesicle_exosome",vesicle_exosome_df)

saveWorkbook(
  wb,
  "../results/Localization_tables.xlsx",
  overwrite=TRUE
)

#Plot driver landscape
top20_df <- summarized_df |>
  arrange(desc(CDS)) |>
  head(20)

p1 <- ggplot(summarized_df,
             aes(x = logFC, y = log1p(degree), color = CDS)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_point(data = top20_df,
             size = 3.5,
             shape = 21,
             color = "black",
             aes(fill = CDS)) +
  geom_text_repel(
    data = top20_df,
    aes(label = Symbol),
    size = 3,
    color = "black"
  ) +
  scale_color_gradient(low = "grey80", high = "red",
                       aesthetics = c("color", "fill"),
                       name = "CDS") +
  labs(
    x = "Log2 Fold Change",
    y = "Network Degree (log1p)",
    color = "CDS",
    title = "Driver Candidate Landscape"
  ) +
  theme_classic() +
  theme(text = element_text())

p1

ggsave(
  "../results/DriverScore_landscape_plot.png",
  p1,
  width = 6,
  height = 5,
  dpi = 300
)


#Novelty analysis
filtered_themed <- filtered_themed |>
  mutate(novelty_label = case_when(
    PubMed_hits < 10   ~ "Understudied",
    PubMed_hits < 100  ~ "Emerging",
    PubMed_hits < 1000 ~ "Established",
    TRUE               ~ "Well-studied"
  ))
filtered_themed$novelty_label = factor(filtered_themed$novelty_label, levels = c("Understudied", "Emerging", "Established", "Well-studied"))
filtered_themed$driver =
  filtered_themed$CDS > quantile(filtered_themed$CDS, 0.75, na.rm = TRUE)

p2 =
  ggplot(filtered_themed,
         aes(x=PubMed_hits,
             y=CDS)) +
  geom_point(aes(color=novelty_label),
             size=3,
             alpha=0.8) +
  scale_x_log10() +
  geom_hline(yintercept= quantile(filtered_themed$CDS, 0.75),
             linetype="dashed") +
  scale_color_manual(values=c("#F4C430","#C49A28", "#8B6914", "#4A3600")) +
  theme_classic() +
  geom_text_repel(
    data = subset(filtered_themed, driver),
    aes(label = Symbol),
    size = 3,
    force = 3,
    max.overlaps = 15
  ) +
  theme(text = element_text()) +
  labs(
    title = "Novelty and Priority Assessment of Candidate Proteins" 
  )
  

p2

ggsave(
  "../results/Literature_novelty_plot.png",
  p2,
  width=6,
  height=5,
  dpi = 300
)

p12 = p1 + p2 + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A')
p12

ggsave(
  "../results/Literature_novelty_plot.png",
  p12,
  width=12,
  height=6,
  dpi = 600
)

#pheatmap
candidate_proteins <- filtered_themed$UNIPROT
match_rows <- sapply(rownames(imputed_matrix), function(x) {
  ids <- strsplit(x, ";")[[1]]
  any(ids %in% candidate_proteins)
})

logFC <- imputed_matrix[rownames(imputed_matrix) %in% candidate_proteins, ]
id_to_symbol <- setNames(filtered_themed$Symbol, filtered_themed$UNIPROT)
rownames(logFC) <- id_to_symbol[rownames(logFC)]
scaledlogFC = t(scale(t(logFC)))

png("../results/Heatmap_candidates.png", 
    width = 8, height = 10, units = "in", res = 300)
hp_sorted <- Heatmap(scaledlogFC,
                     name = "Z-score",
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     column_order = 1:6,
                     cluster_rows = TRUE,
                     show_row_dend = TRUE,
                     cluster_columns = FALSE,
                     row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                     column_names_gp = gpar(fontsize = 10),
                     rect_gp = gpar(col = "white", lwd = 1),
                     top_annotation = HeatmapAnnotation(
                       Group = c(rep("Tumor", 3), rep("Normal", 3)),
                       col = list(Group = c(Normal = "#1A9E8F", Tumor = "#7B2D8B"))
                     ))
draw(hp_sorted)
dev.off()

#Volcano plot
# Join symbol
imputed_result <- imputed_result |>
  left_join(filtered_themed |> dplyr::select(UNIPROT, Symbol), by = "UNIPROT")
# Category
imputed_result <- imputed_result |>
  mutate(category = case_when(
    adj.P.Val >= 0.05 | abs(logFC) < 1 ~ "Non-significant",
    UNIPROT %in% filtered_themed$UNIPROT ~ "Candidate",
    TRUE ~ "DE only"
  ))

# Plot
volcano <- ggplot(imputed_result, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = category), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(
    "Non-significant" = "grey70",
    "DE only"         = "#6BAED6",
    "Candidate"       = "#D6604D"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(
    data = subset(imputed_result, category == "Candidate"),
    aes(label = Symbol),
    size = 3,
    max.overlaps = 20
  ) +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       color = "Category")

volcano

ggsave(
  "../results/Volcano Plot.png",
  volcano,
  width=6,
  height=5,
  dpi = 600
)

