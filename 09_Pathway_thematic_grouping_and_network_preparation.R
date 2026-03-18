############################################################
# 09 Pathway thematic grouping and network preparation
############################################################

library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

gsea_reactome_res =
  read.csv("../results/reactome_gsea_full.csv")


############################################################
# Group pathways into biological themes
############################################################

gsea_reactome_res_group = gsea_reactome_res

gsea_reactome_res_group$Biological_Theme = NA


# Vesicle trafficking

gsea_reactome_res_group$Biological_Theme[
  grepl("vesicle|endosomal|ESCRT|exosome",
        gsea_reactome_res_group$Description,
        ignore.case = TRUE)
] <- "Vesicle trafficking"


# Motility signaling

gsea_reactome_res_group$Biological_Theme[
  grepl("RHO|motility|migration|cytoskeleton",
        gsea_reactome_res_group$Description,
        ignore.case = TRUE)
] <- "Motility signaling"


# Cell adhesion

gsea_reactome_res_group$Biological_Theme[
  grepl("integrin|adhesion|junction|ECM",
        gsea_reactome_res_group$Description,
        ignore.case = TRUE)
] <- "Cell adhesion"


# Immune evasion

gsea_reactome_res_group$Biological_Theme[
  grepl("interferon|cytokine|antigen|MHC|PD|CTLA",
        gsea_reactome_res_group$Description,
        ignore.case = TRUE)
] <- "Immune evasion"


############################################################
# Focused pathways
############################################################

gsea_reactome_res_focused =
  gsea_reactome_res_group[
    !is.na(gsea_reactome_res_group$Biological_Theme),
  ]


############################################################
# Prepare plotting dataframe
############################################################

df_plot <-
  gsea_reactome_res_focused %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES)) %>%
  mutate(
    Description =
      gsub("Signaling by ", "", Description)
  ) %>%
  mutate(
    Description =
      gsub("Regulation of ", "", Description)
  ) %>%
  mutate(
    Description =
      str_wrap(Description, 35)
  ) %>%
  mutate(
    Description =
      forcats::fct_reorder(Description, NES)
  )


############################################################
# Pathway enrichment dotplot
############################################################

p1 =
  ggplot(
    df_plot,
    aes(
      x = NES,
      y = Description,
      size = setSize,
      color = p.adjust
    )
  ) +
  geom_point() +
  scale_color_gradient(
    low = "red",
    high = "blue",
    name = "Adjusted P-value"
  ) +
  scale_size_continuous(
    range = c(3,10),
    name = "Gene count"
  ) +
  labs(
    title = "Top enriched Reactome pathways",
    x = "Normalized Enrichment Score",
    y = "Pathway description"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title =
      element_text(hjust = 0.5,
                   face = "bold")
  )


ggsave(
  "../results/reactome_pathway_dotplot.png",
  p1,
  width = 8,
  height = 6
)


############################################################
# Theme-level summary plot
############################################################

p2 =
  ggplot(
    df_plot,
    aes(
      x = Biological_Theme,
      y = NES,
      size = setSize,
      color = p.adjust
    )
  ) +
  geom_point() +
  theme_minimal()


ggsave(
  "../results/reactome_theme_plot.png",
  p2,
  width = 6,
  height = 4
)


############################################################
# Extract genes from core enrichment
############################################################

gene_lists =
  strsplit(
    gsea_reactome_res_focused$core_enrichment,
    "/"
  )

genes =
  unique(
    unlist(gene_lists)
  )


############################################################
# Convert ENTREZ -> SYMBOL
############################################################

symbols =
  bitr(
    genes,
    fromType = "ENTREZID",
    toType   = "SYMBOL",
    OrgDb    = org.Hs.eg.db
  )

gene_symbols =
  symbols$SYMBOL


############################################################
# Export gene list for STRING network
############################################################

gene_lists =
  strsplit(
    gsea_reactome_res$core_enrichment, 
    "/"
  )
