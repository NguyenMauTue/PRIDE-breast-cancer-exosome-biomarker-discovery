############################################################
# 15 Final tables and visualization
############################################################

library(ggplot2)
library(ggrepel)
library(openxlsx)

filtered_themed =
  read.csv("../results/BiomarkerCandidates_themed.csv")


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
  "results/Localization_tables.xlsx",
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
  scale_color_gradient(low = "grey80", high = "red") +
  scale_fill_gradient(low = "grey80", high = "red") +
  labs(
    x = "Log2 Fold Change",
    y = "Network Degree (log1p)",
    color = "CDS",
    fill = "CDS",
    title = "Driver Candidate Landscape"
  ) +
  theme_classic()

p1

ggsave(
  "../results/DriverScore_landscape_plot.png",
  p1,
  width = 6,
  height = 5
)


#Novelty analysis
filtered_themed$novel =
  filtered_themed$PubMed_hits < quantile(filtered_themed$PubMed_hits, 0.25) &
  filtered_themed$CDS > quantile(filtered_themed$CDS, 0.75)

p2 =
  ggplot(filtered_themed,
         aes(x=PubMed_hits,
             y=CDS)) +
  geom_point(aes(color=novel),
             size=3,
             alpha=0.8) +
  scale_x_log10() +
  geom_hline(yintercept= quantile(filtered_themed$CDS, 0.75),
             linetype="dashed") +
  geom_vline(
    xintercept=
      quantile(filtered_themed$PubMed_hits, 0.25),
    linetype="dashed"
  ) +
  scale_color_manual(values=c("grey70","red")) +
  theme_classic() +
  geom_text_repel(
    data=subset(filtered_themed,novel),
    aes(label=Symbol),
    size=3
  )

p2

ggsave(
  "../results/Literature_novelty_plot.png",
  p2,
  width=6,
  height=5
)
