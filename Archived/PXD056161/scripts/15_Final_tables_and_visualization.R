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
library(circlize)
library(xtable)

filtered_themed <-
  read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_themed.csv"))
imputed_matrix <-
  readRDS(here::here("PXD056161", "data", "imputed_matrix.rds"))
imputed_result <-
  read.csv(here::here("PXD056161", "results", "tables", "differential_expression_imputed.csv"))

# Summary table — đổi Localization → Biological_Module + robustness_label
summarized_df <-
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
      "Biological_Module"
    )
  ]

# Abbreviate robustness_label
summarized_df <- summarized_df |>
  mutate(robustness_label = case_when(
    robustness_label == "robust_candidate"          ~ "R",
    robustness_label == "weight_sensitive_candidate" ~ "WS",
    TRUE ~ robustness_label
  ))

# Split by Biological_Module cho Excel
listModules <- unique(summarized_df$Biological_Module)

# Excel export
wb <- createWorkbook()

addWorksheet(wb, "All")
writeData(wb, "All", summarized_df)

for (mod in listModules) {
  clean_name <- gsub("[^A-Za-z0-9]", "_", mod)
  clean_name <- substr(clean_name, 1, 31)  # Excel sheet name max 31 chars
  addWorksheet(wb, clean_name)
  writeData(wb, clean_name, subset(summarized_df, Biological_Module == mod))
}

saveWorkbook(
  wb,
  here::here("PXD056161", "results", "tables", "Module_tables.xlsx"),
  overwrite = TRUE
)

# LaTeX export
latex_df <- filtered_themed |>
  arrange(desc(CDS)) |>
  mutate(
    robustness_short = case_when(
      robustness_label == "robust_candidate"           ~ "R",
      robustness_label == "weight_sensitive_candidate" ~ "WS",
      TRUE ~ robustness_label
    ),
    logFC   = round(logFC, 3),
    adj.P.Val = formatC(adj.P.Val, format = "e", digits = 1),
    CDS     = round(CDS, 4)
  ) |>
  dplyr::select(Symbol, logFC, adj.P.Val, CDS, Biological_Module, robustness_short)

colnames(latex_df) <- c(
  "Gene",
  "$\\log_2$FC",
  "FDR",
  "CDS",
  "Module",
  "Rob."
)

print(
  xtable(
    latex_df,
    tabular.environment = "longtable",
    floating = FALSE,
    caption = "Full list of driver candidates ranked by AHP-weighted Composite Driver Score (CDS). 
    Biological modules assigned via GO term hierarchy (GO.db): ECM = ECM \\& Cell Adhesion 
    (GO:0031012, GO:0030198, GO:0007160); MS = Motility \\& Signaling (GO:0016477, GO:0000165, 
    GO:0007265, GO:0035023); VT = Vesicle Trafficking (GO:0016192, GO:0036258). 
    Robustness: R = robust candidate; WS = weight-sensitive candidate.",
    label = "tab:candidates",
    digits = c(0, 0, 3, 0, 4, 0, 0)
  ),
  include.rownames  = FALSE,
  sanitize.colnames.function = identity,
  sanitize.text.function = identity,
  file = here::here("PXD056161", "results", "tables", "Supplementary_Table1.tex")
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
  here::here("PXD056161", "results", "figures", "DriverScore_landscape_plot.png"),
  p1,
  width = 6,
  height = 5,
  dpi = 300
)

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
  here::here("PXD056161", "results", "figures", "Volcano_Plot.png"),
  volcano,
  width=6,
  height=5,
  dpi = 600
)

