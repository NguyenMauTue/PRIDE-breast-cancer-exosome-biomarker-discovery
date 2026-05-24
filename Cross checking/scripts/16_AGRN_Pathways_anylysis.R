########################################################
# Script 16: AGRN Pathway Context & Co-protein Analysis
# Sentinel proteins — satellite proteins around the hub
# Goal: Which pathways contain AGRN? Are co-occurring
#        proteins dysregulated in the same direction?
########################################################

library(ReactomePA)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(tibble)
library(clusterProfiler)
library(reactome.db)

# ── 0. Load dependencies ──────────────────────────────
de_results   <- read.csv("../results/limma_network_table.csv")        # from script 06
gsea_results <- readRDS("../results/reactome_gsea_object.rds")      # from script 08

# ── 1. AGRN Entrez ID ─────────────────────────────────
agrn_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "O00468",
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "UNIPROT"
)

# ── 2. All Reactome pathways containing AGRN ──────────
# Thay getPathway bằng enrichPathway rồi filter
# Query tất cả pathways chứa AGRN bằng cách 
# check reactome pathway gene sets trực tiếp
# Lấy tất cả pathways → genes mapping
pathway2gene <- as.data.frame(reactomePATHID2EXTID)
colnames(pathway2gene) <- c("pathway_id", "entrez")

# Filter pathways chứa AGRN
agrn_pathways <- pathway2gene %>%
  filter(entrez == agrn_entrez$ENTREZID)

# Lấy pathway names
pathway_names <- as.data.frame(reactomePATHID2NAME)
colnames(pathway_names) <- c("pathway_id", "pathway_name")

agrn_pathways <- agrn_pathways %>%
  left_join(pathway_names, by = "pathway_id") %>%
  filter(grepl("Homo sapiens", pathway_name))  # only human pathways

# ── 3. For each pathway: pull all member genes ────────
# then intersect with our DE dataset

# All proteins in DE results → Entrez
dataset_uniprot <- de_results$UNIPROT

dataset_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = dataset_uniprot,
  columns = c("ENTREZID", "SYMBOL", "UNIPROT"),
  keytype = "UNIPROT"
) %>% 
  filter(!is.na(ENTREZID)) %>%
  distinct(UNIPROT, .keep_all = TRUE)

# Pathway gene lists from Reactome
# Thay toàn bộ lapply block bằng cái này — không cần getGene
pathway_names <- as.data.frame(reactomePATHID2NAME)


satellite_list <- lapply(agrn_pathways$pathway_id, function(pid) {
  
  # Lấy tất cả genes trong pathway từ pathway2gene trực tiếp
  pathway_genes <- pathway2gene %>%
    filter(pathway_id == pid)
  
  if (nrow(pathway_genes) == 0) return(NULL)
  
  # Intersect với dataset
  co_proteins <- dataset_entrez %>%
    filter(ENTREZID %in% pathway_genes$entrez)
  
  if (nrow(co_proteins) == 0) return(NULL)
  
  # Pull DE stats
  co_de <- de_results %>%
    filter(UNIPROT %in% co_proteins$UNIPROT) %>%
    left_join(co_proteins %>% dplyr::select(UNIPROT, SYMBOL), by = "UNIPROT") %>%
    dplyr::select(UNIPROT, SYMBOL, logFC, adj.P.Val) %>%
    distinct(UNIPROT, .keep_all = TRUE) %>%  # thêm dòng này
    mutate(
      pathway_id   = pid,
      pathway_name = agrn_pathways$pathway_name[agrn_pathways$pathway_id == pid][1],
      is_AGRN      = UNIPROT == "O00468",
      direction    = case_when(logFC > 0 ~ "up", logFC < 0 ~ "down", TRUE ~ "ns"),
      significant  = adj.P.Val < 0.05
    )
  
  return(co_de)
})

satellite_df <- bind_rows(satellite_list) %>%
  arrange(.data$pathway_name, desc(abs(logFC)))

satellite_df <- bind_rows(satellite_list) %>%
  distinct(UNIPROT, pathway_id, .keep_all = TRUE) %>%  # dup trong cùng pathway
  arrange(.data$pathway_name, desc(abs(logFC)))
  
# ── 4. Directionality summary per pathway ────────────
# Key question: are satellites going the SAME direction as AGRN (up)?

dir_summary <- satellite_df %>%
  filter(!is_AGRN) %>%
  group_by(pathway_name) %>%
  summarise(
    n_proteins   = n(),
    n_up         = sum(direction == "up"),
    n_down       = sum(direction == "down"),
    n_sig        = sum(significant),
    pct_same_dir = round(n_up / n_proteins * 100, 1),  # AGRN is upregulated
    .groups = "drop"
  ) %>%
  arrange(desc(pct_same_dir))

# ── 5. Save outputs ───────────────────────────────────
saveRDS(satellite_df,  "../data/agrn_satellite_proteins.rds")
write.csv(satellite_df,  "../data/agrn_satellite_proteins.csv",  row.names = FALSE)
write.csv(dir_summary,   "../data/agrn_pathway_dir_summary.csv", row.names = FALSE)

# ── 6. Plot: logFC of co-proteins per pathway ────────
# Dot plot — one row per protein, colored by direction, AGRN highlighted

library(ggplot2)
library(tidyr)
library(dplyr)

# Prep data — long format cho stacked bar
dir_plot <- dir_summary %>%
  arrange(pct_same_dir) %>%
  mutate(pathway_short = gsub("Homo sapiens: ", "", pathway_name),
         pathway_short = factor(pathway_short, levels = pathway_short)) %>%
  dplyr::select(pathway_short, n_up, n_down, pct_same_dir) %>%
  pivot_longer(cols = c(n_up, n_down), 
               names_to = "direction", 
               values_to = "n") %>%
  mutate(direction = ifelse(direction == "n_up", "up", "down"))

# Plot
max_n <- max(dir_summary$n_proteins)  # = 45
p16 <- ggplot(dir_plot, aes(x = n, y = pathway_short, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +
  geom_point(data = dir_plot %>% distinct(pathway_short, pct_same_dir),
             aes(x = pct_same_dir / 100 * max_n,  
                 y = pathway_short),
             inherit.aes = FALSE,
             color = "#4A7C59", size = 2.5) +
  geom_line(data = dir_plot %>% distinct(pathway_short, pct_same_dir),
            aes(x = pct_same_dir / 100 * max_n,   
                y = as.numeric(pathway_short)),
            inherit.aes = FALSE,
            color = "#C9A84C", linewidth = 0.8) +
  scale_fill_manual(values = c(up = "#E8A0B4", down = "#4A6FA5")) +
  scale_x_continuous(
    name = "n proteins",
    breaks = seq(0, max_n, by = 5),
    limits = c(0, max_n),
    sec.axis = sec_axis(~ . / max_n * 100,
                        name = "% same direction as AGRN",
                        breaks = seq(0, 100, by = 20),
                        labels = function(x) paste0(x, "%"))
  ) +
  labs(y = NULL, fill = "direction",
       title = "AGRN satellite proteins: directionality by pathway",
       subtitle = "Bar = protein count | Line = % concordant with AGRN (upregulated)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12),
        legend.position = "bottom")

ggsave("../results/16_AGRN_satellite_directionality.png",
       p16, width = 9, height = 7, dpi = 300)

cat("\nScript 16 complete. Outputs:\n")
cat("  ../data/agrn_satellite_proteins.csv\n")
cat("  ../data/agrn_pathway_dir_summary.csv\n")
cat("  ../plots/16_AGRN_satellite_logFC.png\n")
