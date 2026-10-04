library(tidyverse)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(tidyr)
# ── Panel 5.1 — Rank-level bump chart ──────────────────────────────
source(here::here("R", "Helper", "theme_paper.R"))
ranks <- read.csv(here::here('Ablation_Table', "results", "tables", "ablation_ranks_partial.csv"))

k <- 15  

bump_df <- ranks %>%
  filter(rank_AHPCDS <= k) %>%
  mutate(
    direction = ifelse(
      rank_Centralityonly > rank_FConly,
      "Dropped (FC \u2192 Centrality)",
      "Rose (FC \u2192 Centrality)"
    )
  ) %>%
  dplyr::select(Symbol, direction, rank_AHPCDS, rank_FConly, rank_Centralityonly, rank_EqualMCDA) %>%
  pivot_longer(
    cols = starts_with("rank_"),
    names_to = "method",
    values_to = "rank"
  ) %>%
  mutate(
    method = recode(method,
                    rank_AHPCDS         = "AHP-CDS",
                    rank_FConly          = "FC-only",
                    rank_Centralityonly  = "Centrality-only",
                    rank_EqualMCDA       = "Equal-weight"
    ),
    method = factor(method, levels = c("FC-only", "Centrality-only", "Equal-weight", "AHP-CDS"))
  )

pal_direction <- c(
  "Dropped (FC \u2192 Centrality)" = unname(pal_sig["Down"]),
  "Rose (FC \u2192 Centrality)"    = unname(pal_sig["Up"])
)

ht_51 <- ggplot(bump_df, aes(x = method, y = rank, group = Symbol, color = direction)) +
  geom_line(linewidth = 0.5, alpha = 0.75) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_text_repel(
    data = filter(bump_df, method == "AHP-CDS"),
    aes(label = Symbol),
    direction = "y", hjust = 0, nudge_x = 0.15,
    size = 2.5, family = "Arial", color = "black",
    segment.size = 0.2, segment.color = "grey60",
    box.padding = 0.15, min.segment.length = 0
  ) +
  scale_color_manual(values = pal_direction, name = NULL) +
  scale_y_reverse(breaks = seq(1, 38, 5)) +
  coord_cartesian(ylim = c(33.5, 1), clip = "on") +
  labs(x = NULL, y = "Rank") +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7),
    legend.position = "bottom",
    plot.margin = margin(5, 40, 5, 5)
  )
ht_51 <- ht_51 +
  annotate(
    "segment",
    x = 0.85, xend = 1.15, y = 35, yend = 35,
    linetype = "dashed", color = "grey50", linewidth = 0.3
  ) +
  annotate(
    "text", x = 1, y = 34.5, label = "rank > 35 truncated",
    size = 2.4, family = "Arial", color = "grey50", hjust = 0.5
  )

ht_51

# Panel 5.2 — Criterion-level evaluation: contribution-share, Top-15 vs. rest of pool
contrib_raw <- read.csv(here::here("Ablation_Table", "results", "tables", "RQ2A_contribution_shares.csv"))

panel2a_df <- bind_rows(
  contrib_raw %>% filter(rank_AHP <= 10)               %>% mutate(group = "Top-10"),
  contrib_raw %>% filter(rank_AHP <= 15)                %>% mutate(group = "Top-15"),
  contrib_raw %>% filter(rank_AHP > 15)                 %>% mutate(group = "Rest of pool")
) %>%
  group_by(group) %>%
  summarise(
    `Fold-change` = mean(P_FC),
    `FDR`         = mean(P_FDR),
    `Degree`      = mean(P_Deg),
    `Betweenness` = mean(P_Bet),
    .groups = "drop"
  ) %>%
  # Safety net: force each group's shares to sum to exactly 1. Per-candidate
  # P_* already sum to 1 by construction, but summarise(mean()) can pick up a
  # few ULPs of floating-point drift depending on summation order -- without
  # this, a stacked segment landing a hair past the x=1 hard limit below gets
  # silently dropped by scale_x_continuous(limits=...) (NA'd out, not just
  # visually cropped), which is what caused the missing Fold-change segment.
  mutate(row_sum = `Fold-change` + FDR + Degree + Betweenness) %>%
  mutate(across(c(`Fold-change`, FDR, Degree, Betweenness), ~ .x / row_sum)) %>%
  dplyr::select(-row_sum) %>%
  pivot_longer(-group, names_to = "criterion", values_to = "value") %>%
  mutate(
    group = factor(group, levels = c("Rest of pool", "Top-15", "Top-10")),
    criterion = factor(criterion,
                       levels = c("Fold-change", "FDR", "Degree", "Betweenness"))
  )

# --- Palette (match project convention, e.g. pal_criterion alongside
# existing pal_sig / pal_tier / pal_tier_string) ------------------------
pal_criterion <- c(
  "Fold-change"          = "#E69F00",
  "FDR"         = "#56B4E9",
  "Degree"      = "#009E73",
  "Betweenness" = "#CC79A7"
)

# --- Plot ---------------------------------------------------------------
p_panel2a <- ggplot(panel2a_df, aes(x = value, y = group, fill = criterion)) +
  geom_col(position = position_stack(reverse = TRUE),
           width = 0.55, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white", size = 3, fontface = "plain"
  ) +
  scale_fill_manual(values = pal_criterion, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1)) +  # visual crop only -- does NOT drop/NA data
  labs(
    x = "Mean proportional contribution",
    y = NULL,
    title = NULL  # panel label (A) added via patchwork/cowplot at mega-figure assembly
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

p_panel2a

# Panel 2B — Criterion-level evaluation: LOCO delta-rank continuum
loco_raw <- read.csv(here::here("Ablation_Table", "results", "tables", "RQ2B_LOCO_delta_r_wide.csv"))

library(ggplot2)
library(dplyr)

# --- 1. Khai báo Shape & Màu VIỀN tinh gọn cho đúng 4 Archetypes ------------
shape_archetype <- c(
  "FC-rescued"         = 24, # Tam giác hướng lên (▲)
  "Degree-rescued"     = 23, # Hình thoi (◆)
  "Jointly-rescued"    = 22, # Hình vuông (■)
  "Jointly-suppressed" = 25  # Tam giác ngược (▼)
)


# --- 2. Phân loại chuẩn xác 4 góc phần tư -----------------------------------
loco_df <- loco_raw %>%
  mutate(
    archetype = case_when(
      delta_r_without_FC >= 0 & delta_r_without_Deg < 0  ~ "FC-rescued",
      delta_r_without_FC < 0  & delta_r_without_Deg >= 0 ~ "Degree-rescued",
      delta_r_without_FC >= 0 & delta_r_without_Deg >= 0 ~ "Jointly-rescued",
      delta_r_without_FC < 0  & delta_r_without_Deg < 0  ~ "Jointly-suppressed"
    ),
    archetype = factor(archetype, levels = c("FC-rescued", "Degree-rescued", 
                                             "Jointly-rescued", "Jointly-suppressed"))
  )

# --- 3. Plot Panel 2B -------------------------------------------------------
p_panel2b <- ggplot(loco_df, aes(x = delta_r_without_FC, y = delta_r_without_Deg)) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4, linetype = "dashed") +
  
  geom_point(aes(shape = archetype), 
             fill = "white", size = 2.2, stroke = 0.9, alpha = 0.85) +
  
  scale_shape_manual(values = shape_archetype, name = NULL) +
  
  labs(
    x = expression(Delta*italic(r)~"(rank shift when FC removed)"),
    y = expression(Delta*italic(r)~"(rank shift when Degree removed)")
  ) +
  theme_minimal(base_size = 11) +
  theme_paper()

p_panel2b


# 5.3 Panel 1 — Portfolio-level evaluation: G* sweep vs. null band
gstar_raw <- read.csv(here::here("Ablation_Table", "results", "tables", "RQ5_Gstar_sweep_results.csv"))
null_raw  <- read.csv(here::here("Ablation_Table", "results", "tables", "RQ5_null_distribution.csv"))

gstar_plot_df <- gstar_raw %>%
  filter(method %in% c("AHP_CDS", "RWR", "Equal_weight")) %>%
  mutate(method = recode(method,
                         "AHP_CDS" = "AHP-CDS",
                         "RWR" = "RWR",
                         "Equal_weight" = "Equal-weight"),
         method = factor(method, levels = c("AHP-CDS", "RWR", "Equal-weight")))

# --- Palette --------------------------------------------------------
# AHP-CDS is the focal method -> darkest/most saturated + thicker line;
# RWR and Equal-weight (comparators) get lighter/muted tones so AHP-CDS
# reads as the throughline at a glance.
pal_method <- c(
  "AHP-CDS"      = "#440154",
  "RWR"          = "#21908C",
  "Equal-weight" = "#F98E09"
)

# --- Plot -------------------------------------------------------------
p_531 <- ggplot() +
  geom_ribbon(data = null_raw, aes(x = k, ymin = null_p2_5, ymax = null_p97_5),
              fill = "grey70", alpha = 0.35) +
  geom_line(data = null_raw, aes(x = k, y = null_mean),
            color = "grey45", linetype = "dashed", linewidth = 0.4) +
  geom_line(data = gstar_plot_df, aes(x = k, y = G_star, color = method,
                                      linewidth = method)) +
  geom_point(data = gstar_plot_df, aes(x = k, y = G_star, color = method),
             size = 1.8) +
  scale_color_manual(values = pal_method, name = NULL) +
  scale_linewidth_manual(values = c("AHP-CDS" = 1.3, "RWR" = 0.6, "Equal-weight" = 0.6),
                         guide = "none") +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 30)) +
  coord_cartesian(ylim = c(0.6, 1.02)) +
  labs(x = "k (top-k depth)", y = expression(G^"*")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
p_531

# 5.3 Panel 2 — Portfolio-level evaluation: tier-composition stacked bar
gstar_raw <- read.csv(here::here("Ablation_Table", "results", "tables", "RQ5_Gstar_sweep_results.csv"))
tier_df <- gstar_raw %>%
  filter(method %in% c("AHP_CDS", "RWR", "Equal_weight")) %>%
  mutate(
    method = recode(method,
                    "AHP_CDS" = "AHP-CDS",
                    "RWR" = "RWR",
                    "Equal_weight" = "Equal-weight"),
    method = factor(method, levels = c("AHP-CDS", "RWR", "Equal-weight")),
    `Tier 2/3` = tier2 + tier3
  ) %>%
  select(method, k, `Tier 1` = tier1, `Tier 2/3`, `Tier 4` = tier4) %>%
  pivot_longer(cols = c(`Tier 1`, `Tier 2/3`, `Tier 4`),
               names_to = "tier", values_to = "count") %>%
  mutate(
    tier = factor(tier, levels = c("Tier 1", "Tier 2/3", "Tier 4")),
    k = factor(k, levels = c(5, 10, 15, 20, 25, 30))
  )

# --- Palette (pal_tier, 3-bucket scheme shared with 5.3/5.4) -----------
pal_tier <- c(
  "Tier 1"   = "#001959",
  "Tier 2/3" = "#708A41",
  "Tier 4"   = "#FBB680"
)

# --- Plot -------------------------------------------------------------
p_532 <- ggplot(tier_df, aes(x = method, y = count, fill = tier)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7,
           color = "white", linewidth = 0.3) +
  facet_grid(~ k, switch = "x") +
  scale_fill_manual(values = pal_tier, name = NULL) +
  labs(x = "k (top-k depth)", y = "Candidate count") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.ticks.x = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom"
  )
p_532

# ------------------- STRING network panel (igraph + ggraph) --------------
library(here)
library(igraph)
library(ggraph)
library(dplyr)
source(here::here("R", "Helper","string_api_helper.R"))  # build_string_network()

# ── Case-study proteins get emphasis via size/stroke, NOT a new color —
# they keep their pal_tier fill so the tier-color logic stays legible
# across the whole figure.
case_study_ids <- c("COL18A1", "GSTP1", "VCL", "CASP14", "HLA-DRA", "TSG101", "VAMP3")

# ── Build the network via the STRING API helper (return_graph=TRUE) ──────
# candidate_genes: your 143-protein pool's gene symbols
# candidate_tiers: lookup table with columns Symbol, tier (values must
#   match pal_tier's names: "Tier 1"/"Tier 2/3"/"Tier 4")
candidate_genes <- read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_annotated.csv"))  
candidate_tiers <- read.csv(here::here("Ablation_Table", "results", "tables", "candidates_tiered_CTD.csv"))         
net <- build_string_network(candidate_genes, required_score = 700, return_graph = TRUE)
g   <- net$graph

# ── Attach tier metadata to nodes by Symbol ───────────────────────────────
# candidates_tiered_CTD.csv stores tier as numeric (1/2/3/4) — convert to
# the string labels pal_tier expects. Tier 2 and 3 collapse into the same
# "Tier 2/3" bucket used everywhere else in the figure.
node_tiers <- data.frame(Symbol = V(g)$name) %>%
  left_join(candidate_tiers, by = "Symbol") %>%
  mutate(tier_label = case_when(
    tier == 1 ~ "Tier 1",
    tier == 2 ~ "Tier 2",
    tier == 3 ~ "Tier 3",
    tier == 4 ~ "Tier 4",
    TRUE ~ NA_character_   # not in the candidate table (e.g. a STRING
    # first-order neighbor outside the scored pool)
  ))

V(g)$tier <- node_tiers$tier_label
V(g)$is_case_study <- V(g)$name %in% case_study_ids
V(g)$node_size      <- ifelse(V(g)$is_case_study, 4.5, 2.5)
V(g)$stroke_width   <- ifelse(V(g)$is_case_study, 0.8, 0.3)

n_nodes_full <- vcount(g)
n_edges_full <- ecount(g)

# ── Subset to the case-study neighborhood BEFORE layout — laying out all
# 232 nodes and fading the rest wastes most of the panel on scattered,
# unreadable background nodes. This is now the induced subgraph of
# case-study proteins + their direct STRING neighbors; the full network
# (n_nodes_full/n_edges_full above) can go in Supplementary if needed.
case_study_idx <- which(V(g)$is_case_study)
neighbor_idx   <- unique(unlist(ego(g, order = 1, nodes = case_study_idx)))
g_focus        <- induced_subgraph(g, neighbor_idx)

n_nodes <- vcount(g_focus)
n_edges <- ecount(g_focus)

# ── Plot ──────────────────────────────────────────────────────────────────
set.seed(30032026)  # layout is stochastic — fix seed for a reproducible figure
p_string <- ggraph(g_focus, layout = "stress") +
  geom_edge_link(color = "grey70", width = 0.25, alpha = 0.4) +
  geom_node_point(aes(fill = tier, size = node_size, stroke = stroke_width),
                  shape = 21, color = "white") +
  scale_fill_manual(values = pal_tier_string, name = "CTD tier", na.value = "grey85") +
  scale_size_identity() +
  scale_continuous_identity(aesthetics = "stroke", guide = "none") +
  geom_node_text(
    data = function(x) dplyr::filter(x, is_case_study),
    aes(label = name),
    size = 2.4, family = "Arial", color = "grey20",
    repel = TRUE, segment.size = 0.2, max.overlaps = Inf
  ) +
  labs(caption = sprintf("Case-study neighborhood: %d nodes, %d edges (STRING score \u2265 700); full network %d nodes/%d edges",
                         n_nodes, n_edges, n_nodes_full, n_edges_full)) +
  theme_graph_paper() +
  theme(plot.caption = element_text(size = 7, color = "grey30", hjust = 1))

p_string

#________________________________AHP Case-study ___________________________________ #
library(ggplot2)
library(ggrepel)
library(dplyr)
library(colorspace)
source(here::here("R", "Helper", "theme_paper.R"))

CDS_df <- read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_annotated.csv"))


pal_rank_border <- c(
  "Rose"    = "#E66101", 
  "Dropped" = "#2B83BA"  
)


shape_archetype <- c(
  "FC-rescued"         = 24, 
  "Degree-rescued"     = 23, 
  "Jointly-rescued"    = 22, 
  "Jointly-suppressed" = 25  
)

driver_landscape <- function(candidates,
                             logfc_col  = "logFC",
                             degree_col = "degree",
                             cds_col    = "CDS",
                             gene_col   = "Symbol",
                             archetypes = list(
                               `HLA-DRA` = list(label = "FC-rescued",         rank_change = "Rose"),
                               GSTP1     = list(label = "FC-rescued",         rank_change = "Dropped"),
                               CASP14    = list(label = "FC-rescued",         rank_change = "Rose"),
                               TSG101    = list(label = "Degree-rescued",    rank_change = "Rose"),
                               VCL       = list(label = "Degree-rescued",    rank_change = "Dropped"),
                               COL18A1   = list(label = "Jointly-rescued",   rank_change = "Rose"),
                               VAMP3     = list(label = "Jointly-suppressed",rank_change = "Dropped")
                             )) {
  
  candidates$log_degree <- log1p(candidates[[degree_col]])
  candidates$archetype  <- candidates[[gene_col]]
  candidates$archetype[!candidates$archetype %in% names(archetypes)] <- NA
  
  p <- ggplot(candidates, aes(x = .data[[logfc_col]],
                              y = log_degree,
                              fill = .data[[cds_col]])) +
    geom_point(data = candidates[is.na(candidates$archetype), ],
               shape = 21, size = 1.5, alpha = 0.4, color = "transparent")
  
  archetype_data <- candidates %>% 
    filter(!is.na(archetype)) %>% 
    rowwise() %>% 
    mutate(
      label_type  = archetypes[[archetype]]$label,
      rank_status = archetypes[[archetype]]$rank_change,
      
      shape_val   = shape_archetype[[label_type]],
      border_col  = pal_rank_border[[rank_status]],
      
      lbl_text    = paste0(archetype, "\n(", label_type, ")"),
      
      nd_x = case_when(
        archetype == "CASP14" ~ -0.8,
        archetype == "GSTP1"  ~  0.8,
        TRUE ~ 0
      ),
      nd_y = case_when(
        archetype %in% c("CASP14", "GSTP1") ~ 0.4,
        TRUE ~ 0.45
      )
    ) %>% 
    ungroup()
  
  p <- p +
    geom_point(data = archetype_data, 
               aes(shape = I(shape_val), color = I(border_col)),
               size = 3.8, stroke = 1.5, alpha = 0.95) +
    
    geom_text_repel(data = archetype_data,
                    aes(label = lbl_text, color = I(border_col)),
                    size = 3.2, family = "Arial", fontface = "bold",
                    nudge_x = archetype_data$nd_x,
                    nudge_y = archetype_data$nd_y,
                    hjust = 0.5,
                    segment.color = archetype_data$border_col, 
                    segment.size = 0.4,
                    min.segment.length = 0, 
                    box.padding = 0.3,
                    point.padding = 0.3)
  
  p +
    scale_fill_gradient(low = "#E0E0E0", high = "#1A1A1A", name = "CDS",
                        breaks = c(0.25, 0.50, 0.75),
                        guide = guide_colorbar(barwidth = 8, barheight = 0.6)) +
    labs(x = "log2 Fold-Change", y = "log1p(Degree)") +
    theme_paper()
}

scatter_plot <- driver_landscape(CDS_df)
scatter_plot

# 5.4 Panel 3 — Case studies: contribution fingerprint
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)

case_studies <- c("COL18A1", "GSTP1", "CASP14", "HLA-DRA", "VCL", "TSG101", "VAMP3")

archetype_info <- tibble::tribble(
  ~Symbol,   ~archetype,           ~rank_change, ~shape_icon,
  "COL18A1", "Jointly-rescued",    "Rose",       "■",
  "HLA-DRA", "FC-rescued",         "Rose",       "▲",
  "CASP14",  "FC-rescued",         "Rose",       "▲",
  "GSTP1",   "FC-rescued",         "Dropped",    "▲",
  "TSG101",  "Degree-rescued",    "Rose",       "◆",
  "VCL",     "Degree-rescued",    "Dropped",    "◆",
  "VAMP3",   "Jointly-suppressed", "Dropped",    "▼"
)


fingerprint_df <- contrib_raw %>%
  filter(Symbol %in% case_studies) %>%
  left_join(archetype_info, by = "Symbol") %>%
  mutate(
    y_label = paste0(shape_icon, "  ", Symbol),
    rank_change = factor(rank_change, levels = c("Rose", "Dropped"))
  ) %>%
  select(Symbol, y_label, rank_AHP, rank_change,
         FC = P_FC, FDR = P_FDR, Degree = P_Deg, Betweenness = P_Bet) %>%
  pivot_longer(cols = c(FC, FDR, Degree, Betweenness),
               names_to = "criterion", values_to = "value") %>%
  mutate(criterion = factor(criterion, levels = c("FC", "FDR", "Degree", "Betweenness")))


pal_criterion <- c(
  "FC"          = "#E69F00",
  "FDR"         = "#56B4E9",
  "Degree"      = "#009E73",
  "Betweenness" = "#CC79A7"
)


p_heatmap_colored <- ggplot(fingerprint_df, aes(x = criterion, y = y_label)) +
  geom_tile(aes(fill = criterion, alpha = value), color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", value)),
            color = "black", size = 3.2, fontface = "bold") +
  
  facet_wrap2(
    ~ rank_change, 
    scales = "free_y",
    strip = strip_themed(
      background_x = list(
        element_rect(fill = "#D55E00", color = NA),
        element_rect(fill = "#0072B2", color = NA)  
      )
    )
  ) +
  
  
  scale_fill_manual(values = pal_criterion, name = "Criterion") +
  scale_alpha_continuous(range = c(0.25, 1.0), guide = "none") +
  
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme_paper() 

p_heatmap_colored


