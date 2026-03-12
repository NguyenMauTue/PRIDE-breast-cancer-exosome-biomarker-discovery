############################################################
# 10 Protein-protein interaction network (STRING)
############################################################

# Gene symbols exported from pathway enrichment

gene_symbols =
  read.table("results/genes_for_STRING.txt",
             stringsAsFactors = FALSE)[,1]

# Network construction was performed using the STRING database
# https://string-db.org

# Procedure:
# 1. Upload gene list (Homo sapiens)
# 2. Interaction score threshold: 0.7 (high confidence)
# 3. Remove disconnected nodes
# 4. Export interaction network (TSV)

# The downloaded network file is saved as:

# ../data/string_network.tsv


############################################################
# Extract STRING PPI network
############################################################

ppi = read.delim("../data/string_interactions.tsv")

ppi_high = ppi[ppi$combined_score > 0.7, ]

library(igraph)

g = graph_from_data_frame(
  ppi_high[,c("X.node1","node2")],
  directed = FALSE
)


############################################################
# Network size
############################################################

vcount(g)   # number of proteins
ecount(g)   # number of interactions


############################################################
# Connected components
############################################################

comp <- components(g)

comp$csize


############################################################
# Network centrality
############################################################

deg = degree(g)

bet = betweenness(g)


############################################################
# Network summary table
############################################################

network_summary = data.frame(
  Symbol = names(deg),
  degree = deg,
  betweenness = bet
)
############################################################
# Save results
############################################################

write.csv(network_summary,
          "../results/network_summary.csv",
          row.names = FALSE)
