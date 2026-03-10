############################################################
# 10 Protein-protein interaction network (STRING)
############################################################

# Gene symbols exported from pathway enrichment

gene_symbols =
  read.table("../results/genes_for_STRING.txt",
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
