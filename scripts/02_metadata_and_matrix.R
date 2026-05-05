############################################################
# 02 Construct LFQ matrix and sample metadata
############################################################

library(data.table)

filtered_protein_group =
  fread("../data/filtered_protein_groups.csv")

# Identify LFQ columns
lfq_cols = grep("^LFQ intensity",
                colnames(filtered_protein_group),
                value = TRUE)

clean_names = gsub("LFQ intensity ", "", lfq_cols)

# Parse metadata
info = strsplit(clean_names, "_")

metadata = data.frame(
  Sample = clean_names,
  Group = sapply(info, function(x)
    ifelse(x[2] == "MCF10", "Normal", "Tumor")),
  BioRep = sapply(info, function(x) x[3]),
  TechRep = sapply(info, function(x) x[4]),
  BioRep_full = sapply(info, function(x)
    paste0(ifelse(x[2] == "MCF10", "Normal", "Tumor"), x[3]))
)

# Build LFQ matrix
lfq_matrix = as.matrix(
  filtered_protein_group[, lapply(.SD, as.numeric),
                         .SDcols = lfq_cols]
)

rownames(lfq_matrix) =
  make.unique(filtered_protein_group$`Protein IDs`)

colnames(lfq_matrix) = clean_names

# Convert zeros to NA
lfq_matrix[lfq_matrix == 0] = NA

# Log transform
log_lfq_matrix = log2(lfq_matrix)

saveRDS(log_lfq_matrix,
        "../data/log_lfq_matrix.rds")

write.csv(metadata,
          "../data/sample_metadata.csv",
          row.names = FALSE)
