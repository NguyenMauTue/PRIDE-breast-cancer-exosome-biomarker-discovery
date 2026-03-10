############################################################
# 01 Load and filter MaxQuant proteinGroups dataset
############################################################
#Loading library
library(data.table)

# Load dataset
protein_group = fread(
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD056161/proteinGroups.txt"
)

# Quality filtering
filtered_protein_group = protein_group[
  Peptides >= 2 &
  `Sequence coverage [%]` >= 5 &
  Reverse != "+" &
  `Potential contaminant` != "+" &
  `Only identified by site` != "+"
]

# Save filtered dataset
fwrite(filtered_protein_group,
       "../data/filtered_protein_groups.csv")
