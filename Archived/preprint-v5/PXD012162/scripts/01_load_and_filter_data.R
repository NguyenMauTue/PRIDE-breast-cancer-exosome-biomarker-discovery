############################################################
# 01 Load and filter MaxQuant proteinGroups dataset
############################################################
#Loading library
library(readxl)
library(data.table)

# Load dataset
protein_group <- read_excel(here::here("PXD012162", "data", "pr5c00553_si_001.xlsx"),
                            sheet = "Tab. S1",
                            col_types = "text",
                            na = "")

protein_group$Peptides <- as.numeric(protein_group$Peptides)
protein_group$`Sequence.coverage.[%]` <- as.numeric(protein_group$`Sequence.coverage.[%]`)
protein_group[is.na(protein_group)] <- ""


# Quality filtering
filtered_protein_group = protein_group[
  protein_group$Peptides >= 2 &
  protein_group$`Sequence.coverage.[%]` >= 5 &
  protein_group$Reverse != "+" &
  protein_group$`Potential.contaminant` != "+" &
  protein_group$`Only.identified.by.site` != "+",
]

# Save filtered dataset
fwrite(filtered_protein_group,
       here::here("PXD012162", "data", "filtered_protein_groups.csv"))
  