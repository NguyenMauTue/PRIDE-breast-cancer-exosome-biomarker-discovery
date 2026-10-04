############################################################
# run_all.R — Master pipeline script
# Prerequisites: renv::restore() to install dependencies
############################################################
options(timeout = 3600)
scripts <- c(
  "01_load_and_filter_data.R",
  "02_metadata_and_matrix.R",
  "03_quality_control.R",
  "04_missingness_analysis.R",
  "05_imputation.R",
  "06_differential_expression_analysis.R",
  "07_robustness_analysis.R",
  "08_pathway_enrichment_analysis.R",
  "09_Extract_genes_for_STRING.R",
  "10_Protein-protein_interaction_network.R",
  "11_Network-expression_integration.R",
  "12_CDS_via_AHP.R",
  "13_Functional_annotation.R",
  "14_Biological_theme_classification.R",
  "15_Final_tables_and_visualization.R",
  "16_AGRN_Pathways_anylysis.R"
)

for (script in scripts) {
  message("\n>>> Running: ", script)
  tryCatch(
    source(script, chdir = TRUE),
    error = function(e) {
      message("ERROR in ", script, ": ", e$message)
      stop(e)
    }
  )
  message("<<< Done: ", script)
}

# Export session info after all packages loaded
sink("../results/sessionInfo.txt")
print(sessionInfo())
sink()

message("\nPipeline complete. Results in results/")
