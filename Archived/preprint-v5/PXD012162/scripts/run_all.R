############################################################
# run_all.R — PXD012162 pipeline runner
# Sources scripts 01-15 in order. Requires AHP-CDS.Rproj to be
# open (or working directory anywhere inside the project) so
# here::here() resolves correctly.
############################################################

experiment <- "PXD012162"

scripts <- c(
  "01_load_and_filter_data.R",
  "02_metadata_and_matrix.R",
  "03_quality_control.R",
  "04_missingness_analysis.R",
  "05_imputation.R",
  "06_Differential_expression_analysis.R",
  "07_Robustness_analysis.R",
  "08_Pathway_enrichment_analysis.R",
  "09_Extract_genes_for_STRING.R",
  "10_Protein-protein_interaction_network.R",
  "11_Network-expression_integration.R",
  "12_CDS_via_AHP.R",
  "14_Biological_theme_classification.R",
  "15_Final_tables_and_visualization.R"
)

run_script <- function(experiment, script_name) {
  message(strrep("-", 60))
  message(sprintf("[%s] Running %s ...", experiment, script_name))
  message(strrep("-", 60))

  path <- here::here(experiment, "scripts", script_name)

  result <- tryCatch(
    {
      source(path, local = new.env())
      TRUE
    },
    error = function(e) {
      message("\n*** FAILED: ", script_name, " ***")
      message("Error: ", conditionMessage(e))
      FALSE
    }
  )

  if (!result) {
    stop(sprintf(
      "Pipeline stopped at %s (%s). Fix the error above, then re-run from this script onward.",
      script_name, experiment
    ), call. = FALSE)
  }

  message(sprintf("[%s] Done: %s\n", experiment, script_name))
}

# ── Steps 01-15 ──────────────────────────────────────────
for (s in scripts) run_script(experiment, s)


message(strrep("=", 60))
message(sprintf("[%s] PIPELINE COMPLETE", experiment))
message(strrep("=", 60))
