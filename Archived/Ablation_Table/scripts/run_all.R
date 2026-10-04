############################################################
# run_all.R — Ablation_Table pipeline runner
# Sources scripts 04-16 in order (skips 06a, which is an
# unused exploration script, not part of the pipeline).
#
# PREREQUISITE: both PXD056161/scripts/run_all.R and
# PXD012162/scripts/run_all.R must already have been run —
# several scripts here read their results/data directly
# (limma_network_table.csv, BiomarkerCandidates_themed.csv,
# network_summary.csv, differential_expression_imputed.csv,
# imputed_matrix.rds, string_interactions.tsv from PXD056161;
# Module_tables.xlsx from PXD012162).
############################################################

experiment <- "Ablation_Table"

scripts <- c(
  "04_ablation_conditions.R",
  "07_ren2019_gr.R",
  "RQ1_ranking_difference.R",
  "RQ2_mechanism_analysis.R",
  "RQ5_tier_audit_CTD.R"
)

# ── Preflight checks ─────────────────────────────────────
#required_upstream <- c(

#)

#missing <- required_upstream[!file.exists(required_upstream)]

if (length(missing) > 0) {
  message("Missing required input file(s):")
  for (m in missing) message("  - ", m)
  stop(
    "Run PXD056161/scripts/run_all.R and PXD012162/scripts/run_all.R first — Ablation_Table depends on both experiments' outputs.",
    call. = FALSE
  )
}

# ── Run scripts ───────────────────────────────────────────
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

for (s in scripts) run_script(experiment, s)

message(strrep("=", 60))
message(sprintf("[%s] PIPELINE COMPLETE", experiment))
message(strrep("=", 60))
