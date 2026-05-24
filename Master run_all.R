############################################################
# Master run_all.R
# Runs complete pipeline: main + cross-validation
# Prerequisites: renv::restore() to install dependencies
############################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(timeout = 3600)

# ── 1. Main pipeline ────────────────────────────────────
message("\n", strrep("=", 50))
message("MAIN PIPELINE")
message(strrep("=", 50))
source("scripts/run all.R", chdir = TRUE)

# ── 2. Cross-validation pipeline ────────────────────────
message("\n", strrep("=", 50))
message("CROSS-VALIDATION PIPELINE")
message(strrep("=", 50))
source("Cross checking/scripts/run all.R", chdir = TRUE)

# ── 3. Visualization ────────────────────────
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
message("\n", strrep("=", 50))
message("Final_visualization")
message(strrep("=", 50))
source("Visualization.R")
# ── 3. Done ─────────────────────────────────────────────
message("\n", strrep("=", 50))
message("ALL PIPELINES COMPLETE.")
message("Results in: results/ and Cross checking/results/ and main folder")
message(strrep("=", 50))
