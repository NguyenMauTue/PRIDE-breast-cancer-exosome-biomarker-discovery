############################################################
# Master run_all.R
# Runs the full AHP-CDS pipeline: PXD056161 -> Ablation_Table
# Prerequisites: renv::restore() to install dependencies.
# Open AHP-CDS.Rproj first so here::here() resolves correctly
# (no setwd() needed).
############################################################
# ── 0. Preamp ──────────────────────
source(here::here("R", "Preamp", "fetch_annotation_table.R"))

# ── 1. PXD056161 ──────────────────────
message("\n", strrep("=", 60))
message("PXD056161 PIPELINE")
message(strrep("=", 60))
source(here::here("PXD056161", "scripts", "run_all.R"))

# ── 2. Ablation_Table (depends on PXD056161 outputs) ────
message("\n", strrep("=", 60))
message("ABLATION TABLE PIPELINE")
message(strrep("=", 60))
source(here::here("Ablation_Table", "scripts", "run_all.R"))

# ── Done ──────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("ALL PIPELINES COMPLETE.")
message("Results in: PXD056161/results/, PXD012162/results/, Ablation_Table/results/")
message(strrep("=", 60))
