# ============================================================
# AHP-CDS Case-Study Selection Audit
# ============================================================
# Purpose:
#   Apply a predefined, deterministic rule to the full candidate
#   pool and identify candidates whose baseline direction and
#   FC-LOCO direction are discordant.
#
# IMPORTANT:
#   This script DOES NOT select Case Study 7.
#   It produces an auditable classification and shortlist for
#   researcher inspection.
#
# Required input columns:
#   candidate       : unique candidate/protein identifier
#   rank_fc         : rank under FC-only baseline
#   rank_centrality : rank under Centrality-only baseline
#   delta_r_no_fc   : precomputed Delta r without FC
#   c_fc            : FC contribution to composite score
#   c_deg           : Degree contribution to composite score
#
# Baseline rule:
#   baseline_delta = rank_centrality - rank_fc
#   > 0  => Dropped
#   < 0  => Rose
#
# LOCO rule:
#   delta_r_no_fc > 0 => FC-rescued
#   delta_r_no_fc < 0 => network-rescued
#
# No near-neutral category is used.
# Therefore delta_r_no_fc == 0 is classified as "zero_change"
# and is NOT treated as rescued or discordant.
#
# Contribution rule:
#   c_deg == 0 AND c_fc > 0
#       => FC-dominant
#   c_fc == 0 AND c_deg > 0
#       => network-dominant
#   c_fc > 0 AND c_deg > 0 AND
#       abs(c_fc - c_deg) <= Q1(abs(c_fc - c_deg))
#       => convergent
#   otherwise
#       => mixed_nonconvergent
#
# Q1 is calculated from the ENTIRE candidate pool using
# abs(c_fc - c_deg), before any case-study filtering.
#
# ============================================================

# ----------------------------
# Settings
# ----------------------------

RANK_FILE       <- here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv")
CONTRIBUTION_FILE <- here::here("Ablation_Table", "results", "tables", "RQ2A_contribution_shares.csv")
LOCO_FILE       <- here::here("Ablation_Table", "results", "tables", "RQ2B_LOCO_without_FC_delta_r_full.csv")

MERGED_OUTPUT   <- here::here("Ablation_Table", "data","candidate_pool_merged.csv")
AUDIT_OUTPUT    <- here::here("Ablation_Table", "results", "tables","case7_audit_results.csv")
SHORTLIST_OUTPUT <- here::here("Ablation_Table", "results", "tables","discordant_candidates.csv")

# Floating-point tolerance for values that are mathematically
# zero but may be represented as tiny floating-point numbers.
# This is NOT a scientific/biological threshold.
ZERO_TOL <- 1e-12


# ----------------------------
# Helpers
# ----------------------------

stop_msg <- function(...) {
  stop(paste0(...), call. = FALSE)
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop_msg(
      "File not found: ", path,
      "\nCheck the filename/path at the top of the script."
    )
  }
  
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop_msg(
      "Could not read ", path, ": ", e$message
    )
  )
}

require_columns <- function(dat, cols, file_label) {
  missing_cols <- setdiff(cols, names(dat))
  if (length(missing_cols) > 0) {
    stop_msg(
      file_label, " is missing column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

check_unique_key <- function(dat, key, file_label) {
  if (any(is.na(dat[[key]]) | dat[[key]] == "")) {
    stop_msg(file_label, ": missing/empty identifier in ", key)
  }
  
  if (anyDuplicated(dat[[key]])) {
    dup <- unique(dat[[key]][duplicated(dat[[key]])])
    stop_msg(
      file_label, ": duplicate identifiers in ", key,
      ": ", paste(dup, collapse = ", ")
    )
  }
}

as_numeric_checked <- function(x, column_name) {
  if (is.numeric(x)) return(x)
  
  converted <- suppressWarnings(as.numeric(x))
  
  if (any(is.na(converted) & !is.na(x))) {
    stop_msg(
      "Column '", column_name,
      "' contains non-numeric values."
    )
  }
  
  converted
}

classify_contribution <- function(c_fc, c_deg, q1, zero_tol = ZERO_TOL) {
  
  fc_zero  <- abs(c_fc) <= zero_tol
  deg_zero <- abs(c_deg) <= zero_tol
  
  if (!fc_zero && deg_zero) {
    return("FC-dominant")
  }
  
  if (fc_zero && !deg_zero) {
    return("network-dominant")
  }
  
  if (!fc_zero && !deg_zero &&
      abs(c_fc - c_deg) <= q1) {
    return("convergent")
  }
  
  if (fc_zero && deg_zero) {
    return("no_FC_or_Degree_contribution")
  }
  
  return("mixed_nonconvergent")
}


# ----------------------------
# Read the three sources
# ----------------------------

rank_dat <- read_required_csv(RANK_FILE)
contrib_dat <- read_required_csv(CONTRIBUTION_FILE)
loco_dat <- read_required_csv(LOCO_FILE)

require_columns(
  rank_dat,
  c("UNIPROT", "Symbol", "CDS", "rank_AHPCDS",
    "rank_FConly", "rank_Centralityonly"),
  RANK_FILE
)

require_columns(
  contrib_dat,
  c("UNIPROT", "C_FC", "C_Deg", "rank_AHP"),
  CONTRIBUTION_FILE
)

require_columns(
  loco_dat,
  c("id", "rank_AHP", "rank_baseline", "delta_r"),
  LOCO_FILE
)

check_unique_key(rank_dat, "UNIPROT", RANK_FILE)
check_unique_key(contrib_dat, "UNIPROT", CONTRIBUTION_FILE)
check_unique_key(loco_dat, "id", LOCO_FILE)


# ----------------------------
# Select source columns
# ----------------------------

rank_keep <- rank_dat[
  c("UNIPROT", "Symbol", "CDS", "rank_AHPCDS",
    "rank_FConly", "rank_Centralityonly")
]

contrib_keep <- contrib_dat[
  c("UNIPROT", "C_FC", "C_Deg", "rank_AHP")
]

loco_keep <- loco_dat[
  c("id", "rank_AHP", "rank_baseline", "delta_r")
]

names(loco_keep)[names(loco_keep) == "id"] <- "UNIPROT"


# ----------------------------
# Merge shuffled sources
# ----------------------------

merged <- merge(
  rank_keep,
  contrib_keep,
  by = "UNIPROT",
  all = TRUE,
  sort = FALSE
)

merged <- merge(
  merged,
  loco_keep,
  by = "UNIPROT",
  all = TRUE,
  suffixes = c("_contrib", "_loco"),
  sort = FALSE
)

expected_n <- nrow(rank_keep)

if (nrow(merged) != expected_n) {
  stop_msg(
    "Merged candidate count is ", nrow(merged),
    " but source 1 contains ", expected_n,
    ". Check identifiers and source files."
  )
}

# No source should have lost candidates.
if (any(is.na(merged$Symbol)) ||
    any(is.na(merged$CDS)) ||
    any(is.na(merged$rank_FConly)) ||
    any(is.na(merged$rank_Centralityonly)) ||
    any(is.na(merged$C_FC)) ||
    any(is.na(merged$C_Deg)) ||
    any(is.na(merged$rank_AHP_contrib)) ||
    any(is.na(merged$rank_AHP_loco)) ||
    any(is.na(merged$rank_baseline)) ||
    any(is.na(merged$delta_r))) {
  
  stop_msg(
    "The merge produced missing values. ",
    "This usually indicates identifiers that do not match across sources."
  )
}


# ----------------------------
# Validate cross-source fields
# ----------------------------

if (any(merged$rank_AHP_contrib != merged$rank_AHP_loco)) {
  bad <- merged$UNIPROT[
    merged$rank_AHP_contrib != merged$rank_AHP_loco
  ]
  
  stop_msg(
    "rank_AHP mismatch between contribution and LOCO files: ",
    paste(bad, collapse = ", ")
  )
}

if (any(abs(merged$CDS - merged$CDS[
  match(merged$UNIPROT, rank_dat$UNIPROT)
]) > ZERO_TOL)) {
  stop_msg("Unexpected CDS mismatch after merge.")
}

# Use the AHP rank from source 1 as the canonical AHP-CDS rank.
# Keep source-2/source-3 rank_AHP only for validation.
merged$rank_AHP <- merged$rank_AHPCDS
merged$rank_baseline_LOCO <- merged$rank_baseline
merged$delta_r_no_fc <- merged$delta_r

merged <- merged[
  c("UNIPROT", "Symbol", "CDS", "rank_AHPCDS",
    "rank_FConly", "rank_Centralityonly",
    "rank_AHP", "rank_baseline_LOCO", "delta_r_no_fc",
    "C_FC", "C_Deg")
]

write.csv(
  merged,
  MERGED_OUTPUT,
  row.names = FALSE,
  quote = TRUE
)


# ----------------------------
# Numeric validation
# ----------------------------

numeric_cols <- c(
  "CDS",
  "rank_AHPCDS",
  "rank_FConly",
  "rank_Centralityonly",
  "rank_AHP",
  "rank_baseline_LOCO",
  "delta_r_no_fc",
  "C_FC",
  "C_Deg"
)

for (col in numeric_cols) {
  merged[[col]] <- as_numeric_checked(merged[[col]], col)
}

if (any(is.na(merged[, numeric_cols]))) {
  stop_msg("Missing numeric values remain after merge.")
}

if (any(merged$C_FC < -ZERO_TOL) ||
    any(merged$C_Deg < -ZERO_TOL)) {
  stop_msg("Negative C_FC or C_Deg detected. Check contribution input.")
}


# ----------------------------
# Baseline classification
# ----------------------------
# Baseline is Centrality-only relative to FC-only:
#
#   baseline_delta_rank =
#       rank_Centralityonly - rank_FConly
#
# Smaller rank = better.
#
# > 0 : Centrality-only rank is worse than FC-only -> Dropped
# < 0 : Centrality-only rank is better than FC-only -> Rose
# = 0 : unchanged

merged$baseline_delta_rank <-
  merged$rank_Centralityonly - merged$rank_FConly

merged$baseline_direction <- ifelse(
  merged$baseline_delta_rank > 0,
  "Dropped",
  ifelse(
    merged$baseline_delta_rank < 0,
    "Rose",
    "Unchanged"
  )
)


# ----------------------------
# LOCO classification
# ----------------------------
#
# delta_r_no_fc is supplied directly from the LOCO dataset.
#
# > 0 : removing FC makes rank worse -> FC-rescued
# < 0 : removing FC makes rank better -> network-rescued
# = 0 : zero_change
#
# No near-neutral threshold is used.

merged$loco_direction <- ifelse(
  merged$delta_r_no_fc > 0,
  "FC-rescued",
  ifelse(
    merged$delta_r_no_fc < 0,
    "network-rescued",
    "zero_change"
  )
)


# ----------------------------
# Discordance
# ----------------------------

merged$discordant <- (
  (merged$baseline_direction == "Dropped" &
     merged$loco_direction == "network-rescued") |
    (merged$baseline_direction == "Rose" &
       merged$loco_direction == "FC-rescued")
)

merged$discordance_type <- ifelse(
  merged$baseline_direction == "Dropped" &
    merged$loco_direction == "network-rescued",
  "Dropped + network-rescued",
  ifelse(
    merged$baseline_direction == "Rose" &
      merged$loco_direction == "FC-rescued",
    "Rose + FC-rescued",
    NA_character_
  )
)


# ----------------------------
# Contribution classification
# ----------------------------

merged$contribution_difference <-
  abs(merged$C_FC - merged$C_Deg)

# IMPORTANT:
# Q1 is calculated across the FULL candidate pool.
q1_contribution_difference <- as.numeric(
  quantile(
    merged$contribution_difference,
    probs = 0.25,
    na.rm = FALSE,
    type = 7,
    names = FALSE
  )
)

merged$contribution_class <- vapply(
  seq_len(nrow(merged)),
  function(i) {
    classify_contribution(
      merged$C_FC[i],
      merged$C_Deg[i],
      q1_contribution_difference
    )
  },
  character(1)
)


# ----------------------------
# Manual-inspection ordering
# ----------------------------
# This is ONLY an ordering aid.
# It is NOT an automatic case-study selection.

merged$manual_inspection_priority <- NA_integer_
discordant_idx <- which(merged$discordant)

if (length(discordant_idx) > 0) {
  
  ordered_idx <- discordant_idx[
    order(
      -abs(merged$delta_r_no_fc[discordant_idx]),
      -merged$contribution_difference[discordant_idx],
      merged$Symbol[discordant_idx],
      merged$UNIPROT[discordant_idx]
    )
  ]
  
  merged$manual_inspection_priority[ordered_idx] <-
    seq_along(ordered_idx)
}


# ----------------------------
# Write audit outputs
# ----------------------------

write.csv(
  merged,
  AUDIT_OUTPUT,
  row.names = FALSE,
  quote = TRUE
)

discordant <- merged[merged$discordant, , drop = FALSE]

if (nrow(discordant) > 0) {
  discordant <- discordant[
    order(discordant$manual_inspection_priority),
    ,
    drop = FALSE
  ]
}

write.csv(
  discordant,
  SHORTLIST_OUTPUT,
  row.names = FALSE,
  quote = TRUE
)


# ----------------------------
# Console report
# ----------------------------

cat("\n")
cat("============================================================\n")
cat("AHP-CDS CASE-STUDY SELECTION AUDIT\n")
cat("============================================================\n\n")

cat("Merged candidate pool size: ", nrow(merged), "\n", sep = "")

cat(
  "Q1 of |C_FC - C_Deg| across FULL candidate pool: ",
  format(q1_contribution_difference, digits = 8),
  "\n\n",
  sep = ""
)

cat("Contribution classification:\n")
print(table(merged$contribution_class))
cat("\n")

cat("Baseline direction (Centrality-only vs FC-only):\n")
print(table(merged$baseline_direction))
cat("\n")

cat("LOCO direction (precomputed Delta r without FC):\n")
print(table(merged$loco_direction))
cat("\n")

cat("Discordant candidates: ", sum(merged$discordant), "\n", sep = "")

if (sum(merged$discordant) == 0) {
  
  cat("\nNo discordant candidates were found.\n")
  cat("This script does NOT recommend adding Case Study 7.\n")
  
} else {
  
  cat("\nDiscordance types:\n")
  print(table(merged$discordance_type, useNA = "ifany"))
  cat("\n")
  
  cat("Candidates for MANUAL INSPECTION (not automatically selected):\n\n")
  
  display_cols <- c(
    "manual_inspection_priority",
    "UNIPROT",
    "Symbol",
    "baseline_direction",
    "loco_direction",
    "discordance_type",
    "delta_r_no_fc",
    "C_FC",
    "C_Deg",
    "contribution_difference",
    "contribution_class"
  )
  
  print(
    discordant[, display_cols, drop = FALSE],
    row.names = FALSE
  )
}

cat("\n")
cat("Output files:\n")
cat("  - ", MERGED_OUTPUT, "\n", sep = "")
cat("  - ", AUDIT_OUTPUT, "\n", sep = "")
cat("  - ", SHORTLIST_OUTPUT, "\n", sep = "")

cat("\n")
cat("IMPORTANT:\n")
cat("  * Source rows may be shuffled; merging is by UNIPROT/id.\n")
cat("  * Cross-source rank_AHP consistency is checked.\n")
cat("  * No near-neutral threshold is used.\n")
cat("  * No candidate is selected automatically as Case Study 7.\n")
cat("  * Manual researcher judgment remains the final selection step.\n")
cat("============================================================\n\n")

