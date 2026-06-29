# =============================================================================
# run_all.R
#
# Single entry point for the analysis pipeline. Run from the repository root,
# either interactively or as:
#
#     Rscript run_all.R
#
# Stage 00 (sample construction) is run separately and ONCE only when you have
# the restricted HILDA release. In that case set DATA_MODE=real and build the
# three analysis panels first:
#
#     Rscript 00_load_hilda.R          # real data: builds hilda_panel_data_*.rds
#
# By default DATA_MODE is "synthetic": this script builds synthetic stand-in
# panels automatically (via 00b_make_synthetic_hilda.R) so a fresh clone runs
# end to end without the restricted data. Synthetic results are illustrative,
# not the paper's estimates.
#
# This script then runs, in ONE R session, in this order:
#
#   1. PACKAGES.R        install / version-check all dependencies
#   2. config.R          DATA_MODE and output paths
#   2b. data step        ensure the three panels exist (build synthetic if not)
#   3. 01_master.R       Parts A-F inline, then it internally sources
#                          02b_part_o_aligned.R      (estimator-aligned betaX,
#                                                     Table A4 / Table 13)
#                          02_cs_event_study.R       (child-penalty event study)
#                          07b_implied_dynamic_cost_concave.R (late-horizon
#                                                     return x gap share, App. A6)
#                          08_gelbach_decomposition.R (experience-channel decomposition,
#                                                      Table 10 / Appendix A6)
#                          04_part_g_ftbb_iv.R       (FTB-B reform IV, Part G)
#                          05_part_h_policy.R        (policy economics, Part H)
#                          03_training_event_study.R (training event study)
#                        The internal order is deliberate: 07b reads the event
#                        study and 08 reads 07b, so 02 -> 07b -> 08 run in that
#                        sequence; Part H needs Part G's in-memory objects and the
#                        CS results file, so CS and the event studies run before Part H.
#   4. 06_figures.R      all publication figures; reuses 01_master.R's in-memory
#                        objects, so it MUST share this session.
#
# The scripts read and write by bare filename, so the working directory must be
# the folder that holds the scripts and the hilda_panel_data_*.rds inputs. This
# file sets that working directory to its own location automatically.
# =============================================================================

# ── Resolve this file's directory and use it as the working directory ─────────
this_file <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)        # works when run via Rscript
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1]) else NA_character_
}, error = function(e) NA_character_)

if (is.na(this_file) && !is.null(sys.frames()[[1]]$ofile)) {
  this_file <- sys.frames()[[1]]$ofile             # works when source()d
}

repo_root <- if (!is.na(this_file)) normalizePath(dirname(this_file)) else getwd()
setwd(repo_root)
cat(sprintf("Working directory: %s\n", getwd()))

# ── Stage 1: dependencies (first, so config.R and the data step can use dplyr) ─
cat("\n==== Stage 1/3: PACKAGES.R ====\n")
source("PACKAGES.R")

# ── Configuration: DATA_MODE and output paths ─────────────────────────────────
source("config.R")   # sets DATA_MODE ("synthetic" default | "real") and paths

# ── Data: ensure the three analysis panels are present, per DATA_MODE ─────────
# In "synthetic" mode (the default) the panels are generated on the fly by
# 00b_make_synthetic_hilda.R, so a fresh clone runs end to end with no manual
# step and without the restricted HILDA release. In "real" mode the panels must
# already have been built by 00_load_hilda.R from the restricted release.
required_data <- c(
  "hilda_panel_data_W12_W24_slim.rds",  # Parts A-F
  "hilda_panel_data_extended.rds",      # Parts G/H
  "hilda_panel_data_training.rds"       # training event study
)
missing_data <- required_data[!file.exists(required_data)]

if (length(missing_data)) {
  if (identical(DATA_MODE, "synthetic")) {
    cat("\n==== Building SYNTHETIC panels (DATA_MODE = 'synthetic') ====\n")
    cat("Results on synthetic data are illustrative, NOT the paper's estimates.\n")
    source("00b_make_synthetic_hilda.R")
  } else {
    stop(paste0(
      "DATA_MODE = 'real', but analysis panel(s) are missing:\n  ",
      paste(missing_data, collapse = "\n  "),
      "\nBuild them first with:  Rscript 00_load_hilda.R\n",
      "These are derived from the restricted HILDA release and are not in the ",
      "repository.\nSee Section 2 (Data availability) of README.md.\n",
      "To run without the restricted data instead, use DATA_MODE=synthetic."
    ))
  }
  still_missing <- required_data[!file.exists(required_data)]
  if (length(still_missing))
    stop("Panels still missing after the build step:\n  ",
         paste(still_missing, collapse = "\n  "))
}

# ── Stage 2: full analysis (Parts A-H; 01 sources 02-05 in dependency order) ──
cat("\n==== Stage 2/3: 01_master.R ====\n")
source("01_master.R")

# ── Stage 3: figures (same session, reuses 01_master.R objects) ───────────────
cat("\n==== Stage 3/3: 06_figures.R ====\n")
source("06_figures.R")

cat("\nPipeline complete. Figures are in output/figures/.\n")

# ── Result-completeness self-check ───────────────────────────────────────────
# Announce, in the log, whether every expected result object was produced. This
# turns a silent partial run into an explicit PASS/MISSING line per result.
.expected <- c(
  "cs_event_study_results.rds",                  # event study (P2)
  "implied_dynamic_cost_concave_results.rds",    # 07b late-horizon share
  "gelbach_decomposition_results.rds",           # Table 10 / App. A6
  "ftbb_reform1_iv_results.rds",                 # Part G
  "policy_economics_results.rds",                # Part H
  "training_event_study_results.rds",            # training
  "robustness_master.rds",                        # robustness battery
  "betaX_lee_bounds.csv",                         # Part O / Table 13
  "se_ledger_A3_aligned.csv"                      # Part O / Table A4
)
cat("\n==== Result-completeness check ====\n")
.missing <- character(0)
for (f in .expected) {
  ok <- file.exists(f)
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "MISSING", f))
  if (!ok) .missing <- c(.missing, f)
}
if (length(.missing)) {
  cat(sprintf("\n  %d expected result(s) MISSING — the run is PARTIAL.\n", length(.missing)))
} else {
  cat("\n  All expected results present.\n")
}
