# =============================================================================
# config.R
#
# Central path configuration for the replication package. Sourced by
# 00_load_hilda.R (and available to the other scripts). Keeping every path in
# one file is what lets the package run on a machine other than the authors'.
#
# Run contract: run all scripts from the directory that contains them, so the
# bare-name reads and writes (for example readRDS("hilda_panel_data_*.rds"))
# resolve consistently. run_all.R sets that working directory automatically.
# =============================================================================

# ── Namespace protection (avoid MASS/stats masking dplyr verbs) ───────────────
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag
between <- dplyr::between   # masked by a later-attached package (e.g. ivcrtest)

# ── Data mode: which panels the pipeline runs on ──────────────────────────────
# "synthetic" (default): use the synthetic panels written by
#   00b_make_synthetic_hilda.R, so the package runs end to end without the
#   restricted HILDA release. Results on synthetic data are illustrative and are
#   NOT the paper's estimates.
# "real": use the panels built from the restricted HILDA release by
#   00_load_hilda.R. Only the authors and approved replicators can use this.
# Override without editing this file via the DATA_MODE environment variable,
# e.g.  DATA_MODE=real Rscript run_all.R
DATA_MODE <- tolower(Sys.getenv("DATA_MODE", unset = "synthetic"))
if (!DATA_MODE %in% c("synthetic", "real"))
  stop(sprintf("config.R: DATA_MODE must be 'synthetic' or 'real', got '%s'.",
               DATA_MODE))

# ── Raw HILDA data (restricted; not in the repository) ────────────────────────
# Folder holding the raw Combined_<wave><release>0c.dta wave files. Only stage
# 00 (00_load_hilda.R) needs this. Set it in one of two ways:
#   1. the HILDA_DATA_DIR environment variable (preferred; nothing to edit), or
#   2. by placing the .dta files in data/raw/ under the working directory.
data_dir <- Sys.getenv("HILDA_DATA_DIR",
                       unset = file.path(getwd(), "data", "raw"))

# ── Derived panels and results ────────────────────────────────────────────────
# The loader writes, and the analysis scripts read, the panels by bare name in
# the working directory, so the analysis directory is the working directory.
analysis_dir <- getwd()

# ── Output tree for tables and figures ────────────────────────────────────────
output_dir <- file.path(getwd(), "output")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
for (d in c(output_dir, fig_dir, tab_dir))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ── Sanity message (not fatal) ────────────────────────────────────────────────
# Only stage 00 (real mode) needs the raw data; the rest of the pipeline runs
# from the derived panels. A missing raw directory matters only in real mode,
# and even then is a warning, not an error.
if (DATA_MODE == "real" && !dir.exists(data_dir))
  message(sprintf(
    "config.R: raw data directory not found at '%s'.\n  Set HILDA_DATA_DIR or place the .dta files there before running 00_load_hilda.R.",
    data_dir))

cat(sprintf("config.R: DATA_MODE = %s\n", DATA_MODE))
if (DATA_MODE == "real") cat(sprintf("config.R: data_dir = %s\n", data_dir))