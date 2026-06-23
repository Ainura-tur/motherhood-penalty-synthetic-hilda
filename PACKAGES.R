# =============================================================================
# PACKAGES.R
#
# Installs and version-checks every package the analysis pipeline needs.
# Run this ONCE before run_all.R (or before sourcing MASTER_hh.R directly).
# The master script fails fast with an instruction to run this file if a
# required package is missing, so this is the first step of any fresh setup.
#
# Reference environment (from the run log captured on 2026-06-15):
#   R version 4.5.2 (2025-10-31)
# The versions pinned below are those used to produce the published results.
# Newer versions will usually work; the pins document what was actually run.
# =============================================================================

# ── Namespace protection (avoid MASS/stats masking dplyr verbs) ───────────────
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

# ── CRAN dependencies (package = version used in the reference run) ───────────
# tidyverse pulls dplyr, ggplot2, tidyr, purrr, readr, tibble, stringr,
# forcats, lubridate and scales, so those are not listed separately.
cran_pkgs <- c(
  tidyverse    = "2.0.0",   # data wrangling + ggplot2 + scales
  fixest       = "0.13.2",  # high-dimensional FE estimation, sunab event study
  lmtest       = "0.9-40",  # coefficient tests
  sandwich     = "3.1-1",   # robust / clustered covariance
  patchwork    = "1.3.2",   # figure composition
  DoubleML     = "1.0.2",   # double / debiased ML
  mlr3         = "1.2.0",   # learner backend for DoubleML
  mlr3learners = "0.13.0",  # ranger / glmnet learners for mlr3
  ranger       = "0.17.0",  # random forest nuisance learner
  glmnet       = "4.1-10",  # LASSO nuisance learner
  grf          = "2.5.0",   # causal forest
  ivreg        = "0.6-5",   # IV regression (required by ivcrtest)
  car          = "3.1-3",   # linear hypothesis tests
  AER          = "1.2-15",  # applied econometrics helpers
  MASS         = "7.3-65",  # (loaded transitively; pinned for reproducibility)
  future       = "1.67.0",  # parallel backend
  haven        = "2.5.5"    # read_dta() for the raw HILDA waves (00_load_hilda.R)
)

repos <- "https://cloud.r-project.org"

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing '%s' ...", pkg))
    install.packages(pkg, repos = repos)
  }
}

invisible(lapply(names(cran_pkgs), install_if_missing))

# ── Version report (warn, do not stop, on a mismatch) ─────────────────────────
cat("\nDependency check (reference run used R 4.5.2):\n")
cat(sprintf("  R: installed %s\n", getRversion()))
for (pkg in names(cran_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-13s MISSING (install failed)\n", pkg))
    next
  }
  have <- as.character(utils::packageVersion(pkg))
  want <- cran_pkgs[[pkg]]
  flag <- if (identical(have, want)) "ok" else sprintf("ref %s", want)
  cat(sprintf("  %-13s %-10s %s\n", pkg, have, flag))
}

# ── ivcrtest: the Contextual Robustness (CR) test (Ichimura 2025) ─────────────
# ivcrtest is NOT on CRAN. It is the custom package that provides
# ivcrtest::iv_cr_test(), used in PART C of MASTER_hh.R. The reference run
# used ivcrtest 0.1.0, which depends on ivreg (installed above).
#
# Install it from its source before running the analysis, for example:
#   # remotes::install_github("<owner>/ivcrtest")        # if hosted on GitHub
#   # install.packages("ivcrtest_0.1.0.tar.gz", repos = NULL, type = "source")
# Replace the source above with the actual location of the package.
if (!requireNamespace("ivcrtest", quietly = TRUE)) {
  cat("\n  ivcrtest      MISSING (custom package, not on CRAN)\n")
  cat("    Install ivcrtest 0.1.0 from source before running the analysis;\n")
  cat("    see the note at the bottom of PACKAGES.R. PART C (the CR test)\n")
  cat("    cannot run without it.\n")
} else {
  cat(sprintf("  %-13s %-10s %s\n", "ivcrtest",
              as.character(utils::packageVersion("ivcrtest")), "ref 0.1.0"))
}

cat("\nPACKAGES.R complete.\n")
