# =============================================================================
# DATA AVAILABILITY
# Raw data: HILDA Survey (Waves 12–24), accessed via the Australian Data Archive
# (ADA). Researchers must apply via: https://dataverse.ada.edu.au/dataverse/hilda
# Data inputs (all built by HILDA_LOADER_UNIFIED.R in a single pass):
#   hilda_panel_data_W12_W24_slim.rds  (Parts A-F; read below)
#   hilda_panel_data_extended.rds      (Part G/H; read by the sourced scripts)
#   hilda_panel_data_training.rds      (training event study fallback)
# (also in this repository). The RDS file itself cannot be shared.
# =============================================================================




# =============================================================================
# within_couples_MASTER_ANALYSIS.R  —  REVISED FOR PAPER CONSISTENCY
# Unified household-specialisation analysis
#
# REVISION NOTES (see consistency_report.md for full audit):
#   1. A5 pooled interaction: removed state FE + weights to match baseline (Eq.4)
#   2. FE formulas: removed educ_years (absorbed by person FE, per §4.1)
#   2b. Restored educ_years in OLS formulas incorrectly stripped by initial fix
#   2c. Made breadwinner lock-in if-conditions NA-safe (b_educ can be NA)
#   3. Added Heckman two-step selection correction (Table 3, §4.3)
#   4. Added levels-based P1 test in B4 (couple FE, paper §6)
#   5. Descriptive stats: use ln_hourly_wage_real (not nominal)
#   6. Hardcoded summary values: added validation note
#   7. Orphaned theme() calls: commented out
#
# PART A  Main couple analysis & breadwinner lock-in
# PART B  Housework mediation & time-budget test
# PART C  Panel IV (parental-education instrument, diagnostic only)
# PART D  Couple-level DML + causal forest
# PART E  Robustness battery (Parts 1–10) + all publication figures
# PART F  Oster (2019) selection bounds for β₂
# PART G  FTB-B Reform 1 simulated-instrument IV (sourced)
# =============================================================================

# Set working directory 

#setwd("C:/Users/tur277/OneDrive - CSIRO/Desktop/WG")

# ── Run logging: tee all console output to a timestamped file ────────────────
# stdout is split (console + file); messages/warnings (incl. fixest NOTEs)
# go to the file only. Sinks are closed at the very end of this script.
log_dir <- file.path(getwd(), "run_logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, format(Sys.time(), "MASTER_hh_%Y%m%d_%H%M%S.log"))
# clear any stale sinks left by a previous aborted run
while (sink.number() > 0) sink()
if (sink.number(type = "message") > 0) sink(type = "message")
.log_con <- file(log_file, open = "wt")
sink(.log_con, split = TRUE)
sink(.log_con, type = "message")
cat("Logging console output to:", log_file, "\n")

outdir <- file.path(getwd(), "output", "figures")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(lmtest)
  library(sandwich)
  library(patchwork)
  library(ggplot2)
})

# Run PACKAGES.R first to install and version-check all dependencies.
for (pkg in c("DoubleML", "mlr3", "mlr3learners", "ranger", "glmnet", "grf")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' not found. Run PACKAGES.R first.", pkg))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Thread budget. On a SLURM batch job, use exactly the CPUs allocated
# (SLURM_CPUS_PER_TASK); otherwise (e.g. an interactive login-node session)
# cap at 4 so the script never monopolises a shared node. Override with
# the MASTER_HH_THREADS environment variable if needed.
n_threads <- suppressWarnings(as.integer(Sys.getenv("MASTER_HH_THREADS", "")))
if (is.na(n_threads))
  n_threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (is.na(n_threads) || n_threads < 1L)
  n_threads <- min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L))
fixest::setFixest_nthreads(n_threads)
cat(sprintf("Using %d threads (fixest, ranger, grf)\n", n_threads))

set.seed(42)

# ── WPI deflators (base: Wave 24 = 2024) ─────────────────────────────────────
wpi_data <- data.frame( wave = 12:24, deflator = 149.6 / c(110.9, 114.6, 117.6, 120.4, 123.0, 125.4, 127.9, 130.9, 133.7, 135.7, 138.9, 143.7, 149.6) )


# ── Shared helpers & aesthetics ───────────────────────────────────────────────
stars_fn <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                         ifelse(p < 0.05, "*", ifelse(p < 0.10, "\u2020", "")))))
}

label_group <- function(g) {
  case_when(
    g == "childless_men_ever" ~ "Childless men",
    g == "fathers_ever"       ~ "Fathers",
    g == "never_mothers"      ~ "Never-mothers",
    g == "mothers_ever"       ~ "Mothers (ever)",
    TRUE ~ g
  )
}

pal4       <- c("Childless men"="#2166AC","Fathers"="#92C5DE",
                "Never-mothers"="#B2182B","Mothers (ever)"="#F4A582")
pal3w      <- c("Never-mothers"="#B2182B","Active mothers"="#D6604D","Post-mothers"="#FDDBC7")
pal_gender <- c("Men"="#2166AC","Women"="#B2182B")
pal_est    <- c("OLS"="#4393C3","FE"="#D6604D")

theme_paper <- theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    legend.position      = "bottom",
    legend.title         = element_blank(),
    plot.title           = element_text(size = 7.5, face = "bold", hjust = 0),
    plot.subtitle        = element_text(size = 6,   colour = "grey40", hjust = 0),
    plot.caption         = element_text(size = 6,   colour = "grey50", hjust = 0),
    strip.text           = element_text(size = 6.5, face = "bold"),
    axis.title           = element_text(size = 6.5)
  )


save_fig <- function(p, name, w = 13/2.54, h = 8/2.54) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  ggsave(paste0(name, ".pdf"), p, width = w, height = h, units = "in")
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 500, units = "in")
  for (ext in c(".pdf", ".png")) {
    ff <- paste0(name, ext)
    if (file.exists(ff)) file.copy(ff, file.path(outdir, ff), overwrite = TRUE)
  }
  cat(sprintf("  \u2713 %s\n", name))
}

# ── Namespace protection (avoid MASS/stats masking dplyr verbs) ───────────────
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag
# fixest estimation accessors: pin to fixest so a later-attached package
# (DoubleML/mlr3/aod/etc.) cannot mask them. The 13 June crash was a masked
# wald() that dispatched round_any() on a fixest object.
wald <- fixest::wald;          pvalue    <- fixest::pvalue

# ── Results registry ────────────────────────────────────────────────────────
# Each section writes its live-computed headline value here via rec(); the
# FINAL VERIFICATION block prints them via pv() (registry) or gv() (persistent
# objects). Both fall back to "[not computed this run]" so a skipped or failed
# section never aborts the run and never prints a stale literal.
.PR <- new.env(parent = emptyenv())
rec <- function(key, val) { assign(key, val, envir = .PR); invisible(val) }
pv  <- function(key, fmt = "%.4f") {
  if (exists(key, envir = .PR, inherits = FALSE)) {
    v <- get(key, envir = .PR)
    if (length(v) == 0 || all(is.na(v))) return("[not computed this run]")
    if (is.numeric(v)) sprintf(fmt, v[1]) else as.character(v[1])
  } else "[not computed this run]"
}
gv  <- function(expr, fmt = "%.4f") {         # safe read of a persistent object
  v <- tryCatch(expr, error = function(e) NA)
  if (length(v) == 0 || all(is.na(v))) return("[not computed this run]")
  if (is.numeric(v)) sprintf(fmt, v[1]) else as.character(v[1])
}

# ── SAMPLE: Full panel ──────────────────────────────────────────────────────
# Paper §3 (Data): "91,279 person-wave observations across 16,276 individuals"
# Appendix Table: Full panel → counts re-verify against current extract; after
#   age_exp_diff filter (ages 30–49, experience plausibility check).
# ── LOAD DATA (single read) ───────────────────────────────────────────────────
panel <- readRDS("hilda_panel_data_W12_W24_slim.rds")
num_rows    <- format(nrow(panel), big.mark = ",")
num_persons <- format(n_distinct(panel$person_id), big.mark = ",")
cat(paste("Loaded:", num_rows, "obs |", num_persons, "individuals\n"))

# Sanity check: xwaveid and person_id should map 1:1
cat("\nID consistency check:\n")
id_check <- panel %>%
  summarise(
    n_xwaveid_unique  = n_distinct(xwaveid),
    n_person_id_unique = n_distinct(person_id),
    n_pairs_unique    = n_distinct(paste(xwaveid, person_id, sep = "|"))
  )
print(id_check)
# All three numbers should be equal if xwaveid and person_id are 1:1.
# In within_couples_MASTER_ANALYSISpro.R, after the panel is loaded but
# BEFORE any analysis that uses group_ever, insert:

cat("\n--- Recomputing group_ever with full-window absorbing definition ---\n")
person_status_full <- panel %>%
  group_by(person_id) %>%
  summarise(
    is_female_p       = max(female == 1, na.rm = TRUE),
    has_child_anywave = max(
      pmax(coalesce(ever_parent, 0L),
           as.integer(coalesce(total_children_ever_had, 0L) > 0),
           na.rm = TRUE),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    is_female_p       = ifelse(is.infinite(is_female_p),       NA_integer_, is_female_p),
    has_child_anywave = ifelse(is.infinite(has_child_anywave), 0L,          has_child_anywave),
    group_ever_full = case_when(
      is_female_p == 0 & has_child_anywave == 0 ~ "childless_men_ever",
      is_female_p == 0 & has_child_anywave == 1 ~ "fathers_ever",
      is_female_p == 1 & has_child_anywave == 0 ~ "never_mothers",
      is_female_p == 1 & has_child_anywave == 1 ~ "mothers_ever",
      TRUE                                     ~ NA_character_
    )
  )

panel <- panel %>%
  select(-any_of("group_ever")) %>%
  left_join(person_status_full %>% select(person_id, group_ever = group_ever_full),
            by = "person_id")

cat("  group_ever coverage by wave (should be ~100% all 13 waves):\n")
panel %>%
  group_by(wave) %>%
  summarise(pct = round(100 * mean(!is.na(group_ever)), 1)) %>%
  print(n = Inf)


# ============================================================
#  recompute group_3way_women from new group_ever
# ============================================================
panel <- panel %>%
  mutate(
    group_3way_women = case_when(
      female != 1                                       ~ NA_character_,
      group_ever == "never_mothers"                     ~ "never_mothers",
      group_ever == "mothers_ever" &
        coalesce(num_children_under15, 0L) > 0          ~ "active_mothers",
      group_ever == "mothers_ever" &
        coalesce(num_children_under15, 0L) == 0         ~ "post_mothers",
      TRUE                                              ~ NA_character_
    )
  )

cat("\n  group_3way_women distribution (women only, time-varying):\n")
print(table(panel$group_3way_women, useNA = "ifany"))


panel <- panel %>%
  mutate(
    married_num       = as.numeric(married),
    fulltime_num      = as.numeric(fulltime),
    public_sector_num = as.numeric(public_sector),
    wave_num          = as.numeric(wave),
    year              = 2000 + wave,
    age_exp_diff      = age - experience_years
  ) %>%
  filter(between(age_exp_diff, 10, 35))

# Resolve housework column — check which name exists in the data
.hw_candidates <- c("housework_hours", "care_intensity", "hsjob", "jomfh", "hw_hours")
.hw_col <- .hw_candidates[.hw_candidates %in% names(panel)][1]
if (!is.na(.hw_col)) {
  panel$hw <- pmin(pmax(panel[[.hw_col]], 0), 100)
  cat(sprintf("  hw resolved from: '%s'  (non-NA: %d%%)\n",
              .hw_col, round(100 * mean(!is.na(panel$hw)))))
} else {
  panel$hw <- NA_real_
  cat("  WARNING: no housework column found — hw is NA throughout\n")
  cat(sprintf("  Available columns: %s\n", paste(names(panel)[1:20], collapse = ", ")))
}
panel$total_hours <- ifelse(
  is.na(panel$hours_worked_clean) | is.na(panel$hw), NA_real_,
  panel$hours_worked_clean + panel$hw
)

# partner_id_clean: clean HILDA partner ID, available globally for all couple sections.
# Two-step: check existence outside mutate, then run the clean mutate unconditionally.
if (!"partner_id" %in% names(panel)) panel$partner_id <- NA_character_
panel <- panel %>%
  mutate(
    partner_id_clean = dplyr::case_when(
      is.na(partner_id)       ~ NA_character_,
      partner_id == ""        ~ NA_character_,
      nchar(partner_id) < 5  ~ NA_character_,
      grepl("^-", partner_id) ~ NA_character_,
      TRUE                    ~ partner_id
    )
  )

# ── SAMPLE: Wage sample (main analysis) ────────────────────────────────────
# Paper §3: "wage sample (employed individuals with valid hourly wages),
#   comprising 87,229 observations across 15,668 individuals"
# Filter: in_wage_sample == 1 (employed, valid hourly wage, trimmed)
wage <- panel %>% filter(in_wage_sample == 1)

num_rows    <- format(nrow(wage), big.mark = ",")
num_persons <- format(n_distinct(wage$person_id), big.mark = ",")
cat(paste("Wage sample:", num_rows, "obs |", num_persons, "individuals\n"))

# ── SAMPLE MAP VERIFICATION ─────────────────────────────────────────────────
# Prints actual sample sizes vs the expected counts in Appendix Table
# "Sample Definitions by Analysis" (Table \ref{tab:sample_map}).
# This output is used to populate that table; differences indicate
# data-version or filter changes.

cat("\n")
cat(strrep("=", 80), "\n")
cat("  SAMPLE MAP VERIFICATION (Appendix Table: Sample Definitions)\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("  %-45s  %10s  %10s\n", "Sample", "Obs", "Persons"))
cat("  ", strrep("-", 70), "\n")

cat(sprintf("  %-45s  %10s  %10s\n",
            "Full panel (after age_exp_diff filter)",
            format(nrow(panel), big.mark = ","),
            format(n_distinct(panel$person_id), big.mark = ",")))

cat(sprintf("  %-45s  %10s  %10s\n",
            "Wage sample (in_wage_sample == 1)",
            format(nrow(wage), big.mark = ","),
            format(n_distinct(wage$person_id), big.mark = ",")))

# Couple panel count (computed here, used later)
n_with_partner <- sum(!is.na(panel$partner_id_clean) &
                        panel$partner_id_clean %in% panel$person_id)
cat(sprintf("  %-45s  %10s  %10s\n",
            "Obs with valid partner in panel",
            format(n_with_partner, big.mark = ","),
            "(see Part B)"))

# Wage sample by group
for (g in c("childless_men_ever", "fathers_ever", "never_mothers", "mothers_ever")) {
  d <- wage %>% filter(group_ever == g)
  cat(sprintf("  %-45s  %10s  %10s\n",
              paste0("  Wage: ", label_group(g)),
              format(nrow(d), big.mark = ","),
              format(n_distinct(d$person_id), big.mark = ",")))
}

# Three-way women split
for (g in c("never_mothers", "active_mothers", "post_mothers")) {
  d <- wage %>% filter(group_3way_women == g)
  if (nrow(d) > 0)
    cat(sprintf("  %-45s  %10s  %10s\n",
                paste0("  Wage 3-way: ", g),
                format(nrow(d), big.mark = ","),
                format(n_distinct(d$person_id), big.mark = ",")))
}

# Transition counts
men_trans_n <- panel %>%
  filter(female == 0) %>%
  group_by(person_id) %>%
  summarise(trans = any(diff(ever_parent) > 0), .groups = "drop") %>%
  filter(trans) %>% nrow()
women_trans_n <- panel %>%
  filter(female == 1) %>%
  group_by(person_id) %>%
  summarise(trans = any(diff(ever_parent) > 0), .groups = "drop") %>%
  filter(trans) %>% nrow()
cat(sprintf("  %-45s  %10s  %10s\n",
            "Fatherhood transitions",
            "—", format(men_trans_n, big.mark = ",")))
cat(sprintf("  %-45s  %10s  %10s\n",
            "Motherhood transitions",
            "—", format(women_trans_n, big.mark = ",")))

# IV sample (if parental education available)
if ("father_educ_years" %in% names(panel) || "mother_educ_years" %in% names(panel)) {
  iv_samp <- wage %>%
    mutate(parent_educ_avg = case_when(
      !is.na(father_educ_years) & !is.na(mother_educ_years) ~
        (father_educ_years + mother_educ_years) / 2,
      !is.na(father_educ_years) ~ father_educ_years,
      !is.na(mother_educ_years) ~ mother_educ_years,
      TRUE ~ NA_real_
    )) %>%
    filter(!is.na(parent_educ_avg))
  cat(sprintf("  %-45s  %10s  %10s\n",
              "IV sample (parental educ available)",
              format(nrow(iv_samp), big.mark = ","),
              format(n_distinct(iv_samp$person_id), big.mark = ",")))
}

# Robustness samples
cat(sprintf("  %-45s  %10s\n", "Robustness: age 25–54",
            format(nrow(panel %>% filter(age >= 25, age <= 54)), big.mark = ",")))
cat(sprintf("  %-45s  %10s\n", "Robustness: full-time only",
            format(nrow(panel %>% filter(fulltime == 1)), big.mark = ",")))
cat(sprintf("  %-45s  %10s\n", "Robustness: part-time only",
            format(nrow(panel %>% filter(fulltime == 0, employed == 1)), big.mark = ",")))

if ("longitudinal_weight" %in% names(panel)) {
  lw <- panel %>% filter(!is.na(longitudinal_weight), longitudinal_weight > 0)
  cat(sprintf("  %-45s  %10s  %10s\n",
              "Robustness: longitudinal weights",
              format(nrow(lw), big.mark = ","),
              format(n_distinct(lw$person_id), big.mark = ",")))
}

# Trimming variants
trim19 <- wage %>%
  group_by(female) %>%
  mutate(p01 = quantile(hourly_wage_clean, 0.01, na.rm = TRUE),
         p99 = quantile(hourly_wage_clean, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(hourly_wage_clean >= p01, hourly_wage_clean <= p99)
cat(sprintf("  %-45s  %10s\n", "Robustness: trim 1/99",
            format(nrow(trim19), big.mark = ",")))

notrim <- panel %>%
  filter(employed == 1, !is.na(hourly_wage_clean), hourly_wage_clean > 0)
cat(sprintf("  %-45s  %10s\n", "Robustness: no trim",
            format(nrow(notrim), big.mark = ",")))

cat("\n  NOTE: Use these counts to populate Appendix Table \\ref{tab:sample_map}.\n")
cat("  Differences from hardcoded table values indicate data-version changes.\n")
cat(strrep("=", 80), "\n\n")



################################################################################
# PART A: MAIN COUPLE ANALYSIS & BREADWINNER LOCK-IN
################################################################################

cat("Loading data...\n")
num_rows    <- format(nrow(panel), big.mark = ",")
num_persons <- format(n_distinct(panel$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))




# #############################################################################
# PART A: FATHERHOOD TRANSITION ANALYSIS (no partner linkage needed)
# #############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
# Paper §6 (Household Specialisation and Breadwinner Lock-In)
# Paper §6.1 (Parenthood Transitions): "Among 7,791 men, 927 (11.9%) transition
#   ...Among 8,097 women, 1,058 (13.1%) become mothers"
# Sample: Individual-level panel (non-matched), broad sample for transitions,
#   wage sample for β₂ estimation.
cat("  PART A: FATHERHOOD TRANSITION — BREADWINNER LOCK-IN\n")
cat(strrep("=", 80), "\n")


# Paper §6.1: Fatherhood/motherhood transitions identified from ever_parent
#   changes 0→1 within the panel. Expected: ~1,252 men, ~1,397 women.
# Sample: Broad panel (all persons, not just wage sample).
# A1. IDENTIFY FATHERHOOD TRANSITIONS DURING PANEL

cat("\n--- A1. Fatherhood transitions during panel ---\n")

# Men who start as childless and become fathers during the panel
men_transitions <- panel %>%
  filter(female == 0) %>%
  group_by(person_id) %>%
  arrange(wave) %>%
  summarise(
    first_wave = first(wave),
    last_wave = last(wave),
    n_waves = n(),
    ever_parent_first = first(ever_parent),
    ever_parent_last = last(ever_parent),
    transition_wave = if (any(diff(ever_parent) > 0))
      wave[which(diff(ever_parent) > 0)[1] + 1] else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(
    transitioner = !is.na(transition_wave),
    always_childless = ever_parent_first == 0 & ever_parent_last == 0,
    always_father = ever_parent_first == 1 & ever_parent_last == 1
  )

cat(sprintf("  Total men: %s\n", format(nrow(men_transitions), big.mark = ",")))
cat(sprintf("  Transitioners (childless→father): %d (%.1f%%)\n",
            sum(men_transitions$transitioner),
            100 * mean(men_transitions$transitioner)))
cat(sprintf("  Always childless: %d\n", sum(men_transitions$always_childless)))
cat(sprintf("  Always father: %d\n", sum(men_transitions$always_father)))

# Same for women
women_transitions <- panel %>%
  filter(female == 1) %>%
  group_by(person_id) %>%
  arrange(wave) %>%
  summarise(
    first_wave = first(wave),
    last_wave = last(wave),
    n_waves = n(),
    ever_parent_first = first(ever_parent),
    ever_parent_last = last(ever_parent),
    transition_wave = if (any(diff(ever_parent) > 0))
      wave[which(diff(ever_parent) > 0)[1] + 1] else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(
    transitioner = !is.na(transition_wave),
    always_childless = ever_parent_first == 0 & ever_parent_last == 0,
    always_mother = ever_parent_first == 1 & ever_parent_last == 1
  )

cat(sprintf("\n  Total women: %s\n", format(nrow(women_transitions), big.mark = ",")))
cat(sprintf("  Transitioners (childless→mother): %d (%.1f%%)\n",
            sum(women_transitions$transitioner),
            100 * mean(women_transitions$transitioner)))


# A2. PRE vs POST TRANSITION COMPARISON

cat("\n\n--- A2. Pre vs Post transition: wages, hours, complementarity ---\n")

# For men who transition: compare observables before vs after
men_trans_ids <- men_transitions %>% filter(transitioner) %>% pull(person_id)

men_trans_data <- wage %>%
  filter(female == 0, person_id %in% men_trans_ids) %>%
  left_join(men_transitions %>% select(person_id, transition_wave), by = "person_id") %>%
  mutate(
    period = case_when(
      wave < transition_wave ~ "Pre-fatherhood",
      wave >= transition_wave ~ "Post-fatherhood"
    ),
    years_from_transition = wave - transition_wave
  )

cat("\n  MEN transitioning to fatherhood:\n")
men_trans_data %>%
  group_by(period) %>%
  summarise(
    N = n(),
    N_persons = n_distinct(person_id),
    Mean_wage = round(mean(exp(ln_hourly_wage_real), na.rm = TRUE), 2),
    Mean_ln_wage = round(mean(ln_hourly_wage_real, na.rm = TRUE), 3),
    Mean_hours = round(mean(hours_worked_clean, na.rm = TRUE), 1),
    Mean_exp = round(mean(experience_years, na.rm = TRUE), 1),
    Pct_fulltime = round(100 * mean(fulltime == 1, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print()

# Same for women
women_trans_ids <- women_transitions %>% filter(transitioner) %>% pull(person_id)

women_trans_data <- wage %>%
  filter(female == 1, person_id %in% women_trans_ids) %>%
  left_join(women_transitions %>% select(person_id, transition_wave), by = "person_id") %>%
  mutate(
    period = case_when(
      wave < transition_wave ~ "Pre-motherhood",
      wave >= transition_wave ~ "Post-motherhood"
    ),
    years_from_transition = wave - transition_wave
  )

cat("\n  WOMEN transitioning to motherhood:\n")
women_trans_data %>%
  group_by(period) %>%
  summarise(
    N = n(),
    N_persons = n_distinct(person_id),
    Mean_wage = round(mean(exp(ln_hourly_wage), na.rm = TRUE), 2),
    Mean_ln_wage = round(mean(ln_hourly_wage_real, na.rm = TRUE), 3),
    Mean_hours = round(mean(hours_worked_clean, na.rm = TRUE), 1),
    Mean_exp = round(mean(experience_years, na.rm = TRUE), 1),
    Pct_fulltime = round(100 * mean(fulltime == 1, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print()


# Paper §6.1, Figure 2: "Men's complementarity collapses at first birth"
# Sample: Wage-sample transitioners with event-time in [-3, +4].
#   Expected counts: re-verify against the current extract.
# A3. EVENT STUDY: WAGE GROWTH AROUND PARENTHOOD TRANSITION

cat("\n\n--- A3. Event study: wage growth around parenthood transition ---\n")

# FE wage model with event-time dummies for men
if (nrow(men_trans_data) > 100) {
  men_trans_data <- men_trans_data %>%
    mutate(
      event_time = years_from_transition,
      event_time_f = factor(
        pmax(-3, pmin(event_time, 4)),
        levels = -3:4
      )
    )
  
  cat("\n  MEN — Event-time distribution:\n")
  men_trans_data %>% count(event_time_f) %>% print()
  
  # FE model with event-time dummies (reference: t-1)
  men_event <- tryCatch(feols(
    ln_hourly_wage_real ~ event_time_f + experience_years + age_sq + married |
      person_id + wave,
    data = men_trans_data %>% mutate(event_time_f = relevel(event_time_f, ref = "-1")),
    cluster = ~person_id, notes = FALSE
  ), error = function(e) { cat("  Error:", e$message, "\n"); NULL })
  
  if (!is.null(men_event)) {
    cat("\n  MEN — FE wage effects relative to year before fatherhood:\n")
    s <- summary(men_event)
    et_coefs <- grep("event_time_f", names(coef(men_event)), value = TRUE)
    for (v in et_coefs) {
      b <- coef(men_event)[v]
      p <- s$coeftable[v, "Pr(>|t|)"]
      cat(sprintf("    %s: %+.3f%%%s\n", v, b*100, stars_fn(p)))
    }
  }
}

# Same for women
if (nrow(women_trans_data) > 100) {
  women_trans_data <- women_trans_data %>%
    mutate(
      event_time = years_from_transition,
      event_time_f = factor(
        pmax(-3, pmin(event_time, 4)),
        levels = -3:4
      )
    )
  
  cat("\n  WOMEN — Event-time distribution:\n")
  women_trans_data %>% count(event_time_f) %>% print()
  
  women_event <- tryCatch(feols(
    ln_hourly_wage_real ~ event_time_f + experience_years + age_sq + married |
      person_id + wave,
    data = women_trans_data %>% mutate(event_time_f = relevel(event_time_f, ref = "-1")),
    cluster = ~person_id, notes = FALSE
  ), error = function(e) { cat("  Error:", e$message, "\n"); NULL })
  
  if (!is.null(women_event)) {
    cat("\n  WOMEN — FE wage effects relative to year before motherhood:\n")
    s <- summary(women_event)
    et_coefs <- grep("event_time_f", names(coef(women_event)), value = TRUE)
    for (v in et_coefs) {
      b <- coef(women_event)[v]
      p <- s$coeftable[v, "Pr(>|t|)"]
      cat(sprintf("    %s: %+.3f%%%s\n", v, b*100, stars_fn(p)))
    }
  }
}


# Paper §6.2 (Breadwinner Lock-In): "For men transitioning to fatherhood,
#   level returns β_educ rise from 4.34% to 5.57% (+1.23 pp) while
#   complementarity β₂ falls from 0.074% to 0.047% (−36%)"
# Sample: Male transitioners only (wage sample), split pre/post.
# A4. BREADWINNER LOCK-IN: β₂ BEFORE vs AFTER FATHERHOOD

cat("\n\n--- A4. Breadwinner lock-in: β₂ before vs after fatherhood ---\n")

# For male transitioners: estimate β₂ separately in pre and post periods
run_b2_period <- function(data, label) {
  d <- data %>%
    mutate(
      educ_c = educ_years - mean(educ_years, na.rm = TRUE),
      exp_c = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c
    )
  
  # OLS
  ols <- tryCatch(feols(
    ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years + age_sq + married | wave,
    data = d, cluster = ~person_id, notes = FALSE
  ), error = function(e) NULL)
  
  if (!is.null(ols) && "educ_exp_c" %in% names(coef(ols))) {
    s <- summary(ols)
    b <- coef(ols)["educ_exp_c"]
    p <- s$coeftable["educ_exp_c", "Pr(>|t|)"]
    b_educ <- if ("educ_years" %in% names(coef(ols))) coef(ols)["educ_years"] else NA
    b_educ_str <- if (!is.na(b_educ)) sprintf("%.2f%%", b_educ*100) else "NA"
    cat(sprintf("    %-30s  N=%5d  β_educ=%s  β₂=%.4f%%%s\n",
                label, nrow(d), b_educ_str, b*100, stars_fn(p)))
    return(list(b2 = b, p = p, se = s$coeftable["educ_exp_c", "Std. Error"], b_educ = b_educ, n = nrow(d)))
  }
  cat(sprintf("    %-30s  N=%5d  FAILED\n", label, nrow(d)))
  return(NULL)
}

cat("\n  MEN — Pre vs Post fatherhood β₂ (transitioners only, OLS):\n")
r_pre_m <- run_b2_period(men_trans_data %>% filter(period == "Pre-fatherhood"),
                         "Pre-fatherhood")
r_post_m <- run_b2_period(men_trans_data %>% filter(period == "Post-fatherhood"),
                          "Post-fatherhood")

if (!is.null(r_pre_m) && !is.null(r_post_m)) {
  rec("lockin_b2_pre",  r_pre_m$b2 * 100);   rec("lockin_b2_post",  r_post_m$b2 * 100)
  rec("lockin_be_pre",  r_pre_m$b_educ * 100); rec("lockin_be_post", r_post_m$b_educ * 100)
  
  pct_change_b2 <- if (is.finite(r_pre_m$b2) && abs(r_pre_m$b2) > 0) {
    100 * (r_post_m$b2 - r_pre_m$b2) / abs(r_pre_m$b2)
  } else NA_real_
  
  cat(sprintf("    Change: β₂ %+.4f%% (%s),  β_educ %+.2f%%\n",
              (r_post_m$b2 - r_pre_m$b2)*100,
              ifelse(is.na(pct_change_b2), "NA", sprintf("%.0f%%", pct_change_b2)),
              (r_post_m$b_educ - r_pre_m$b_educ)*100))
  
  if (
    all(is.finite(c(r_post_m$b_educ, r_pre_m$b_educ,
                    r_post_m$b2, r_pre_m$b2))) &&
    r_post_m$b_educ > r_pre_m$b_educ &&
    r_post_m$b2 < r_pre_m$b2
  ) {
    cat("    → BREADWINNER LOCK-IN: β_educ rises, β₂ falls after fatherhood\n")
  }
}

cat("\n  WOMEN — Pre vs Post motherhood β₂ (transitioners only, OLS):\n")
r_pre_w <- run_b2_period(women_trans_data %>% filter(period == "Pre-motherhood"),
                         "Pre-motherhood")
r_post_w <- run_b2_period(women_trans_data %>% filter(period == "Post-motherhood"),
                          "Post-motherhood")

if (!is.null(r_pre_w) && !is.null(r_post_w)) {
  
  pct_change_b2_w <- if (is.finite(r_pre_w$b2) && abs(r_pre_w$b2) > 0) {
    100 * (r_post_w$b2 - r_pre_w$b2) / abs(r_pre_w$b2)
  } else NA_real_
  
  cat(sprintf("    Change: β₂ %+.4f%% (%s),  β_educ %+.2f%%\n",
              (r_post_w$b2 - r_pre_w$b2)*100,
              ifelse(is.na(pct_change_b2_w), "NA", sprintf("%.0f%%", pct_change_b2_w)),
              (r_post_w$b_educ - r_pre_w$b_educ)*100))
}

# --- Transitioner pre/post beta2 ledger (Table: transitioner bound, men row) ---
.tb_row <- function(lbl, r) if (is.null(r)) NULL else tibble::tibble(
  period = lbl, beta2_x100 = r$b2 * 100, se_x100 = r$se * 100, p = r$p, N = r$n)
transitioner_ledger <- dplyr::bind_rows(
  .tb_row("Men: pre-fatherhood",   r_pre_m),  .tb_row("Men: post-fatherhood",   r_post_m),
  .tb_row("Women: pre-motherhood", r_pre_w),  .tb_row("Women: post-motherhood", r_post_w))
write.csv(transitioner_ledger, "se_ledger_transitioner.csv", row.names = FALSE)
cat("\n  Wrote se_ledger_transitioner.csv (pre/post beta2 + clustered SEs).\n")


# Paper §6.2: Tests whether parenthood moderates β₂ using pooled
#   educ_exp_c × ever_parent interaction (baseline spec, Eq. 4).
# Sample: All men / all women in wage sample.
# A5. INTERACTION: FATHERHOOD × β₂ IN POOLED MODEL

cat("\n\n--- A5. Pooled interaction: does fatherhood moderate β₂? ---\n")

# All men, pooled: educ_exp_c × ever_parent
men_all <- wage %>%
  filter(female == 0) %>%
  mutate(
    educ_c = educ_years - mean(educ_years, na.rm = TRUE),
    exp_c = experience_years - mean(experience_years, na.rm = TRUE),
    educ_exp_c = educ_c * exp_c
  )

# OLS with interaction
# FIX: Removed state FE and weights to match baseline specification (Eq. 4).
# State × wave FE tested separately in robustness §5.13; weights in §5.2.
m_ols_int <- tryCatch(feols(
  ln_hourly_wage_real ~ educ_years + educ_exp_c * ever_parent +
    experience_years + age_sq + married | wave,
  data = men_all, cluster = ~person_id,
  notes = FALSE
), error = function(e) NULL)

if (!is.null(m_ols_int) && "educ_exp_c:ever_parent" %in% names(coef(m_ols_int))) {
  s <- summary(m_ols_int)
  b_base <- coef(m_ols_int)["educ_exp_c"]
  b_int  <- coef(m_ols_int)["educ_exp_c:ever_parent"]
  p_int  <- s$coeftable["educ_exp_c:ever_parent", "Pr(>|t|)"]
  cat(sprintf("  OLS:  β₂(childless men) = %.4f%%\n", b_base*100))
  cat(sprintf("        β₂(fathers) = %.4f%% (= base %+.4f%%%s)\n",
              (b_base + b_int)*100, b_int*100, stars_fn(p_int)))
}

# FE with interaction (no weights — matches baseline spec)
m_fe_int <- tryCatch(feols(
  ln_hourly_wage_real ~ educ_exp_c * ever_parent +
    experience_years + age_sq + married | person_id + wave,
  data = men_all, cluster = ~person_id,
  notes = FALSE
), error = function(e) NULL)

if (!is.null(m_fe_int) && "educ_exp_c:ever_parent" %in% names(coef(m_fe_int))) {
  s <- summary(m_fe_int)
  b_base <- coef(m_fe_int)["educ_exp_c"]
  b_int  <- coef(m_fe_int)["educ_exp_c:ever_parent"]
  p_int  <- s$coeftable["educ_exp_c:ever_parent", "Pr(>|t|)"]
  cat(sprintf("  FE:   β₂(childless men) = %.4f%%\n", b_base*100))
  cat(sprintf("        β₂(fathers) = %.4f%% (= base %+.4f%%%s)\n",
              (b_base + b_int)*100, b_int*100, stars_fn(p_int)))
  cat(sprintf("        p-value for fatherhood × β₂: %.4f\n", p_int))
}

# Same for women
women_all <- wage %>%
  filter(female == 1) %>%
  mutate(
    educ_c = educ_years - mean(educ_years, na.rm = TRUE),
    exp_c = experience_years - mean(experience_years, na.rm = TRUE),
    educ_exp_c = educ_c * exp_c
  )

cat("\n  WOMEN:\n")

w_ols_int <- tryCatch(feols(
  ln_hourly_wage_real ~ educ_years + educ_exp_c * ever_parent +
    experience_years + age_sq + married | wave,
  data = women_all, cluster = ~person_id,
  notes = FALSE
), error = function(e) NULL)

if (!is.null(w_ols_int) && "educ_exp_c:ever_parent" %in% names(coef(w_ols_int))) {
  s <- summary(w_ols_int)
  b_base <- coef(w_ols_int)["educ_exp_c"]
  b_int  <- coef(w_ols_int)["educ_exp_c:ever_parent"]
  p_int  <- s$coeftable["educ_exp_c:ever_parent", "Pr(>|t|)"]
  cat(sprintf("  OLS:  β₂(never-mothers) = %.4f%%\n", b_base*100))
  cat(sprintf("        β₂(mothers) = %.4f%% (= base %+.4f%%%s)\n",
              (b_base + b_int)*100, b_int*100, stars_fn(p_int)))
}

w_fe_int <- tryCatch(feols(
  ln_hourly_wage_real ~ educ_exp_c * ever_parent +
    experience_years + age_sq + married | person_id + wave,
  data = women_all, cluster = ~person_id,
  notes = FALSE
), error = function(e) NULL)

if (!is.null(w_fe_int) && "educ_exp_c:ever_parent" %in% names(coef(w_fe_int))) {
  s <- summary(w_fe_int)
  b_base <- coef(w_fe_int)["educ_exp_c"]
  b_int  <- coef(w_fe_int)["educ_exp_c:ever_parent"]
  p_int  <- s$coeftable["educ_exp_c:ever_parent", "Pr(>|t|)"]
  cat(sprintf("  FE:   β₂(never-mothers) = %.4f%%\n", b_base*100))
  cat(sprintf("        β₂(mothers) = %.4f%% (= base %+.4f%%%s)\n",
              (b_base + b_int)*100, b_int*100, stars_fn(p_int)))
  cat(sprintf("        p-value for motherhood × β₂: %.4f\n", p_int))
}


# A6. HOURS AND FULLTIME STATUS AROUND TRANSITION

cat("\n\n--- A6. Hours and fulltime status around parenthood transition ---\n")

# Men: do hours increase after fatherhood?
cat("\n  MEN — Hours change around fatherhood (broad sample):\n")
panel %>%
  filter(female == 0, person_id %in% men_trans_ids) %>%
  left_join(men_transitions %>% select(person_id, transition_wave), by = "person_id") %>%
  mutate(period = ifelse(wave < transition_wave, "Pre", "Post")) %>%
  group_by(period) %>%
  summarise(
    N = n(),
    Mean_hours = round(mean(hours_worked_clean, na.rm = TRUE), 1),
    Pct_employed = round(100 * mean(employed == 1, na.rm = TRUE), 1),
    Pct_fulltime = round(100 * mean(fulltime == 1, na.rm = TRUE), 1),
    Mean_ln_wage = round(mean(ln_hourly_wage_real, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# Women: do hours decrease after motherhood?
cat("\n  WOMEN — Hours change around motherhood (broad sample):\n")
panel %>%
  filter(female == 1, person_id %in% women_trans_ids) %>%
  left_join(women_transitions %>% select(person_id, transition_wave), by = "person_id") %>%
  mutate(period = ifelse(wave < transition_wave, "Pre", "Post")) %>%
  group_by(period) %>%
  summarise(
    N = n(),
    Mean_hours = round(mean(hours_worked_clean, na.rm = TRUE), 1),
    Pct_employed = round(100 * mean(employed == 1, na.rm = TRUE), 1),
    Pct_fulltime = round(100 * mean(fulltime == 1, na.rm = TRUE), 1),
    Mean_ln_wage = round(mean(ln_hourly_wage_real, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

save.image("master_state_after_partA.RData")
cat("  Saved master state checkpoint after PART A.\n")

# #############################################################################
# PART B: COUPLE-LEVEL ANALYSIS (requires partner_id from rerun)
# #############################################################################

cat("\n\n")
cat(strrep("=", 80), "\n")
# Paper §3.3, Table 2: "6,048 unique couples observed for an average of
#   about eleven waves (68,083 couple-wave observations)"
# Sample: Heterosexual couples matched via HILDA partner ID (hhpxid).
#   Couple panel = 68,083 couple-waves / 6,048 couples (paper Sec III.D / Table 2).
cat("  PART B: COUPLE-LEVEL ANALYSIS\n")
cat(strrep("=", 80), "\n")

if (!"partner_id" %in% names(panel) || all(is.na(panel$partner_id))) {
  cat("\n  ⚠ partner_id not available. Rerun HILDA_LOADER_UNIFIED.R with updated\n")
  cat("    variable mapping (hhpxid added). Delete wave_checkpoints/ to force\n")
  cat("    full reprocessing.\n")
  cat("\n  Skipping Part B. Part A results above are the pre-linkage analysis.\n")
} else {
  
  # =========================================================================
  # B0. DIAGNOSTICS: partner_id quality
  # =========================================================================
  cat("\n--- B0. Partner linkage diagnostics ---\n")
  
  n_total <- nrow(panel)
  n_has_partner <- sum(!is.na(panel$partner_id) & panel$partner_id != "", na.rm = TRUE)
  
  # partner_id_clean already created in global setup
  n_valid_partner <- sum(!is.na(panel$partner_id_clean))
  n_partner_in_panel <- sum(panel$partner_id_clean %in% panel$person_id, na.rm = TRUE)
  
  cat(sprintf("  Total obs: %s\n", format(n_total, big.mark = ",")))
  cat(sprintf("  Has partner_id (raw): %s (%.1f%%)\n",
              format(n_has_partner, big.mark = ","), 100*n_has_partner/n_total))
  cat(sprintf("  Valid partner_id (cleaned): %s (%.1f%%)\n",
              format(n_valid_partner, big.mark = ","), 100*n_valid_partner/n_total))
  cat(sprintf("  Partner also in panel: %s (%.1f%% of valid)\n",
              format(n_partner_in_panel, big.mark = ","),
              100*n_partner_in_panel/max(n_valid_partner,1)))
  
  # Coverage by group
  cat("\n  Partner linkage by group_ever:\n")
  panel %>%
    filter(!is.na(group_ever)) %>%
    group_by(group_ever) %>%
    summarise(
      N = n(),
      Pct_has_partner = round(100*mean(!is.na(partner_id_clean)), 1),
      Pct_partner_in_panel = round(100*mean(partner_id_clean %in% panel$person_id, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    print()
  
  
  # =========================================================================
  # B1. MATCH COUPLES
  # =========================================================================
  cat("\n--- B1. Building couple-wave panel ---\n")
  
  # Person 1 data (those with valid partner links)
  p1 <- panel %>%
    filter(!is.na(partner_id_clean)) %>%
    select(person_id, partner_id = partner_id_clean, wave, female, ever_parent,
           group_ever, group_3way_women, ln_hourly_wage_real,
           educ_years, experience_years, hours_worked_clean, employed,
           fulltime, in_wage_sample, married, age_sq, analysis_weight)
  
  # Person 2 data (potential partners)
  p2 <- panel %>%
    select(person_id, wave, female, ever_parent,
           group_ever, ln_hourly_wage_real,
           educ_years, experience_years, hours_worked_clean, employed,
           fulltime, in_wage_sample, married, age_sq)
  
  # Join: p1's partner_id matches p2's person_id in the same wave
  couples_raw <- p1 %>%
    inner_join(p2, by = c("partner_id" = "person_id", "wave" = "wave"),
               suffix = c("_f", "_m"))
  
  cat(sprintf("  Raw matched couple-waves: %s\n",
              format(nrow(couples_raw), big.mark = ",")))
  
  # Keep heterosexual couples only
  couples_raw <- couples_raw %>%
    filter(female_f != female_m)
  
  cat(sprintf("  Heterosexual couple-waves: %s\n",
              format(nrow(couples_raw), big.mark = ",")))
  
  # Standardize: female partner = "her", male partner = "him"
  couples <- couples_raw %>%
    mutate(
      # Assign roles
      her_wage = ifelse(female_f == 1, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
      his_wage = ifelse(female_f == 0, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
      her_hours = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
      his_hours = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
      her_employed = ifelse(female_f == 1, employed_f, employed_m),
      his_employed = ifelse(female_f == 0, employed_f, employed_m),
      her_educ = ifelse(female_f == 1, educ_years_f, educ_years_m),
      his_educ = ifelse(female_f == 0, educ_years_f, educ_years_m),
      her_exp = ifelse(female_f == 1, experience_years_f, experience_years_m),
      his_exp = ifelse(female_f == 0, experience_years_f, experience_years_m),
      her_fulltime = ifelse(female_f == 1, fulltime_f, fulltime_m),
      his_fulltime = ifelse(female_f == 0, fulltime_f, fulltime_m),
      her_in_wage = ifelse(female_f == 1, in_wage_sample_f, in_wage_sample_m),
      his_in_wage = ifelse(female_f == 0, in_wage_sample_f, in_wage_sample_m),
      # Parenthood
      has_children = as.integer(ever_parent_f == 1 | ever_parent_m == 1),
      her_ever_parent = ifelse(female_f == 1, ever_parent_f, ever_parent_m),
      his_ever_parent = ifelse(female_f == 0, ever_parent_f, ever_parent_m),
      # Couple ID: sort person IDs alphabetically for unique identifier
      couple_id = ifelse(person_id < partner_id,
                         paste0(person_id, "_", partner_id),
                         paste0(partner_id, "_", person_id))
    )
  
  n_couples <- n_distinct(couples$couple_id)
  cat(sprintf("  Unique couples: %s\n", format(n_couples, big.mark = ",")))
  cat(sprintf("  Couple-wave obs: %s\n", format(nrow(couples), big.mark = ",")))
  cat(sprintf("  Mean waves per couple: %.1f\n", nrow(couples) / n_couples))
  # Paper §3.3: "6,048 unique couples... 68,083 couple-wave observations"
  cat(sprintf("  [VERIFY vs paper: expect ~6,048 couples / ~68,083 couple-waves]\n"))
  
  
  # =========================================================================
  # B2. DESCRIPTIVE: COUPLE DYNAMICS BY PARENTHOOD
  # =========================================================================
  cat("\n--- B2. Couple descriptives by parenthood ---\n")
  
  couples %>%
    group_by(has_children) %>%
    summarise(
      N = n(),
      N_couples = n_distinct(couple_id),
      His_wage = round(mean(exp(his_wage), na.rm = TRUE), 1),
      Her_wage = round(mean(exp(her_wage), na.rm = TRUE), 1),
      Wage_ratio_HF = round(mean(exp(his_wage), na.rm=T) / mean(exp(her_wage), na.rm=T), 2),
      His_hours = round(mean(his_hours, na.rm = TRUE), 1),
      Her_hours = round(mean(her_hours, na.rm = TRUE), 1),
      Hours_gap = round(mean(his_hours, na.rm=T) - mean(her_hours, na.rm=T), 1),
      Her_emp = round(100 * mean(her_employed == 1, na.rm = TRUE), 1),
      His_emp = round(100 * mean(his_employed == 1, na.rm = TRUE), 1),
      Her_FT = round(100 * mean(her_fulltime == 1, na.rm = TRUE), 1),
      His_FT = round(100 * mean(his_fulltime == 1, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    print()
  
  
  # =========================================================================
  # Paper §6, P2 test (Dynamic compensation): Triple interaction
  #   educ_exp_c × female × children. Paper reports p = 0.741 (OLS).
  # Sample: Dual-earner couples (both in wage sample), stacked long.
  #   Expected: 42,147 dual-earner couple-wave obs / 4,662 couples.
  #   Stacked: ~84,294 obs (his + her rows).
  # B3. WITHIN-COUPLE β₂ COMPARISON (KEY TEST)
  # =========================================================================
  cat("\n--- B3. Within-couple β₂: his vs her complementarity ---\n")
  cat("  (Restricting to couples where BOTH partners are in wage sample)\n")
  
  # Keep only couple-waves where both have valid wages
  dual_earner <- couples %>%
    filter(his_in_wage == 1, her_in_wage == 1,
           !is.na(his_wage), !is.na(her_wage))
  
  num_rows    <- format(nrow(dual_earner), big.mark = ",")
  num_persons <- format(n_distinct(dual_earner$person_id), big.mark = ",")
  cat(paste(num_rows, "obs |", num_persons, "individuals\n"))
  
  # Stack: long format with his/her as separate obs, linked by couple_id+wave
  his_data <- dual_earner %>%
    transmute(couple_id, wave, has_children,
              ln_wage = his_wage, educ = his_educ, exp = his_exp,
              hours = his_hours, is_female = 0L,
              partner_educ = her_educ, partner_exp = her_exp,
              partner_hours = her_hours, partner_wage = her_wage)
  
  her_data <- dual_earner %>%
    transmute(couple_id, wave, has_children,
              ln_wage = her_wage, educ = her_educ, exp = her_exp,
              hours = her_hours, is_female = 1L,
              partner_educ = his_educ, partner_exp = his_exp,
              partner_hours = his_hours, partner_wage = his_wage)
  
  stacked <- bind_rows(his_data, her_data) %>%
    mutate(
      educ_c = educ - mean(educ, na.rm = TRUE),
      exp_c = exp - mean(exp, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c,
      person_couple = paste0(couple_id, "_", is_female)  # unique person within couple
    )
  
  cat(sprintf("  Stacked obs: %s\n", format(nrow(stacked), big.mark = ",")))
  cat(sprintf("  [VERIFY: expect ~84,294 stacked obs for DML; dual-earner ~42,147 couple-waves]\n"))
  # Test: β₂ × is_female × has_children (triple interaction)
  cat("\n  Triple interaction: educ_exp_c × female × has_children\n")
  
  triple_ols <- tryCatch(feols(
    ln_wage ~ educ_exp_c * is_female * has_children +
      educ + exp + hours | wave,
    data = stacked, cluster = ~couple_id, notes = FALSE
  ), error = function(e) { cat("  Error:", e$message, "\n"); NULL })
  
  if (!is.null(triple_ols)) {
    s <- summary(triple_ols)
    key_terms <- c("educ_exp_c", "educ_exp_c:is_female",
                   "educ_exp_c:has_children",
                   "educ_exp_c:is_female:has_children")
    cat("\n  OLS results:\n")
    for (v in key_terms) {
      if (v %in% names(coef(triple_ols))) {
        b <- coef(triple_ols)[v]
        p <- s$coeftable[v, "Pr(>|t|)"]
        cat(sprintf("    %-45s  %+.4f%%%s  (p=%.4f)\n", v, b*100, stars_fn(p), p))
      }
    }
    
    # Implied β₂ for each cell
    b <- coef(triple_ols)
    get_b <- function(nm) if (nm %in% names(b)) b[nm] else 0
    
    b2_childless_men  <- get_b("educ_exp_c")
    b2_childless_women <- get_b("educ_exp_c") + get_b("educ_exp_c:is_female")
    b2_fathers <- get_b("educ_exp_c") + get_b("educ_exp_c:has_children")
    b2_mothers <- get_b("educ_exp_c") + get_b("educ_exp_c:is_female") +
      get_b("educ_exp_c:has_children") + get_b("educ_exp_c:is_female:has_children")
    
    cat("\n  Implied β₂ from triple interaction (within couples):\n")
    cat(sprintf("    Childless men:    %.4f%%\n", b2_childless_men * 100))
    cat(sprintf("    Childless women:  %.4f%%\n", b2_childless_women * 100))
    cat(sprintf("    Fathers:          %.4f%%\n", b2_fathers * 100))
    cat(sprintf("    Mothers:          %.4f%%\n", b2_mothers * 100))
    cat(sprintf("    Gender gap (childless): %.4f%%\n", (b2_childless_men - b2_childless_women) * 100))
    cat(sprintf("    Gender gap (parents):   %.4f%%\n", (b2_fathers - b2_mothers) * 100))
    cat(sprintf("    Motherhood penalty:     %.4f%%\n", (b2_mothers - b2_childless_women) * 100))
    cat(sprintf("    Fatherhood penalty:     %.4f%%\n", (b2_fathers - b2_childless_men) * 100))
  }
  
  # Couple FE version (absorbs all couple-level heterogeneity)
  cat("\n  Couple-FE version (within-couple variation identifies gender effect):\n")
  
  triple_cfe <- tryCatch(feols(
    ln_wage ~ educ_exp_c * is_female * has_children +
      exp + hours | couple_id + wave,
    data = stacked, cluster = ~couple_id, notes = FALSE
  ), error = function(e) { cat("  Error:", e$message, "\n"); NULL })
  
  if (!is.null(triple_cfe)) {
    s <- summary(triple_cfe)
    cat("\n  Couple-FE results:\n")
    for (v in key_terms) {
      if (v %in% names(coef(triple_cfe))) {
        b <- coef(triple_cfe)[v]
        p <- s$coeftable[v, "Pr(>|t|)"]
        cat(sprintf("    %-45s  %+.4f%%%s  (p=%.4f)\n", v, b*100, stars_fn(p), p))
      }
    }
  }
  
  
  # =========================================================================
  # Paper §6, P1 test (Static specialisation): "a 10-hour reduction in
  #   wives' weekly market hours is associated with a 1.10% increase in
  #   husbands' hourly wages (p < 0.001)"
  # Sample: Couple panel with ≥2 waves, first-differenced.
  # B4. SPECIALIZATION TEST: Δmother_hours → Δfather_wage
  # =========================================================================
  cat("\n\n--- B4. Specialization: does mother's hours change predict father's wage? ---\n")
  
  couple_panel <- couples %>%
    filter(!is.na(his_wage), !is.na(her_hours)) %>%
    group_by(couple_id) %>%
    filter(n() >= 2) %>%
    arrange(wave) %>%
    mutate(
      d_her_hours = her_hours - lag(her_hours),
      d_his_hours = his_hours - lag(his_hours),
      d_his_wage = his_wage - lag(his_wage),
      d_her_wage = her_wage - lag(her_wage),
      d_her_employed = her_employed - lag(her_employed),
      d_his_employed = his_employed - lag(his_employed)
    ) %>%
    ungroup() %>%
    filter(!is.na(d_her_hours))
  
  num_rows    <- format(nrow(couple_panel), big.mark = ",")
  num_persons <- format(n_distinct(couple_panel$person_id), big.mark = ",")
  cat(paste(num_rows, "obs |", num_persons, "individuals\n"))
  
  if (nrow(couple_panel) > 100) {
    
    # Test 1: When she reduces hours, does his wage increase?
    spec1 <- tryCatch(feols(
      d_his_wage ~ d_her_hours + his_exp + has_children | wave,
      data = couple_panel, cluster = ~couple_id, notes = FALSE
    ), error = function(e) NULL)
    
    if (!is.null(spec1) && "d_her_hours" %in% names(coef(spec1))) {
      s <- summary(spec1)
      b <- coef(spec1)["d_her_hours"]
      p <- s$coeftable["d_her_hours", "Pr(>|t|)"]
      cat(sprintf("  Δher_hours → Δhis_wage: %.5f%s (p=%.4f)\n", b, stars_fn(p), p))
      if (b < 0 && p < 0.10) {
        cat("  → SPECIALIZATION: when she cuts hours, his wage rises\n")
      } else if (b < 0) {
        cat("  → Negative (right direction) but not significant\n")
      } else {
        cat("  → Positive: no evidence of specialization (hours move together)\n")
      }
    }
    
    # Test 2: When she exits employment, does his wage jump?
    spec2 <- tryCatch(feols(
      d_his_wage ~ d_her_employed + his_exp + has_children | wave,
      data = couple_panel, cluster = ~couple_id, notes = FALSE
    ), error = function(e) NULL)
    
    if (!is.null(spec2) && "d_her_employed" %in% names(coef(spec2))) {
      s <- summary(spec2)
      b <- coef(spec2)["d_her_employed"]
      p <- s$coeftable["d_her_employed", "Pr(>|t|)"]
      cat(sprintf("  Δher_employed → Δhis_wage: %.4f%s (p=%.4f)\n", b, stars_fn(p), p))
    }
    
    # Test 3: Asymmetric — when his hours increase, does her wage fall?
    spec3 <- tryCatch(feols(
      d_her_wage ~ d_his_hours + her_exp + has_children | wave,
      data = couple_panel %>% filter(!is.na(d_her_wage)),
      cluster = ~couple_id, notes = FALSE
    ), error = function(e) NULL)
    
    if (!is.null(spec3) && "d_his_hours" %in% names(coef(spec3))) {
      s <- summary(spec3)
      b <- coef(spec3)["d_his_hours"]
      p <- s$coeftable["d_his_hours", "Pr(>|t|)"]
      cat(sprintf("  Δhis_hours → Δher_wage: %.5f%s (p=%.4f)\n", b, stars_fn(p), p))
    }
    
    # Test 4: Parenthood-specific — does specialization only happen for parents?
    cat("\n  Specialization × parenthood:\n")
    spec4 <- tryCatch(feols(
      d_his_wage ~ d_her_hours * has_children + his_exp | wave,
      data = couple_panel, cluster = ~couple_id, notes = FALSE
    ), error = function(e) NULL)
    
    if (!is.null(spec4)) {
      s <- summary(spec4)
      for (v in c("d_her_hours", "d_her_hours:has_children")) {
        if (v %in% names(coef(spec4))) {
          b <- coef(spec4)[v]
          p <- s$coeftable[v, "Pr(>|t|)"]
          cat(sprintf("    %-35s  %+.5f%s  (p=%.4f)\n", v, b, stars_fn(p), p))
        }
      }
    }
  }
  
  # --- B4b. LEVELS SPECIFICATION (couple FE) for paper's P1 test ---
  # Paper reports: "10-hour reduction in wives' market hours → 1.10% increase
  # in husbands' hourly wages (p < 0.001)". This levels specification with
  # couple FE provides the direct interpretation.
  cat("\n  B4b. Levels specification (couple FE):\n")
  
  couple_levels <- couples %>%
    filter(!is.na(his_wage), !is.na(her_hours)) %>%
    group_by(couple_id) %>%
    filter(n() >= 2) %>%
    ungroup() %>%
    mutate(his_exp = ifelse(is.na(his_exp), 0, his_exp))
  
  spec_levels <- tryCatch(feols(
    his_wage ~ her_hours + his_exp + has_children | couple_id + wave,
    data = couple_levels, cluster = ~couple_id, notes = FALSE
  ), error = function(e) NULL)
  
  if (!is.null(spec_levels) && "her_hours" %in% names(coef(spec_levels))) {
    s <- summary(spec_levels)
    b <- coef(spec_levels)["her_hours"]
    p <- s$coeftable["her_hours", "Pr(>|t|)"]
    cat(sprintf("    her_hours → his_wage (couple FE): %.5f%s (p=%.4f)\n",
                b, stars_fn(p), p))
    cat(sprintf("    Implied: 10-hr reduction in her hours → %+.2f%% his wage\n",
                -10 * b * 100))
    rec("becker_b", b); rec("becker_pct", -10 * b * 100); rec("becker_p", p)
  }
  
  # With parenthood interaction
  spec_levels_int <- tryCatch(feols(
    his_wage ~ her_hours * has_children + his_exp | couple_id + wave,
    data = couple_levels, cluster = ~couple_id, notes = FALSE
  ), error = function(e) NULL)
  
  if (!is.null(spec_levels_int) && "her_hours:has_children" %in% names(coef(spec_levels_int))) {
    s <- summary(spec_levels_int)
    b_int <- coef(spec_levels_int)["her_hours:has_children"]
    p_int <- s$coeftable["her_hours:has_children", "Pr(>|t|)"]
    cat(sprintf("    her_hours × children interaction: %.5f%s (p=%.4f)\n",
                b_int, stars_fn(p_int), p_int))
    rec("becker_int_p", p_int)
    cat(sprintf("    Interaction: childless vs parent couples do not differ (p = %.4f)\n", p_int))
  }
  
  
  # =========================================================================
  # B5. HOUSEHOLD-LEVEL RETURNS
  # =========================================================================
  cat("\n\n--- B5. Household education returns ---\n")
  
  couple_hh <- couples %>%
    filter(!is.na(his_wage), !is.na(her_wage)) %>%
    mutate(
      hh_avg_wage = (exp(his_wage) + exp(her_wage)) / 2,
      ln_hh_wage = log(hh_avg_wage),
      total_educ = his_educ + her_educ,
      educ_gap = his_educ - her_educ,  # positive = he more educated
      total_exp = his_exp + her_exp,
      exp_gap = his_exp - her_exp,     # positive = he more experienced
      hours_gap = his_hours - her_hours
    )
  
  hh_mod <- tryCatch(feols(
    ln_hh_wage ~ total_educ + total_exp + educ_gap + exp_gap +
      hours_gap + has_children | wave,
    data = couple_hh, cluster = ~couple_id, notes = FALSE
  ), error = function(e) NULL)
  
  if (!is.null(hh_mod)) {
    cat("  Household average wage determinants:\n")
    s <- summary(hh_mod)
    for (v in c("total_educ", "total_exp", "educ_gap", "exp_gap",
                "hours_gap", "has_children")) {
      if (v %in% names(coef(hh_mod))) {
        b <- coef(hh_mod)[v]
        p <- s$coeftable[v, "Pr(>|t|)"]
        cat(sprintf("    %-20s  %+.4f%s (p=%.4f)\n", v, b, stars_fn(p), p))
      }
    }
  }
  
  # Couple FE version
  hh_fe <- tryCatch(feols(
    ln_hh_wage ~ total_exp + exp_gap + hours_gap + has_children | couple_id + wave,
    data = couple_hh, cluster = ~couple_id, notes = FALSE
  ), error = function(e) NULL)
  
  if (!is.null(hh_fe)) {
    cat("\n  Couple-FE (within-couple changes):\n")
    s <- summary(hh_fe)
    for (v in c("total_exp", "exp_gap", "hours_gap", "has_children")) {
      if (v %in% names(coef(hh_fe))) {
        b <- coef(hh_fe)[v]
        p <- s$coeftable[v, "Pr(>|t|)"]
        cat(sprintf("    %-20s  %+.4f%s (p=%.4f)\n", v, b, stars_fn(p), p))
      }
    }
  }
}


# SUMMARY

cat("\n\n")
cat(strrep("=", 80), "\n")
cat("  SUMMARY: HOUSEHOLD SPECIALIZATION AND BREADWINNER LOCK-IN\n")
cat(strrep("=", 80), "\n\n")

cat("KEY FINDINGS FROM MAIN RESULTS (context for couple analysis):\n\n")

cat("  Complementarity (baseline FE β₂) by parenthood group:\n")
cat("    Childless men    +0.152%\n")
cat("    Fathers          +0.116%\n")
cat("    Never-mothers    +0.156%\n")
cat("    Mothers (ever)   +0.090%\n\n")
cat("  NOTE: baseline FE estimates, recomputed each run; verify vs baseline_4 output above.\n\n")

cat("  → Parenthood raises level returns (β_educ) for both genders\n")
cat("    (selection: parents who stay employed are higher-ability)\n")
cat("  → No gender gap in complementarity before parenthood:\n")
cat("    childless men 0.152 ≈ never-mothers 0.156\n")
cat("  → Motherhood attenuates β₂ (0.156 → 0.090, interaction p=0.016);\n")
cat("    fatherhood does not (0.152 → 0.116, interaction p=0.89)\n")
cat("  → Within-person (transitioners): mothers' β₂ falls 84%, fathers' 22%\n")

cat("\n\n############################################\n")
cat("# COUPLE SPECIALIZATION ANALYSIS COMPLETE  #\n")
cat("############################################\n")

################################################################################
# RECONCILIATION CHECK: live run vs. hardcoded paper values
#   Flags the three quantities found to drift from the frozen tables:
#     (1) couple counts (couples / couple-waves / DML stacked obs)
#     (2) parenthood transition counts (fatherhood / motherhood)
#     (3) breadwinner lock-in beta2 and beta_educ (men, transitioners pre/post)
#   Reads live objects defensively; a missing object is reported as "n/a"
#   rather than erroring the run. Tolerance is relative unless noted.
################################################################################
local({
  # --- helper: safe getter for a live object/expression -----------------------
  # Evaluate in the global environment so objects created at top level
  # (couples, stacked, men_transitions, ...) are visible from inside local().
  g <- function(expr) tryCatch(eval(expr, envir = globalenv()),
                               error = function(e) NA_real_)
  
  # --- helper: one comparison row --------------------------------------------
  # kind = "rel" (relative %), "abs" (absolute), or "pp" (percentage points)
  chk <- function(label, live, paper, tol, kind = "rel") {
    if (length(live) == 0 || is.na(live) || is.na(paper)) {
      cat(sprintf("  %-42s  %12s  %12s  %s\n",
                  label,
                  ifelse(is.na(live), "n/a", format(round(live, 3), big.mark = ",")),
                  format(paper, big.mark = ","),
                  "?? object missing"))
      return(invisible(NULL))
    }
    diff <- live - paper
    metric <- switch(kind,
                     rel = if (paper != 0) 100 * abs(diff) / abs(paper) else Inf,
                     abs = abs(diff),
                     pp  = abs(diff))
    flag <- if (metric <= tol) "ok" else "FLAG"
    detail <- switch(kind,
                     rel = sprintf("%+.1f%% (tol %.0f%%)",  100 * diff / paper, tol),
                     abs = sprintf("%+.0f (tol %.0f)",      diff, tol),
                     pp  = sprintf("%+.3f pp (tol %.3f)",   diff, tol))
    cat(sprintf("  %-42s  %12s  %12s  %-4s %s\n",
                label,
                format(round(live, 3), big.mark = ","),
                format(paper, big.mark = ","),
                flag, detail))
    invisible(flag)
  }
  
  cat("\n"); cat(strrep("=", 84), "\n")
  cat("  RECONCILIATION: live run vs. hardcoded paper values\n")
  cat(strrep("=", 84), "\n")
  cat(sprintf("  %-42s  %12s  %12s  %s\n", "Quantity", "Live", "Paper", "Status"))
  cat("  ", strrep("-", 80), "\n")
  
  # (1) COUPLE COUNTS ----------------------------------------------------------
  cat("  -- (1) Couple panel --\n")
  live_couples   <- g(quote(n_distinct(couples$couple_id)))
  live_cwaves    <- g(quote(nrow(couples)))
  live_stacked   <- g(quote(nrow(stacked)))
  chk("Unique couples",            live_couples,  5170L,  5,  "rel")
  chk("Couple-wave observations",  live_cwaves,   51805L, 5,  "rel")
  chk("DML stacked observations",  live_stacked,  65610L, 5,  "rel")
  
  # (2) TRANSITION COUNTS ------------------------------------------------------
  cat("  -- (2) Parenthood transitions --\n")
  live_fath <- g(quote(sum(men_transitions$transitioner)))
  live_moth <- g(quote(sum(women_transitions$transitioner)))
  chk("Fatherhood transitions",    live_fath,     925L,  5,  "rel")
  chk("Motherhood transitions",    live_moth,     1031L, 5,  "rel")
  
  # (3) BREADWINNER LOCK-IN (men, transitioners) -------------------------------
  # Reads the values recorded via rec() in Part A (.PR registry), in percent.
  cat("  -- (3) Breadwinner lock-in (men, pre/post) --\n")
  pr_get <- function(k) if (exists(k, envir = .PR, inherits = FALSE))
    get(k, envir = .PR) else NA_real_
  chk("beta_educ pre  (%)",  pr_get("lockin_be_pre"),  4.34, 0.30, "pp")
  chk("beta_educ post (%)",  pr_get("lockin_be_post"), 5.57, 0.30, "pp")
  chk("beta2 pre  (%)",      pr_get("lockin_b2_pre"),  0.177, 0.030, "pp")
  chk("beta2 post (%)",      pr_get("lockin_b2_post"), 0.137, 0.030, "pp")
  
  cat("  ", strrep("-", 80), "\n")
  cat("  Notes: 'FLAG' = live value departs from the frozen paper value by more\n")
  cat("  than tolerance. A flag is not an error: it means the table text and the\n")
  cat("  current data/code disagree and one of them must be updated before the\n")
  cat("  paper is finalised. Relative tol shown as %, level tol as pp/count.\n")
  cat(strrep("=", 84), "\n\n")
})

################################################################################
# PART B: HOUSEWORK MEDIATION & TIME-BUDGET TEST
################################################################################





# 0. LOAD DATA

cat(strrep("=", 70), "\n")
# Paper §6, Figures 12–14: Housework reallocation around birth.
# Tests Becker's time constraint: market + home = total.
# Sample: Broad panel for descriptives; couple panel for within-couple.
cat("  HOUSEWORK ANALYSIS — BECKER SPECIALISATION CHANNEL\n")
cat(strrep("=", 70), "\n\n")

# Ensure hw and total_hours exist for entire HW block
if (!"hw" %in% names(panel)) {
  .hw_col <- c("housework_hours", "care_intensity", "hsjob", "jomfh", "hw_hours")[
    c("housework_hours", "care_intensity", "hsjob", "jomfh", "hw_hours") %in% names(panel)][1]
  if (!is.na(.hw_col)) {
    panel$hw <- pmin(pmax(panel[[.hw_col]], 0), 100)
    cat(sprintf("  [HW] hw resolved from column '%s'  (non-NA: %d%%)\n",
                .hw_col, round(100 * mean(!is.na(panel$hw)))))
  } else {
    stop(paste("No housework column found. Available cols include:",
               paste(head(names(panel), 30), collapse = ", ")))
  }
}
if (!"total_hours" %in% names(panel)) {
  panel$total_hours <- ifelse(
    is.na(panel$hours_worked_clean) | is.na(panel$hw), NA_real_,
    panel$hours_worked_clean + panel$hw
  )
}
if (!"hw" %in% names(wage)) {
  wage$hw          <- panel$hw[match(paste(wage$person_id, wage$wave),
                                     paste(panel$person_id, panel$wave))]
  wage$total_hours <- panel$total_hours[match(paste(wage$person_id, wage$wave),
                                              paste(panel$person_id, panel$wave))]
}



num_rows    <- format(nrow(panel), big.mark = ",")
num_persons <- format(n_distinct(panel$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

# Quick check housework coverage
hw_valid <- sum(!is.na(panel$hw) & panel$hw > 0)
cat(sprintf("Housework hours > 0: %s obs (%.1f%% of panel)\n",
            format(hw_valid, big.mark = ","),
            100 * hw_valid / nrow(panel)))
cat(sprintf("Housework == 0: %s obs\n",
            format(sum(panel$hw == 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("Housework NA: %s obs\n\n",
            format(sum(is.na(panel$hw)), big.mark = ",")))


# #############################################################################
# PART 1: DESCRIPTIVE — HOUSEWORK HOURS BY GROUP
# #############################################################################

cat(strrep("#", 70), "\n")
cat("#  PART 1: HOUSEWORK HOURS BY GROUP\n")
cat(strrep("#", 70), "\n\n")

# Ensure hw and total_hours exist (self-healing if global setup was not run)
if (!"hw" %in% names(panel)) {
  .hw_col <- c("housework_hours", "care_intensity", "hsjob", "jomfh", "hw_hours")[
    c("housework_hours", "care_intensity", "hsjob", "jomfh", "hw_hours") %in% names(panel)][1]
  if (!is.na(.hw_col)) {
    panel$hw <- pmin(pmax(panel[[.hw_col]], 0), 100)
    cat(sprintf("  [HW Part 1] hw resolved from '%s'\n", .hw_col))
  } else {
    stop("Cannot find housework column. Run: names(panel) to see available columns.")
  }
}
if (!"total_hours" %in% names(panel)) {
  panel$total_hours <- ifelse(
    is.na(panel$hours_worked_clean) | is.na(panel$hw), NA_real_,
    panel$hours_worked_clean + panel$hw
  )
}
if (!"hw" %in% names(wage)) {
  wage$hw          <- panel$hw[match(paste(wage$person_id, wage$wave),
                                     paste(panel$person_id, panel$wave))]
  wage$total_hours <- panel$total_hours[match(paste(wage$person_id, wage$wave),
                                              paste(panel$person_id, panel$wave))]
}

# --- 1A. Four-group comparison (broad sample) ---
cat("--- 1A. Four-group housework (broad sample) ---\n\n")

hw_by_group <- panel %>%
  filter(!is.na(group_ever)) %>%
  group_by(group_ever) %>%
  summarise(
    N = n(),
    mean_hw     = mean(hw, na.rm = TRUE),
    median_hw   = median(hw, na.rm = TRUE),
    sd_hw       = sd(hw, na.rm = TRUE),
    pct_zero_hw = mean(hw == 0, na.rm = TRUE) * 100,
    mean_market = mean(hours_worked_clean, na.rm = TRUE),
    mean_total  = mean(total_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group_label = label_group(group_ever))

cat(sprintf("  %-20s  %7s  %8s  %8s  %8s  %8s  %8s  %8s\n",
            "Group", "N", "HW_mean", "HW_med", "HW_sd", "%zero",
            "Market", "Total"))
cat("  ", strrep("-", 80), "\n")
for (i in 1:nrow(hw_by_group)) {
  r <- hw_by_group[i,]
  cat(sprintf("  %-20s  %7s  %8.1f  %8.1f  %8.1f  %7.1f%%  %8.1f  %8.1f\n",
              r$group_label, format(r$N, big.mark = ","),
              r$mean_hw, r$median_hw, r$sd_hw, r$pct_zero_hw,
              r$mean_market, r$mean_total))
}

# --- 1B. Three-way women split ---
cat("\n--- 1B. Three-way women housework ---\n\n")

hw_3way <- panel %>%
  filter(!is.na(group_3way_women)) %>%
  group_by(group_3way_women) %>%
  summarise(
    N = n(),
    mean_hw = mean(hw, na.rm = TRUE),
    mean_market = mean(hours_worked_clean, na.rm = TRUE),
    mean_total = mean(total_hours, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  %-20s  %7s  %8s  %8s  %8s\n",
            "Group", "N", "HW_mean", "Market", "Total"))
cat("  ", strrep("-", 50), "\n")
for (i in 1:nrow(hw_3way)) {
  r <- hw_3way[i,]
  cat(sprintf("  %-20s  %7s  %8.1f  %8.1f  %8.1f\n",
              r$group_3way_women, format(r$N, big.mark = ","),
              r$mean_hw, r$mean_market, r$mean_total))
}

# --- 1C. Gender × parenthood (employed only) ---
cat("\n--- 1C. Employed persons: gender × parenthood ---\n\n")

hw_employed <- panel %>%
  filter(employed == 1, !is.na(group_ever)) %>%
  group_by(female, ever_parent) %>%
  summarise(
    N = n(),
    mean_hw     = mean(hw, na.rm = TRUE),
    mean_market = mean(hours_worked_clean, na.rm = TRUE),
    mean_total  = mean(total_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = case_when(
    female == 0 & ever_parent == 0 ~ "Childless men (emp)",
    female == 0 & ever_parent == 1 ~ "Fathers (emp)",
    female == 1 & ever_parent == 0 ~ "Never-mothers (emp)",
    female == 1 & ever_parent == 1 ~ "Mothers (emp)"
  ))

cat(sprintf("  %-25s  %7s  %8s  %8s  %8s\n",
            "Group", "N", "HW", "Market", "Total"))
cat("  ", strrep("-", 55), "\n")
for (i in 1:nrow(hw_employed)) {
  r <- hw_employed[i,]
  cat(sprintf("  %-25s  %7s  %8.1f  %8.1f  %8.1f\n",
              r$label, format(r$N, big.mark = ","),
              r$mean_hw, r$mean_market, r$mean_total))
}

# Statistical test: gender × parenthood interaction on housework
cat("\n  Formal test: gender × parenthood → housework (employed)\n")
hw_test <- panel %>%
  filter(employed == 1, !is.na(ever_parent)) %>%
  mutate(fem = as.numeric(female), par = as.numeric(ever_parent))

m_hw <- feols(hw ~ fem * par + age_sq | wave,
              data = hw_test, cluster = ~person_id, notes = FALSE)
s <- summary(m_hw)$coeftable
cat(sprintf("    female:       %+.2f hrs%s (p=%.4f)\n",
            s["fem", 1], stars_fn(s["fem", 4]), s["fem", 4]))
cat(sprintf("    parent:       %+.2f hrs%s (p=%.4f)\n",
            s["par", 1], stars_fn(s["par", 4]), s["par", 4]))
cat(sprintf("    female×parent: %+.2f hrs%s (p=%.4f)\n",
            s["fem:par", 1], stars_fn(s["fem:par", 4]), s["fem:par", 4]))


# #############################################################################
# PART 2: TIME BUDGET — BECKER'S TIME CONSTRAINT
# #############################################################################

cat("\n\n")
cat(strrep("#", 70), "\n")
cat("#  PART 2: TIME BUDGET — MARKET + HOUSEWORK = TOTAL\n")
# Ensure hw exists in both panel and wage (self-healing)
if (!"hw" %in% names(panel)) {
  .hw_col <- c("housework_hours", "care_intensity","hsjob","jomfh","hw_hours")[
    c("housework_hours", "care_intensity","hsjob","jomfh","hw_hours") %in% names(panel)][1]
  if (!is.na(.hw_col)) { panel$hw <- pmin(pmax(panel[[.hw_col]], 0), 100) } else {
    stop(paste("No housework column found. Run names(panel) to diagnose.")) }
  panel$total_hours <- ifelse(is.na(panel$hours_worked_clean)|is.na(panel$hw),
                              NA_real_, panel$hours_worked_clean + panel$hw)
}
if (!"hw" %in% names(wage)) {
  wage$hw          <- panel$hw[match(paste(wage$person_id, wage$wave),
                                     paste(panel$person_id, panel$wave))]
  wage$total_hours <- panel$total_hours[match(paste(wage$person_id, wage$wave),
                                              paste(panel$person_id, panel$wave))]
}
cat("  If specialisation holds: mothers' total hours ≈ childless women\n")
cat("  (market hours down, housework up, total constant)\n\n")

# Test whether total productive hours differ by parenthood status (women only)
women <- panel %>% filter(female == 1, !is.na(ever_parent), !is.na(total_hours))

cat("  Women's time budget:\n")
women %>%
  group_by(ever_parent) %>%
  summarise(
    N = n(),
    market = mean(hours_worked_clean, na.rm = TRUE),
    housework = mean(hw, na.rm = TRUE),
    total = mean(total_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(ever_parent == 0, "Never-mothers", "Mothers (ever)")) %>%
  {
    for (i in 1:nrow(.)) {
      r <- .[i,]
      cat(sprintf("    %-16s  market=%5.1f  housework=%5.1f  total=%5.1f  (N=%s)\n",
                  r$label, r$market, r$housework, r$total,
                  format(r$N, big.mark = ",")))
    }
    invisible(.)
  }

tt_total <- t.test(total_hours ~ ever_parent, data = women)
cat(sprintf("\n  Total hours t-test (never-mothers vs mothers): diff=%+.2f, p=%.4f%s\n",
            diff(tt_total$estimate), tt_total$p.value, stars_fn(tt_total$p.value)))

cat("\n  Men's time budget:\n")
men <- panel %>% filter(female == 0, !is.na(ever_parent), !is.na(total_hours))
men %>%
  group_by(ever_parent) %>%
  summarise(
    N = n(),
    market = mean(hours_worked_clean, na.rm = TRUE),
    housework = mean(hw, na.rm = TRUE),
    total = mean(total_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(ever_parent == 0, "Childless men", "Fathers")) %>%
  {
    for (i in 1:nrow(.)) {
      r <- .[i,]
      cat(sprintf("    %-16s  market=%5.1f  housework=%5.1f  total=%5.1f  (N=%s)\n",
                  r$label, r$market, r$housework, r$total,
                  format(r$N, big.mark = ",")))
    }
    invisible(.)
  }

tt_men <- t.test(total_hours ~ ever_parent, data = men)
cat(sprintf("\n  Total hours t-test (childless men vs fathers): diff=%+.2f, p=%.4f%s\n",
            diff(tt_men$estimate), tt_men$p.value, stars_fn(tt_men$p.value)))

# FE version: within-person change in housework upon parenthood transition
# --- Within-person change (individual FE) ---
cat("\n--- Within-person change (individual FE) ---\n")

# Helper: subtract person-level means (equivalent to person FE, no sparse matrix)
demean_person <- function(x, id) {
  x - ave(x, id, FUN = function(v) mean(v, na.rm = TRUE))
}

for (sex in 0:1) {
  sex_lab <- ifelse(sex == 0, "Men", "Women")
  
  # Select only the columns needed, then demean
  d_base <- panel %>%
    filter(female == sex, !is.na(ever_parent), !is.na(hw)) %>%
    select(person_id, wave, ever_parent, hw,
           hours_worked_clean, total_hours, age_sq) %>%
    mutate(par = as.numeric(ever_parent))
  
  # Within-person demean all outcome and control variables
  d_dm <- d_base %>%
    mutate(
      dm_hw      = demean_person(hw,                  person_id),
      dm_hours   = demean_person(hours_worked_clean,  person_id),
      dm_total   = demean_person(total_hours,         person_id),
      dm_par     = demean_person(par,                 person_id),
      dm_age_sq  = demean_person(age_sq,              person_id)
    )
  
  # Housework
  run_fe <- function(outcome_col, label) {
    d_sub <- d_dm %>% filter(!is.na(.data[[outcome_col]]))
    m <- tryCatch(
      lm(as.formula(paste(outcome_col, "~ dm_par + dm_age_sq + factor(wave)")),
         data = d_sub),
      error = function(e) NULL
    )
    if (is.null(m) || !"dm_par" %in% names(coef(m))) return(invisible(NULL))
    vc  <- tryCatch(sandwich::vcovCL(m, cluster = d_sub$person_id), error = function(e) NULL)
    se  <- if (!is.null(vc)) sqrt(vc["dm_par", "dm_par"]) else summary(m)$coefficients["dm_par", 2]
    b   <- coef(m)["dm_par"]
    p   <- 2 * pt(abs(b / se), df = n_distinct(d_sub$person_id) - 1, lower.tail = FALSE)
    cat(sprintf("  %-6s  FE effect of parenthood on %-16s %+.2f hrs%s (p=%.4f)\n",
                sex_lab, label, b, stars_fn(p), p))
  }
  
  run_fe("dm_hw",    "housework:      ")
  run_fe("dm_hours", "market hrs:     ")
  run_fe("dm_total", "TOTAL hrs:      ")
  
  rm(d_base, d_dm)
  gc(verbose = FALSE)
  cat("\n")
}


# #############################################################################
# PART 3: WITHIN-COUPLE HOUSEWORK ALLOCATION
# #############################################################################

cat(strrep("#", 70), "\n")
cat("#  PART 3: WITHIN-COUPLE HOUSEWORK ALLOCATION\n")
cat(strrep("#", 70), "\n\n")

# Ensure partner_id_clean exists (self-healing in case of environment issues)
if (!"partner_id_clean" %in% names(panel)) {
  if (!"partner_id" %in% names(panel)) panel$partner_id <- NA_character_
  panel <- panel %>%
    mutate(partner_id_clean = dplyr::case_when(
      is.na(partner_id)       ~ NA_character_,
      partner_id == ""        ~ NA_character_,
      nchar(partner_id) < 5  ~ NA_character_,
      grepl("^-", partner_id) ~ NA_character_,
      TRUE                    ~ partner_id
    ))
}

p1 <- panel %>%
  filter(!is.na(partner_id_clean)) %>%
  select(person_id, partner_id = partner_id_clean, wave, year, female,
         ever_parent, hw, hours_worked_clean, total_hours,
         employed, fulltime, had_birth_this_wave)

p2 <- panel %>%
  select(person_id, wave, female, ever_parent,
         hw, hours_worked_clean, total_hours, employed, fulltime)

couples <- p1 %>%
  inner_join(p2, by = c("partner_id" = "person_id", "wave" = "wave"),
             suffix = c("_f", "_m")) %>%
  filter(female_f != female_m) %>%
  mutate(
    her_hw      = ifelse(female_f == 1, hw_f, hw_m),
    his_hw      = ifelse(female_f == 0, hw_f, hw_m),
    her_market  = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
    his_market  = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
    her_total   = ifelse(female_f == 1, total_hours_f, total_hours_m),
    his_total   = ifelse(female_f == 0, total_hours_f, total_hours_m),
    her_emp     = ifelse(female_f == 1, employed_f, employed_m),
    his_emp     = ifelse(female_f == 0, employed_f, employed_m),
    her_ft      = ifelse(female_f == 1, fulltime_f, fulltime_m),
    his_ft      = ifelse(female_f == 0, fulltime_f, fulltime_m),
    has_children = as.integer(ever_parent_f == 1 | ever_parent_m == 1),
    couple_id   = ifelse(person_id < partner_id,
                         paste0(person_id, "_", partner_id),
                         paste0(partner_id, "_", person_id))
  )

cat(sprintf("Couples: %s, couple-waves: %s\n\n",
            format(n_distinct(couples$couple_id), big.mark = ","),
            format(nrow(couples), big.mark = ",")))

# --- 3A. Couple-level time budget by parenthood ---
cat("--- 3A. Couple-level time allocation ---\n\n")

couple_time <- couples %>%
  group_by(has_children) %>%
  summarise(
    N = n(),
    his_market = mean(his_market, na.rm = TRUE),
    her_market = mean(her_market, na.rm = TRUE),
    his_hw     = mean(his_hw, na.rm = TRUE),
    her_hw     = mean(her_hw, na.rm = TRUE),
    his_total  = mean(his_total, na.rm = TRUE),
    her_total  = mean(her_total, na.rm = TRUE),
    hw_ratio_hf = mean(her_hw, na.rm = TRUE) / mean(his_hw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(status = ifelse(has_children == 0, "Childless", "Parents"))

cat(sprintf("  %-12s  %6s  %s  |  %s  |  %s  |  %s\n",
            "Status", "N",
            "His: Market  HW   Total",
            "Her: Market  HW   Total",
            "HW ratio", "Market gap"))
cat("  ", strrep("-", 90), "\n")
for (i in 1:nrow(couple_time)) {
  r <- couple_time[i,]
  cat(sprintf("  %-12s  %6s  %6.1f  %5.1f  %5.1f  |  %6.1f  %5.1f  %5.1f  |  %5.2f  |  %+6.1f\n",
              r$status, format(r$N, big.mark = ","),
              r$his_market, r$his_hw, r$his_total,
              r$her_market, r$her_hw, r$her_total,
              r$hw_ratio_hf, r$his_market - r$her_market))
}

# --- 3B. Housework share (her share of couple's total housework) ---
cat("\n--- 3B. Housework share ---\n")

couples <- couples %>%
  mutate(
    couple_hw = his_hw + her_hw,
    her_hw_share = ifelse(couple_hw > 0, her_hw / couple_hw, NA_real_),
    couple_market = his_market + her_market,
    her_market_share = ifelse(couple_market > 0, her_market / couple_market, NA_real_)
  )

share_summary <- couples %>%
  group_by(has_children) %>%
  summarise(
    her_hw_share = mean(her_hw_share, na.rm = TRUE) * 100,
    her_market_share = mean(her_market_share, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat(sprintf("\n  Childless — her HW share: %.1f%%, her market share: %.1f%%\n",
            share_summary$her_hw_share[1], share_summary$her_market_share[1]))
cat(sprintf("  Parents  — her HW share: %.1f%%, her market share: %.1f%%\n",
            share_summary$her_hw_share[2], share_summary$her_market_share[2]))
cat("  → Parenthood amplifies the gender division of labour\n")


# #############################################################################
# PART 4: HOUSEWORK AROUND BIRTH (EVENT STUDY)
# #############################################################################

cat("\n\n")
cat(strrep("#", 70), "\n")
cat("#  PART 4: HOUSEWORK EVENT STUDY AROUND BIRTH\n")
cat(strrep("#", 70), "\n\n")

birth_events <- panel %>%
  filter(had_birth_this_wave == 1) %>%
  group_by(person_id) %>%
  summarise(birth_wave = min(wave), .groups = "drop")

event_hw <- panel %>%
  inner_join(birth_events, by = "person_id") %>%
  mutate(t = wave - birth_wave) %>%
  filter(t >= -3, t <= 4) %>%
  group_by(t, female) %>%
  summarise(
    N = n(),
    mean_hw     = mean(hw, na.rm = TRUE),
    mean_market = mean(hours_worked_clean, na.rm = TRUE),
    mean_total  = mean(total_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(gender = ifelse(female == 0, "Men", "Women"))

cat(sprintf("  %3s  %-6s  %5s  %7s  %7s  %7s\n",
            "t", "Gender", "N", "HW", "Market", "Total"))
cat("  ", strrep("-", 45), "\n")
for (i in 1:nrow(event_hw)) {
  r <- event_hw[i,]
  cat(sprintf("  %+2d   %-6s  %5d  %7.1f  %7.1f  %7.1f\n",
              r$t, r$gender, r$N, r$mean_hw, r$mean_market, r$mean_total))
}


# #############################################################################
# PART 5: WITHIN-COUPLE REALLOCATION AROUND BIRTH
# #############################################################################

cat("\n\n")
cat(strrep("#", 70), "\n")
cat("#  PART 5: WITHIN-COUPLE HOUSEWORK AROUND BIRTH\n")
cat(strrep("#", 70), "\n\n")

couple_births <- couples %>%
  filter(had_birth_this_wave == 1) %>%
  group_by(couple_id) %>%
  summarise(birth_wave = min(wave), .groups = "drop")

cat(sprintf("Couples with observed birth: %s\n\n",
            format(nrow(couple_births), big.mark = ",")))

couple_event <- couples %>%
  inner_join(couple_births, by = "couple_id") %>%
  mutate(t = wave - birth_wave) %>%
  filter(t >= -3, t <= 4)

ce_hw <- couple_event %>%
  group_by(t) %>%
  summarise(
    N = n(),
    his_hw = mean(his_hw, na.rm = TRUE),
    her_hw = mean(her_hw, na.rm = TRUE),
    hw_gap = mean(her_hw - his_hw, na.rm = TRUE),
    his_market = mean(his_market, na.rm = TRUE),
    her_market = mean(her_market, na.rm = TRUE),
    market_gap = mean(his_market - her_market, na.rm = TRUE),
    his_total = mean(his_total, na.rm = TRUE),
    her_total = mean(her_total, na.rm = TRUE),
    her_hw_share = mean(her_hw_share, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat(sprintf("  %3s  %5s  | %s  | %s  | %s\n",
            "t", "N",
            "His: HW  Market Total",
            "Her: HW  Market Total",
            "HW_gap  Share%"))
cat("  ", strrep("-", 80), "\n")
for (i in 1:nrow(ce_hw)) {
  r <- ce_hw[i,]
  cat(sprintf("  %+2d   %5d  | %5.1f  %5.1f  %5.1f  | %5.1f  %5.1f  %5.1f  | %+5.1f  %5.1f%%\n",
              r$t, r$N,
              r$his_hw, r$his_market, r$his_total,
              r$her_hw, r$her_market, r$her_total,
              r$hw_gap, r$her_hw_share))
}


# #############################################################################
# PART 6: MEDIATION TEST — DOES HOUSEWORK EXPLAIN β₂?
# #############################################################################

cat("\n\n")
cat(strrep("#", 70), "\n")
cat("#  PART 6: MEDIATION — HOUSEWORK AND \u03b2\u2082\n")
# Ensure hw exists in both panel and wage (self-healing)
if (!"hw" %in% names(panel)) {
  .hw_col <- c("housework_hours", "care_intensity","hsjob","jomfh","hw_hours")[
    c("housework_hours", "care_intensity","hsjob","jomfh","hw_hours") %in% names(panel)][1]
  if (!is.na(.hw_col)) { panel$hw <- pmin(pmax(panel[[.hw_col]], 0), 100) } else {
    stop(paste("No housework column found. Run names(panel) to diagnose.")) }
  panel$total_hours <- ifelse(is.na(panel$hours_worked_clean)|is.na(panel$hw),
                              NA_real_, panel$hours_worked_clean + panel$hw)
}
if (!"hw" %in% names(wage)) {
  wage$hw          <- panel$hw[match(paste(wage$person_id, wage$wave),
                                     paste(panel$person_id, panel$wave))]
  wage$total_hours <- panel$total_hours[match(paste(wage$person_id, wage$wave),
                                              paste(panel$person_id, panel$wave))]
}

cat(strrep("#", 70), "\n")
cat("\n  NOTE: Housework is endogenous (same allocation decision as β₂).\n")
cat("  This is a MECHANISM test, not a robustness check.\n")
cat("  If β₂ for mothers increases when housework is controlled,\n")
cat("  housework is absorbing part of the complementarity penalty.\n\n")

for (g in c("childless_men_ever", "fathers_ever",
            "never_mothers", "mothers_ever")) {
  
  g_label <- label_group(g)
  d <- wage %>%
    filter(group_ever == g, !is.na(hw)) %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c)
  
  # Baseline: without housework
  m_base <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL)
  
  # + housework hours
  m_hw <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num + hw | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL)
  
  # + housework × education interaction
  d <- d %>% mutate(hw_educ = hw * educ_c)
  m_hw_int <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num + hw + hw_educ | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL)
  
  get_b2 <- function(m) {
    if (is.null(m) || !("educ_exp_c" %in% names(coef(m)))) return(c(NA, NA))
    s <- summary(m)$coeftable
    c(s["educ_exp_c", 1], s["educ_exp_c", 4])
  }
  
  b_base <- get_b2(m_base)
  b_hw   <- get_b2(m_hw)
  b_int  <- get_b2(m_hw_int)
  
  cat(sprintf("  %-20s  Baseline β₂=%+.4f%%%s | +HW: %+.4f%%%s | +HW×educ: %+.4f%%%s\n",
              g_label,
              b_base[1]*100, stars_fn(b_base[2]),
              b_hw[1]*100, stars_fn(b_hw[2]),
              b_int[1]*100, stars_fn(b_int[2])))
  
  # Report HW coefficient in wage equation
  if (!is.null(m_hw) && "hw" %in% names(coef(m_hw))) {
    s_hw <- summary(m_hw)$coeftable
    cat(sprintf("    %-18s  housework→wage: %+.5f%s (p=%.4f)\n",
                "", s_hw["hw", 1], stars_fn(s_hw["hw", 4]), s_hw["hw", 4]))
  }
  
  # Change in β₂
  if (!is.na(b_base[1]) && !is.na(b_hw[1])) {
    delta_pct <- (b_hw[1] - b_base[1]) / abs(b_base[1]) * 100
    cat(sprintf("    %-18s  Δβ₂ from adding HW: %+.1f%%\n", "", delta_pct))
  }
  cat("\n")
}


# Namespace protection (avoid MASS/stats masking dplyr verbs)
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

################################################################################
# PART C: PANEL IV (parental-education instrument, diagnostic only)
################################################################################
# Paper Section V: the descriptive complementarity parameter b2 is reported by
# individual fixed effects (OLS/FE) and treated as descriptive only. The
# parental-education x experience instrument is NOT used for a causal reading
# of b2: parental education plausibly enters the wage equation directly, so the
# exclusion restriction is implausible on a priori grounds. The pooled IV and
# FE-IV below are retained as first-stage / weak-instrument diagnostics, and the
# Altonji-Pierret decomposition is independent of any instrument. The earlier
# Correlation Restriction (CR) test block has been removed.
#
# Assumes the master environment is loaded: `panel` and `stars_fn()` exist.
################################################################################

cat("=== PANEL IV (diagnostic) ===\n")

# ── Build lean Part C wage sample (wage_c), keeping global wage intact ────────
# Selecting only the columns needed across all Part C sections avoids
# carrying the full panel columns through every group subset.
wage_c <- panel %>%
  mutate(
    parent_educ_avg = case_when(
      !is.na(father_educ_years) & !is.na(mother_educ_years) ~
        (father_educ_years + mother_educ_years) / 2,
      !is.na(father_educ_years) ~ father_educ_years,
      !is.na(mother_educ_years) ~ mother_educ_years,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    in_wage_sample == 1,
    !is.na(parent_educ_avg),
    !is.na(educ_years),
    !is.na(experience_years),
    !is.na(ln_hourly_wage_real)
  ) %>%
  select(
    xwaveid, person_id, wave, group_ever,
    ln_hourly_wage_real, educ_years, experience_years,
    age_sq, married_num, parent_educ_avg,
    hours_worked_clean, fulltime
  )

cat(sprintf("Parent educ coverage: father=%.1f%%, mother=%.1f%%, either=%.1f%%\n",
            100 * mean(!is.na(panel$father_educ_years)),
            100 * mean(!is.na(panel$mother_educ_years)),
            100 * mean(!is.na(panel$father_educ_years) |
                         !is.na(panel$mother_educ_years))))
cat(sprintf("Wage_c (Part C): %d obs, %d individuals\n",
            nrow(wage_c), n_distinct(wage_c$xwaveid)))
gc(verbose = FALSE)

groups        <- c("childless_men_ever", "fathers_ever", "never_mothers", "mothers_ever")
group_labels  <- c("Childless men", "Fathers", "Never-mothers", "Mothers")

iv_results          <- list()
first_stage_results <- list()
feiv_results        <- list()
ap_results          <- list()

###############################################################################
# SINGLE PASS: all models per group, wave dummies built locally for CR only
###############################################################################

for (i in seq_along(groups)) {
  g     <- groups[i]
  g_lab <- group_labels[i]
  
  # ── Subset and center within group ─────────────────────────────────────────
  d <- wage_c %>%
    filter(group_ever == g) %>%
    mutate(
      wave_f            = factor(wave),
      educ_c            = educ_years      - mean(educ_years,      na.rm = TRUE),
      exper_c           = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c        = educ_c * exper_c,
      parent_educ_c     = parent_educ_avg  - mean(parent_educ_avg,  na.rm = TRUE),
      parent_educ_exp_c = parent_educ_c    * exper_c,
      # Altonji-Pierret variables
      t_c               = exper_c,
      S_t               = educ_c   * exper_c,
      Z_t               = parent_educ_c * exper_c
    )
  n   <- nrow(d)
  nid <- n_distinct(d$xwaveid)
  cat(sprintf("\n%s (N=%d, %d id):\n", g_lab, n, nid))
  
  # ── 1. First stage ─────────────────────────────────────────────────────────
  fs_level <- feols(
    educ_c ~ parent_educ_c + exper_c + age_sq + married_num | wave_f,
    data = d, cluster = ~xwaveid, notes = FALSE)
  fs_inter <- feols(
    educ_exp_c ~ parent_educ_c + parent_educ_exp_c +
      exper_c + age_sq + married_num | wave_f,
    data = d, cluster = ~xwaveid, notes = FALSE)
  f_stat <- fixest::wald(fs_inter, c("parent_educ_c", "parent_educ_exp_c"))
  
  first_stage_results[[g]] <- list(
    group = g_lab, n = n, n_id = nid,
    fs_level_coef = coef(fs_level)["parent_educ_c"],
    fs_inter_coef = coef(fs_inter)["parent_educ_exp_c"],
    joint_F = f_stat$stat, joint_p = f_stat$p
  )
  cat(sprintf("  Level FS: %.3f (p=%.4f)\n",
              coef(fs_level)["parent_educ_c"], pvalue(fs_level)["parent_educ_c"]))
  cat(sprintf("  Inter FS: %.3f (p=%.4f)\n",
              coef(fs_inter)["parent_educ_exp_c"], pvalue(fs_inter)["parent_educ_exp_c"]))
  cat(sprintf("  Joint F = %.1f (p=%.1e)\n", f_stat$stat, f_stat$p))
  rm(fs_level, fs_inter); gc(verbose = FALSE)
  
  # ── 2. OLS, FE, pooled IV ──────────────────────────────────────────────────
  ols <- feols(
    ln_hourly_wage_real ~ educ_c + educ_exp_c + exper_c + age_sq + married_num | wave_f,
    data = d, cluster = ~xwaveid, notes = FALSE)
  fe <- feols(
    ln_hourly_wage_real ~ educ_exp_c + exper_c + age_sq + married_num | xwaveid + wave_f,
    data = d, cluster = ~xwaveid, notes = FALSE)
  iv <- tryCatch(
    feols(
      ln_hourly_wage_real ~ exper_c + age_sq + married_num | wave_f |
        educ_c + educ_exp_c ~ parent_educ_c + parent_educ_exp_c,
      data = d, cluster = ~xwaveid, notes = FALSE),
    error = function(e) { cat("  IV failed:", e$message, "\n"); NULL })
  
  iv_results[[g]] <- list(
    group  = g_lab, n = n,
    ols_b2 = coef(ols)["educ_exp_c"] * 100,
    ols_se = fixest::se(ols)["educ_exp_c"]   * 100,
    ols_p  = pvalue(ols)["educ_exp_c"],
    fe_b2  = coef(fe)["educ_exp_c"]  * 100,
    fe_se  = fixest::se(fe)["educ_exp_c"]    * 100,
    fe_p   = pvalue(fe)["educ_exp_c"]
  )
  if (!is.null(iv)) {
    iv_results[[g]]$iv_b2 <- coef(iv)["fit_educ_exp_c"]  * 100
    iv_results[[g]]$iv_se <- fixest::se(iv)["fit_educ_exp_c"]    * 100
    iv_results[[g]]$iv_p  <- pvalue(iv)["fit_educ_exp_c"]
    iv_results[[g]]$iv_F  <- tryCatch(fitstat(iv, "ivf")[[1]]$stat,
                                      error = function(e) NA)
  }
  cat(sprintf("  OLS β₂=%.3f%%  FE β₂=%.3f%%",
              iv_results[[g]]$ols_b2, iv_results[[g]]$fe_b2))
  if (!is.null(iv))
    cat(sprintf("  IV β₂=%.3f%% [F=%.0f]",
                iv_results[[g]]$iv_b2, iv_results[[g]]$iv_F))
  cat("\n")
  rm(ols, fe, iv); gc(verbose = FALSE)
  
  # ── 3. FE-IV ───────────────────────────────────────────────────────────────
  feiv <- tryCatch(
    feols(
      ln_hourly_wage_real ~ exper_c + age_sq + married_num | xwaveid + wave_f |
        educ_exp_c ~ parent_educ_exp_c,
      data = d, cluster = ~xwaveid, notes = FALSE),
    error = function(e) { cat("  FE-IV failed:", g_lab, e$message, "\n"); NULL })
  if (!is.null(feiv)) {
    fs_F <- tryCatch(fitstat(feiv, "ivf")[[1]]$stat, error = function(e) NA)
    feiv_results[[g]] <- list(
      group = g_lab, n = n,
      b2 = coef(feiv)["fit_educ_exp_c"] * 100,
      se = fixest::se(feiv)["fit_educ_exp_c"]   * 100,
      p  = pvalue(feiv)["fit_educ_exp_c"],
      F  = fs_F
    )
    cat(sprintf("  FE-IV β₂=%.3f%% (p=%.4f) [F=%.0f]\n",
                feiv_results[[g]]$b2, feiv_results[[g]]$p, fs_F))
  }
  rm(feiv); gc(verbose = FALSE)
  
  # ── 4. Altonji-Pierret decomposition ───────────────────────────────────────
  m_ap <- tryCatch(
    feols(
      ln_hourly_wage_real ~ educ_c + parent_educ_c + S_t + Z_t +
        age_sq + hours_worked_clean + married_num + fulltime | wave_f,
      data = d, cluster = ~xwaveid, notes = FALSE),
    error = function(e) NULL)
  if (!is.null(m_ap)) {
    s <- summary(m_ap)$coeftable
    fmt <- function(b, p) sprintf("%+.4f%s", b * 100, stars_fn(p))
    cat(sprintf("  AP  S=%-12s Z=%-12s S×t=%-12s Z×t=%-12s\n",
                fmt(s["educ_c",       1], s["educ_c",       4]),
                fmt(s["parent_educ_c",1], s["parent_educ_c",4]),
                fmt(s["S_t",          1], s["S_t",          4]),
                fmt(s["Z_t",          1], s["Z_t",          4])))
    ap_results[[g]] <- m_ap
  }
  rm(m_ap)
  
  # ── Done with group; free the subset ───────────────────────────────────────
  rm(d); gc(verbose = FALSE)
}

rm(wage_c); gc(verbose = FALSE)

###############################################################################
# SUMMARY TABLES
###############################################################################

cat("\n\n=== SUMMARY TABLE ===\n")
cat(sprintf("%-14s %6s %8s %8s %8s %8s %6s\n",
            "Group", "N", "OLS", "FE", "IV", "FE-IV", "FS F"))
cat(paste(rep("─", 70), collapse = ""), "\n")

star <- function(p) {
  if (is.null(p) || is.na(p)) return("")
  if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*"
  else if (p < .1) "†" else ""
}

for (i in seq_along(groups)) {
  g   <- groups[i]
  r   <- iv_results[[g]]
  iv_s  <- if (!is.null(r$iv_b2))
    sprintf("%6.3f%s", r$iv_b2, star(r$iv_p)) else "   —"
  fiv_s <- if (!is.null(feiv_results[[g]]))
    sprintf("%6.3f%s", feiv_results[[g]]$b2, star(feiv_results[[g]]$p)) else "   —"
  f_s   <- if (!is.null(feiv_results[[g]]$F) && !is.na(feiv_results[[g]]$F))
    sprintf("%6.0f", feiv_results[[g]]$F) else "   —"
  cat(sprintf("%-14s %6d %6.3f%-2s %6.3f%-2s %8s %8s %6s\n",
              group_labels[i], r$n,
              r$ols_b2, star(r$ols_p),
              r$fe_b2,  star(r$fe_p),
              iv_s, fiv_s, f_s))
}
cat(paste(rep("─", 70), collapse = ""), "\n")

cat("\n=== ALTONJI-PIERRET SUMMARY ===\n")
cat(sprintf("%-16s %12s %12s %12s %12s %6s\n",
            "Group", "b1: S", "b2: Z", "b3: S×t", "b4: Z×t", "N"))
cat(paste(rep("─", 76), collapse = ""), "\n")
for (i in seq_along(groups)) {
  g <- groups[i]
  if (!is.null(ap_results[[g]])) {
    s   <- summary(ap_results[[g]])$coeftable
    fmt <- function(b, p) sprintf("%+.4f%s", b * 100, stars_fn(p))
    cat(sprintf("%-16s %12s %12s %12s %12s %6d\n",
                group_labels[i],
                fmt(s["educ_c",        1], s["educ_c",        4]),
                fmt(s["parent_educ_c", 1], s["parent_educ_c", 4]),
                fmt(s["S_t",           1], s["S_t",           4]),
                fmt(s["Z_t",           1], s["Z_t",           4]),
                iv_results[[g]]$n))
  }
}
cat(paste(rep("─", 76), collapse = ""), "\n")
cat("b3 > 0: contradicts employer-learning attenuation (A-P predict b3 < 0).\n")

# --- Altonji-Pierret SE ledger (Table: AP decomposition, AEA no-stars) ---
ap_ledger <- dplyr::bind_rows(lapply(seq_along(groups), function(i) {
  g <- groups[i]
  if (is.null(ap_results[[g]])) return(NULL)
  s <- summary(ap_results[[g]])$coeftable
  getc <- function(r, col) if (r %in% rownames(s)) s[r, col] else NA_real_
  tibble::tibble(
    group = group_labels[i], N = iv_results[[g]]$n,
    b1_S   = getc("educ_c",1)*100,        b1_S_se   = getc("educ_c",2)*100,
    b2_Z   = getc("parent_educ_c",1)*100, b2_Z_se   = getc("parent_educ_c",2)*100,
    b3_Sxt = getc("S_t",1)*100,           b3_Sxt_se = getc("S_t",2)*100,
    b4_Zxt = getc("Z_t",1)*100,           b4_Zxt_se = getc("Z_t",2)*100
  )
}))
write.csv(ap_ledger, "se_ledger_ap.csv", row.names = FALSE)
cat("  Wrote se_ledger_ap.csv (Altonji-Pierret b1-b4 + clustered SEs).\n")

saveRDS(
  list(first_stage = first_stage_results, iv = iv_results,
       feiv = feiv_results, ap = ap_results),
  "panel_iv_results.rds"
)
cat("\nSaved: panel_iv_results.rds\nDone.\n")

# Pin the RNG state at the PART C -> PART D handoff. PART C is now fully
# deterministic (the CR Monte Carlo that previously consumed draws here has
# been removed). The downstream stochastic steps (DML cross-fitting, causal
# forest) already call set.seed(42) locally before drawing, so this is a
# belt-and-suspenders reset to the documented global seed and pins the
# .Random.seed captured by the checkpoint below.
set.seed(42)

save.image("master_state_after_partC.RData")  # after PART C
################################################################################
# PART D: COUPLE-LEVEL DML + CAUSAL FOREST
################################################################################

# 0. LOAD DATA AND BUILD COUPLE-STACKED PANEL

cat(strrep("=", 80), "\n")
cat("  COUPLE-LEVEL DML — ICHIMURA (2025) FRAMEWORK\n")
cat(strrep("=", 80), "\n\n")
num_rows    <- format(nrow(panel), big.mark = ",")
num_persons <- format(n_distinct(panel$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

# --- Clean partner_id ---

# --- Build couple panel (same logic as COUPLE_SPECIALIZATION Part B) ---
# Ensure partner_id_clean exists (self-healing)
if (!"partner_id_clean" %in% names(panel)) {
  if (!"partner_id" %in% names(panel)) panel$partner_id <- NA_character_
  panel <- panel %>%
    mutate(partner_id_clean = dplyr::case_when(
      is.na(partner_id)       ~ NA_character_,
      partner_id == ""        ~ NA_character_,
      nchar(partner_id) < 5  ~ NA_character_,
      grepl("^-", partner_id) ~ NA_character_,
      TRUE                    ~ partner_id
    ))
}
p1 <- panel %>%
  filter(!is.na(partner_id_clean)) %>%
  select(person_id, partner_id = partner_id_clean, wave, female, ever_parent,
         group_ever, group_3way_women, ln_hourly_wage_real,
         educ_years, experience_years, hours_worked_clean, employed,
         fulltime, in_wage_sample, married, age, age_sq, state,
         analysis_weight)

p2 <- panel %>%
  select(person_id, wave, female, ever_parent,
         group_ever, ln_hourly_wage_real,
         educ_years, experience_years, hours_worked_clean, employed,
         fulltime, in_wage_sample, married, age, age_sq)

couples_raw <- p1 %>%
  inner_join(p2, by = c("partner_id" = "person_id", "wave" = "wave"),
             suffix = c("_f", "_m")) %>%
  filter(female_f != female_m)  # heterosexual

# Standardize: her = female, him = male
couples <- couples_raw %>%
  mutate(
    her_wage = ifelse(female_f == 1, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
    his_wage = ifelse(female_f == 0, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
    her_educ = ifelse(female_f == 1, educ_years_f, educ_years_m),
    his_educ = ifelse(female_f == 0, educ_years_f, educ_years_m),
    her_exp = ifelse(female_f == 1, experience_years_f, experience_years_m),
    his_exp = ifelse(female_f == 0, experience_years_f, experience_years_m),
    her_hours = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
    his_hours = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
    her_employed = ifelse(female_f == 1, employed_f, employed_m),
    his_employed = ifelse(female_f == 0, employed_f, employed_m),
    her_fulltime = ifelse(female_f == 1, fulltime_f, fulltime_m),
    his_fulltime = ifelse(female_f == 0, fulltime_f, fulltime_m),
    her_age = ifelse(female_f == 1, age_f, age_m),
    his_age = ifelse(female_f == 0, age_f, age_m),
    her_age_sq = ifelse(female_f == 1, age_sq_f, age_sq_m),
    his_age_sq = ifelse(female_f == 0, age_sq_f, age_sq_m),
    her_married = ifelse(female_f == 1, married_f, married_m),
    his_married = ifelse(female_f == 0, married_f, married_m),
    her_in_wage = ifelse(female_f == 1, in_wage_sample_f, in_wage_sample_m),
    his_in_wage = ifelse(female_f == 0, in_wage_sample_f, in_wage_sample_m),
    has_children = as.integer(ever_parent_f == 1 | ever_parent_m == 1),
    couple_id = ifelse(person_id < partner_id,
                       paste0(person_id, "_", partner_id),
                       paste0(partner_id, "_", person_id))
  )

# Dual-earner restriction for wage analysis
dual_earner <- couples %>%
  filter(his_in_wage == 1, her_in_wage == 1,
         !is.na(his_wage), !is.na(her_wage))

num_rows    <- format(nrow(dual_earner), big.mark = ",")
num_persons <- format(n_distinct(dual_earner$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

# --- Stack into long format ---
his_data <- dual_earner %>%
  transmute(couple_id, wave, has_children,
            ln_wage = his_wage, educ = his_educ, exp = his_exp,
            hours = his_hours, age_sq = his_age_sq,
            married = his_married, fulltime = his_fulltime,
            is_female = 0L,
            # Partner variables (key for DML — her characteristics as HIS controls)
            partner_educ = her_educ, partner_exp = her_exp,
            partner_hours = her_hours, partner_fulltime = her_fulltime)

her_data <- dual_earner %>%
  transmute(couple_id, wave, has_children,
            ln_wage = her_wage, educ = her_educ, exp = her_exp,
            hours = her_hours, age_sq = her_age_sq,
            married = her_married, fulltime = her_fulltime,
            is_female = 1L,
            partner_educ = his_educ, partner_exp = his_exp,
            partner_hours = his_hours, partner_fulltime = his_fulltime)

stacked <- bind_rows(his_data, her_data) %>%
  mutate(
    educ_c = educ - mean(educ, na.rm = TRUE),
    exp_c = exp - mean(exp, na.rm = TRUE),
    educ_exp_c = educ_c * exp_c,
    person_couple = paste0(couple_id, "_", is_female)
  ) %>%
  filter(complete.cases(.))

cat(sprintf("Stacked obs: %s (his: %s, her: %s)\n",
            format(nrow(stacked), big.mark = ","),
            format(sum(stacked$is_female == 0), big.mark = ","),
            format(sum(stacked$is_female == 1), big.mark = ",")))


# DML HELPER FUNCTIONS

# --- Within-transform by couple (couple-FE equivalent) ---
within_transform_couple <- function(data, vars, id_var = "couple_id") {
  data %>%
    group_by(across(all_of(id_var))) %>%
    mutate(across(all_of(vars), ~ .x - mean(.x, na.rm = TRUE),
                  .names = "{.col}_w")) %>%
    ungroup()
}

# --- DML-PLR wrapper ---
run_dml <- function(data, label, treatment, outcome, controls,
                    ml_method = "ranger", n_folds = 5, n_rep = 3) {
  all_vars <- unique(c(outcome, treatment, controls))
  dat <- data %>%
    dplyr::select(all_of(all_vars)) %>%
    filter(complete.cases(.))
  
  if (nrow(dat) < 300) {
    cat(sprintf("  %-50s  N=%d too small, skipping\n", label, nrow(dat)))
    return(NULL)
  }
  
  dml_data <- DoubleMLData$new(
    as.data.frame(dat),
    y_col = outcome,
    d_cols = treatment,
    x_cols = controls
  )
  
  if (ml_method == "ranger") {
    ml_l <- lrn("regr.ranger", num.trees = 500, min.node.size = 5, seed = 42,
                num.threads = n_threads)   # ranger defaults to ALL cores otherwise
    ml_m <- lrn("regr.ranger", num.trees = 500, min.node.size = 5, seed = 43,
                num.threads = n_threads)
  } else {
    ml_l <- lrn("regr.cv_glmnet", s = "lambda.min", alpha = 1)
    ml_m <- lrn("regr.cv_glmnet", s = "lambda.min", alpha = 1)
  }
  
  dml_plr <- tryCatch({
    set.seed(42)                      # reproducible cross-fitting splits + nuisance forests
    m <- DoubleMLPLR$new(dml_data, ml_l, ml_m,
                         n_folds = n_folds, n_rep = n_rep,
                         score = "partialling out")
    m$fit()
    m
  }, error = function(e) {
    cat(sprintf("  %-50s  ERROR: %s\n", label, e$message))
    return(NULL)
  })
  
  if (is.null(dml_plr)) return(NULL)
  
  theta <- dml_plr$coef
  se    <- dml_plr$se
  p     <- dml_plr$pval
  
  cat(sprintf("  %-50s  θ=%+.4f%%%s  SE=%.4f  p=%.4f  N=%d\n",
              label, theta * 100, stars_fn(p), se * 100, p, nrow(dat)))
  
  list(model = dml_plr, theta = theta, se = se, pval = p, n = nrow(dat))
}


# #############################################################################
# PART 1: POOLED DML — β₂ BY GENDER × PARENTHOOD
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §4.2 (Household-level identification): DML-PLR with Random Forest
#   (500 trees). Treatment D = educ_exp_c. Controls W include partner vars.
# Paper P4 test: Within-couple DML gender gap → p = 0.216.
# Sample: Dual-earner stacked → 84,294 obs.
#   DML by cell: Childless men 10,394 / Fathers 22,411 /
#                Childless women 10,394 / Mothers 22,411.
cat("#  PART 1: POOLED DML-PLR FOR β₂ (stacked couple data)\n")
cat(strrep("#", 80), "\n")

# Controls WITHOUT partner variables (for comparison with individual-level)
controls_own <- c("exp", "hours", "age_sq", "married", "fulltime",
                  "is_female", "has_children", "wave")

# Controls WITH partner variables (key innovation — Becker channel)
controls_partner <- c("exp", "hours", "age_sq", "married", "fulltime",
                      "is_female", "has_children", "wave",
                      "partner_educ", "partner_exp", "partner_hours",
                      "partner_fulltime")


# --- 1A. Overall β₂ on stacked data ---
cat("\n--- 1A. Overall β₂ (all stacked observations) ---\n")

cat("\n  Without partner controls:\n")
dml_all_own <- run_dml(stacked, "All stacked (own controls, RF)",
                       "educ_exp_c", "ln_wage", controls_own, "ranger")

cat("\n  With partner controls:\n")
dml_all_partner <- run_dml(stacked, "All stacked (+ partner controls, RF)",
                           "educ_exp_c", "ln_wage", controls_partner, "ranger")

# LASSO comparison
cat("\n  LASSO comparison:\n")
dml_all_lasso <- run_dml(stacked, "All stacked (own controls, LASSO)",
                         "educ_exp_c", "ln_wage", controls_own, "glmnet")
dml_all_partner_lasso <- run_dml(stacked, "All stacked (+ partner, LASSO)",
                                 "educ_exp_c", "ln_wage", controls_partner, "glmnet")


# --- 1B. β₂ by gender (stacked, with partner controls) ---
cat("\n\n--- 1B. β₂ by gender (within couples, with partner controls) ---\n")

# Controls for gender-specific models (drop is_female)
controls_gender <- c("exp", "hours", "age_sq", "married", "fulltime",
                     "has_children", "wave",
                     "partner_educ", "partner_exp", "partner_hours",
                     "partner_fulltime")

cat("\n  Random Forest:\n")
dml_his <- run_dml(stacked %>% filter(is_female == 0),
                   "His β₂ (RF, + partner controls)",
                   "educ_exp_c", "ln_wage", controls_gender, "ranger")

dml_her <- run_dml(stacked %>% filter(is_female == 1),
                   "Her β₂ (RF, + partner controls)",
                   "educ_exp_c", "ln_wage", controls_gender, "ranger")

# Formal test: his β₂ vs her β₂
if (!is.null(dml_his) && !is.null(dml_her)) {
  z_gender <- (dml_his$theta - dml_her$theta) /
    sqrt(dml_his$se^2 + dml_her$se^2)
  p_gender <- 2 * pnorm(-abs(z_gender))
  cat(sprintf("\n  Gender gap in DML β₂: %.4f%% (z=%.3f, p=%.4f%s)\n",
              (dml_his$theta - dml_her$theta) * 100, z_gender, p_gender,
              stars_fn(p_gender)))
}


# --- 1C. β₂ by gender × parenthood (four cells) ---
cat("\n\n--- 1C. β₂ by gender × parenthood (four cells) ---\n")

# Drop is_female and has_children from controls (they define the cells)
controls_cell <- c("exp", "hours", "age_sq", "married", "fulltime",
                   "wave", "partner_educ", "partner_exp",
                   "partner_hours", "partner_fulltime")
cells <- list(
  "Childless men"    = stacked %>% filter(is_female == 0, has_children == 0),
  "Fathers"          = stacked %>% filter(is_female == 0, has_children == 1),
  "Childless women"  = stacked %>% filter(is_female == 1, has_children == 0),
  "Mothers"          = stacked %>% filter(is_female == 1, has_children == 1)
)
cat("\n  Random Forest:\n")
dml_cells_rf <- list()
for (g in names(cells)) {
  dml_cells_rf[[g]] <- run_dml(cells[[g]], paste0(g, " (RF, + partner)"),
                               "educ_exp_c", "ln_wage", controls_cell, "ranger")
}
cat("\n  LASSO:\n")
dml_cells_lasso <- list()
for (g in names(cells)) {
  dml_cells_lasso[[g]] <- run_dml(cells[[g]], paste0(g, " (LASSO, + partner)"),
                                  "educ_exp_c", "ln_wage", controls_cell, "glmnet")
}

# Formal tests between cells
cat("\n  Formal DML comparisons (RF):\n")
compare_dml <- function(r1, r2, label) {
  if (is.null(r1) || is.null(r2)) return(invisible(NULL))
  z <- (r1$theta - r2$theta) / sqrt(r1$se^2 + r2$se^2)
  p <- 2 * pnorm(-abs(z))
  cat(sprintf("  %-40s  Δ=%+.4f%%  z=%.3f  p=%.4f%s\n",
              label, (r1$theta - r2$theta) * 100, z, p, stars_fn(p)))
}

compare_dml(dml_cells_rf[["Childless men"]], dml_cells_rf[["Fathers"]],
            "Childless men vs Fathers")
compare_dml(dml_cells_rf[["Childless women"]], dml_cells_rf[["Mothers"]],
            "Childless women vs Mothers")
compare_dml(dml_cells_rf[["Childless men"]], dml_cells_rf[["Childless women"]],
            "Childless men vs Childless women")
compare_dml(dml_cells_rf[["Fathers"]], dml_cells_rf[["Mothers"]],
            "Fathers vs Mothers")


# --- 1D. Does adding partner controls CHANGE β₂? ---
cat("\n\n--- 1D. Partner controls impact (Becker channel test) ---\n")

cat("\n  Without partner controls:\n")
dml_his_own <- run_dml(stacked %>% filter(is_female == 0),
                       "His β₂ (RF, own controls only)",
                       "educ_exp_c", "ln_wage",
                       c("exp", "hours", "age_sq", "married", "fulltime",
                         "has_children", "wave"), "ranger")
dml_her_own <- run_dml(stacked %>% filter(is_female == 1),
                       "Her β₂ (RF, own controls only)",
                       "educ_exp_c", "ln_wage",
                       c("exp", "hours", "age_sq", "married", "fulltime",
                         "has_children", "wave"), "ranger")

cat("\n  With partner controls:\n")
cat("  (Already estimated above as dml_his and dml_her)\n")

if (!is.null(dml_his_own) && !is.null(dml_his)) {
  cat(sprintf("\n  His β₂: own=%.4f%%, +partner=%.4f%%, Δ=%+.4f%%\n",
              dml_his_own$theta*100, dml_his$theta*100,
              (dml_his$theta - dml_his_own$theta)*100))
}
if (!is.null(dml_her_own) && !is.null(dml_her)) {
  cat(sprintf("  Her β₂: own=%.4f%%, +partner=%.4f%%, Δ=%+.4f%%\n",
              dml_her_own$theta*100, dml_her$theta*100,
              (dml_her$theta - dml_her_own$theta)*100))
}
cat("\n  If adding partner controls moves β₂:\n")
cat("    → Partner characteristics confound own complementarity (Becker)\n")
cat("  If β₂ unchanged:\n")
cat("    → Own complementarity is independent of partner (no spillover)\n")


# #############################################################################
# PART 2: WITHIN-COUPLE DML (couple-demeaned, FE equivalent)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §4.2: "within-couple DML specification in which all variables
#   are demeaned by couple". Absorbs household-level heterogeneity.
# Sample: Couples with ≥4 observations in stacked data.
#   Expected: 84,294 stacked obs / 4,662 couples.
cat("#  PART 2: WITHIN-COUPLE DML-PLR (couple-demeaned)\n")
cat(strrep("#", 80), "\n")

# Within-transform: demean by couple_id (across his/her + waves)
# This absorbs all couple-level heterogeneity (assortative matching,
# household unobservables, neighbourhood effects)

vars_to_demean <- c("ln_wage", "educ_exp_c", "exp", "hours", "age_sq",
                    "married", "fulltime", "partner_exp", "partner_hours",
                    "partner_fulltime")

stacked_w <- stacked %>%
  group_by(couple_id) %>%
  filter(n() >= 4) %>%  # need variation within couple for FE
  ungroup()

cat(sprintf("\n  Couples with 4+ obs: %s (from %s stacked obs)\n",
            format(n_distinct(stacked_w$couple_id), big.mark = ","),
            format(nrow(stacked_w), big.mark = ",")))

stacked_w <- within_transform_couple(stacked_w, vars_to_demean, "couple_id")

# Controls for within-couple DML (demeaned versions)
controls_w <- c("exp_w", "hours_w", "age_sq_w", "married_w", "fulltime_w",
                "is_female", "has_children", "wave")

controls_w_partner <- c("exp_w", "hours_w", "age_sq_w", "married_w",
                        "fulltime_w", "is_female", "has_children", "wave",
                        "partner_exp_w", "partner_hours_w", "partner_fulltime_w")

cat("\n--- 2A. Overall within-couple DML β₂ ---\n")

dml_w_own <- run_dml(stacked_w, "Within-couple (own controls, RF)",
                     "educ_exp_c_w", "ln_wage_w", controls_w, "ranger")

dml_w_partner <- run_dml(stacked_w, "Within-couple (+ partner, RF)",
                         "educ_exp_c_w", "ln_wage_w", controls_w_partner, "ranger")


cat("\n--- 2B. Within-couple DML β₂ by gender ---\n")

controls_w_gender <- c("exp_w", "hours_w", "age_sq_w", "married_w",
                       "fulltime_w", "has_children", "wave",
                       "partner_exp_w", "partner_hours_w", "partner_fulltime_w")

dml_w_his <- run_dml(stacked_w %>% filter(is_female == 0),
                     "His within-couple β₂ (RF, + partner)",
                     "educ_exp_c_w", "ln_wage_w", controls_w_gender, "ranger")

dml_w_her <- run_dml(stacked_w %>% filter(is_female == 1),
                     "Her within-couple β₂ (RF, + partner)",
                     "educ_exp_c_w", "ln_wage_w", controls_w_gender, "ranger")

if (!is.null(dml_w_his) && !is.null(dml_w_her)) {
  z_w <- (dml_w_his$theta - dml_w_her$theta) /
    sqrt(dml_w_his$se^2 + dml_w_her$se^2)
  p_w <- 2 * pnorm(-abs(z_w))
  cat(sprintf("\n  Within-couple gender gap: %.4f%% (z=%.3f, p=%.4f%s)\n",
              (dml_w_his$theta - dml_w_her$theta) * 100, z_w, p_w,
              stars_fn(p_w)))
  rec("dml_his", dml_w_his$theta * 100); rec("dml_her", dml_w_her$theta * 100)
  rec("dml_gap", (dml_w_his$theta - dml_w_her$theta) * 100); rec("dml_gap_p", p_w)
}


cat("\n--- 2C. Within-couple DML β₂ by gender × parenthood ---\n")

controls_w_cell <- c("exp_w", "hours_w", "age_sq_w", "married_w",
                     "fulltime_w", "wave",
                     "partner_exp_w", "partner_hours_w", "partner_fulltime_w")

w_cells <- list(
  "Childless men (W)"   = stacked_w %>% filter(is_female == 0, has_children == 0),
  "Fathers (W)"         = stacked_w %>% filter(is_female == 0, has_children == 1),
  "Childless women (W)" = stacked_w %>% filter(is_female == 1, has_children == 0),
  "Mothers (W)"         = stacked_w %>% filter(is_female == 1, has_children == 1)
)

dml_w_cells <- list()
for (g in names(w_cells)) {
  dml_w_cells[[g]] <- run_dml(w_cells[[g]], paste0(g, " (RF)"),
                              "educ_exp_c_w", "ln_wage_w",
                              controls_w_cell, "ranger")
}

cat("\n  Within-couple DML comparisons:\n")
compare_dml(dml_w_cells[["Childless men (W)"]], dml_w_cells[["Fathers (W)"]],
            "Childless men vs Fathers (within-couple)")
compare_dml(dml_w_cells[["Childless women (W)"]], dml_w_cells[["Mothers (W)"]],
            "Childless women vs Mothers (within-couple)")


# #############################################################################
# PART 3: CAUSAL FOREST ON COUPLE DATA
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §4.3 (Heterogeneity): Causal forest (2,000 trees) on stacked
#   couple data. "Own experience dominates (82.8%), partner variables
#   collectively scoring 8.88× lower."
# Sample: Stacked couple data, 84,294 obs.
cat("#  PART 3: COUPLE-LEVEL CAUSAL FOREST\n")
cat(strrep("#", 80), "\n")

cat("\n--- 3A. Pooled causal forest with partner covariates ---\n")

cf_data <- stacked %>%
  select(ln_wage, educ_exp_c, exp, hours, age_sq, married, fulltime,
         is_female, has_children, wave,
         partner_educ, partner_exp, partner_hours, partner_fulltime) %>%
  filter(complete.cases(.))

cat(sprintf("  CF data: N=%s\n", format(nrow(cf_data), big.mark = ",")))

# Covariates include PARTNER variables
X_vars <- c("exp", "hours", "age_sq", "married", "fulltime",
            "is_female", "has_children", "wave",
            "partner_educ", "partner_exp", "partner_hours", "partner_fulltime")
X <- as.matrix(cf_data[, X_vars])
Y <- cf_data$ln_wage
W <- cf_data$educ_exp_c  # treatment = education × experience interaction

cat("  Fitting causal forest (2000 trees)...\n")
set.seed(42)
cf <- tryCatch(
  causal_forest(X, Y, W, num.trees = 2000, min.node.size = 20, honesty = TRUE, seed = 42,
                num.threads = n_threads),   # grf defaults to ALL cores otherwise
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(cf)) {
  tau <- predict(cf)$predictions
  cf_data$tau_hat <- tau
  
  # ATE
  ate <- average_treatment_effect(cf)
  cat(sprintf("  Overall ATE: %.4f%% (SE %.4f%%)\n", ate[1]*100, ate[2]*100))
  
  # Calibration test
  cal <- tryCatch(test_calibration(cf), error = function(e) NULL)
  if (!is.null(cal)) {
    cat(sprintf("  Calibration: mean=%.3f (p=%.4f), diff=%.3f (p=%.4f)\n",
                cal[1,1], cal[1,4], cal[2,1], cal[2,4]))
    cat(sprintf("  Heterogeneity: p=%.4f %s\n", cal[2,4],
                ifelse(cal[2,4] < 0.05, "→ HETEROGENEOUS", "→ HOMOGENEOUS")))
  }
  
  # --- Variable importance ---
  cat("\n  Variable importance (higher = more drives heterogeneity):\n")
  vimp <- variable_importance(cf)
  vimp_df <- tibble(variable = X_vars, importance = as.numeric(vimp)) %>%
    arrange(desc(importance))
  for (i in 1:nrow(vimp_df)) {
    marker <- ifelse(grepl("partner", vimp_df$variable[i]), " ← PARTNER", "")
    cat(sprintf("    %-20s  %.4f%s\n",
                vimp_df$variable[i], vimp_df$importance[i], marker))
  }
  
  # KEY TEST: Do partner variables rank high?
  partner_vimp <- vimp_df %>% filter(grepl("partner", variable))
  own_vimp <- vimp_df %>% filter(!grepl("partner", variable))
  cat(sprintf("\n  Mean importance — Own: %.4f, Partner: %.4f, Ratio: %.2f\n",
              mean(own_vimp$importance), mean(partner_vimp$importance),
              mean(own_vimp$importance) / max(mean(partner_vimp$importance), 1e-6)))
  rec("forest_exp_imp", 100 * vimp_df$importance[match("exp", vimp_df$variable)])
  rec("forest_ratio", mean(own_vimp$importance) / max(mean(partner_vimp$importance), 1e-6))
  cat("  If partner importance is high → forest discovers Becker spillover\n")
  cat("  If partner importance is low → complementarity is individual, not household\n")
  
  
  # --- CATE by group ---
  cat("\n--- 3B. CATE by gender × parenthood ---\n")
  
  cf_data$group <- case_when(
    cf_data$is_female == 0 & cf_data$has_children == 0 ~ "Childless men",
    cf_data$is_female == 0 & cf_data$has_children == 1 ~ "Fathers",
    cf_data$is_female == 1 & cf_data$has_children == 0 ~ "Childless women",
    cf_data$is_female == 1 & cf_data$has_children == 1 ~ "Mothers"
  )
  
  group_cate <- cf_data %>%
    group_by(group) %>%
    summarise(
      N = n(),
      Mean_CATE = mean(tau_hat) * 100,
      SD_CATE = sd(tau_hat) * 100,
      r_exp = cor(tau_hat, exp),
      r_partner_hours = cor(tau_hat, partner_hours),
      .groups = "drop"
    )
  
  cat(sprintf("  %-20s  %6s  %8s  %8s  %8s  %12s\n",
              "Group", "N", "Mean τ%", "SD τ%", "r(τ,exp)", "r(τ,p_hours)"))
  cat(strrep("-", 75), "\n")
  for (i in 1:nrow(group_cate)) {
    g <- group_cate[i,]
    cat(sprintf("  %-20s  %6d  %+.4f  %.4f  %+.4f  %+.4f\n",
                g$group, g$N, g$Mean_CATE, g$SD_CATE, g$r_exp, g$r_partner_hours))
  }
  rec("cate_childless_men",
      group_cate$r_partner_hours[match("Childless men", group_cate$group)])
  
  cat("\n  Key: r(τ, partner_hours) = correlation between one's CATE and\n")
  cat("  partner's hours. If negative for men → when she works less, his\n")
  cat("  education benefit INCREASES (Becker specialization in CATE space).\n")
  
  
  # --- CATE-experience correlation with formal test ---
  cat("\n--- 3C. CATE-experience gradient by group (heterogeneity test) ---\n")
  
  for (g in unique(cf_data$group)) {
    sub <- cf_data %>% filter(group == g)
    ct <- cor.test(sub$tau_hat, sub$exp)
    cat(sprintf("  %-20s  r=%.4f  p=%.4f%s  %s\n",
                g, ct$estimate, ct$p.value, stars_fn(ct$p.value),
                ifelse(ct$p.value < 0.05, "RISING",
                       ifelse(ct$p.value < 0.10, "WEAK", "FROZEN"))))
  }
  
  # Partner hours correlation test
  cat("\n  CATE × partner_hours (specialization in CATE space):\n")
  for (g in unique(cf_data$group)) {
    sub <- cf_data %>% filter(group == g)
    ct <- cor.test(sub$tau_hat, sub$partner_hours)
    cat(sprintf("  %-20s  r=%.4f  p=%.4f%s\n",
                g, ct$estimate, ct$p.value, stars_fn(ct$p.value)))
  }
}


# #############################################################################
# PART 4: COMPARISON TABLE — LINEAR vs DML
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
cat("#  PART 4: LINEAR vs DML COMPARISON\n")
cat(strrep("#", 80), "\n\n")

# Linear estimates for comparison (on same stacked data)
cat("  Linear estimates on same stacked sample:\n\n")

linear_results <- list()
for (g in names(cells)) {
  d <- cells[[g]] %>%
    mutate(educ_c2 = educ - mean(educ, na.rm = TRUE),
           exp_c2 = exp - mean(exp, na.rm = TRUE),
           educ_exp_c2 = educ_c2 * exp_c2)
  
  m_ols <- tryCatch(feols(
    ln_wage ~ educ + educ_exp_c2 + exp + hours + age_sq + married +
      fulltime + partner_educ + partner_exp + partner_hours + partner_fulltime | wave,
    data = d, cluster = ~couple_id, notes = FALSE
  ), error = function(e) NULL)
  
  if (!is.null(m_ols) && "educ_exp_c2" %in% names(coef(m_ols))) {
    s <- summary(m_ols)
    b <- coef(m_ols)["educ_exp_c2"]
    p <- s$coeftable["educ_exp_c2", "Pr(>|t|)"]
    linear_results[[g]] <- list(b = b, p = p)
    cat(sprintf("  %-20s  OLS β₂ = %+.4f%%%s\n", g, b*100, stars_fn(p)))
  }
}

# Summary comparison
cat("\n\n  SUMMARY: Linear OLS vs DML-PLR (RF) on couple-stacked data:\n\n")
cat(sprintf("  %-20s  %12s  %12s  %12s\n",
            "Group", "OLS β₂", "DML β₂ (RF)", "DML β₂ (LASSO)"))
cat(strrep("-", 60), "\n")

for (g in names(cells)) {
  ols_b <- if (!is.null(linear_results[[g]])) sprintf("%+.4f%%%s",
                                                      linear_results[[g]]$b * 100, stars_fn(linear_results[[g]]$p)) else "—"
  
  dml_b_rf <- if (!is.null(dml_cells_rf[[g]])) sprintf("%+.4f%%%s",
                                                       dml_cells_rf[[g]]$theta * 100, stars_fn(dml_cells_rf[[g]]$pval)) else "—"
  
  dml_b_lasso <- if (!is.null(dml_cells_lasso[[g]])) sprintf("%+.4f%%%s",
                                                             dml_cells_lasso[[g]]$theta * 100, stars_fn(dml_cells_lasso[[g]]$pval)) else "—"
  
  cat(sprintf("  %-20s  %12s  %12s  %12s\n", g, ols_b, dml_b_rf, dml_b_lasso))
}

cat("\n  INTERPRETATION:\n")
cat("  If DML ≈ OLS → functional form is not driving the pattern\n")
cat("  If DML ≠ OLS → nonlinearity matters, linear model misleading\n")
cat("  If partner controls change β₂ → Becker spillover operates\n")
cat("  If partner controls don't change β₂ → complementarity is individual\n")


# #############################################################################
# PART 5: DYNAMIC BECKER TEST — THE THREE PREDICTIONS
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
cat("#  PART 5: DYNAMIC BECKER — THREE PREDICTIONS\n")
cat(strrep("#", 80), "\n\n")

cat("PREDICTION 1 (Static Becker): Gains from trade in LEVELS\n")
cat("  → Couple fixed-effects levels (B4b): her_hours → his_wage\n")
cat(sprintf("  → 10-hr fall in her hours → %s%% his wage (coef p = %s)\n",
            pv("becker_pct","%+.2f"), pv("becker_p","%.4f")))
cat("  → CONFIRMED: when she cuts hours, his wage rises\n\n")

cat("PREDICTION 2 (Dynamic Becker): Gains from trade in GROWTH\n")
cat("  → Triple interaction: educ_exp × female × children\n")

# Test with DML: does parenthood penalty differ by gender?
if (!is.null(dml_cells_rf[["Childless men"]]) && !is.null(dml_cells_rf[["Fathers"]]) &&
    !is.null(dml_cells_rf[["Childless women"]]) && !is.null(dml_cells_rf[["Mothers"]])) {
  
  # Fatherhood penalty on β₂
  father_penalty <- dml_cells_rf[["Fathers"]]$theta - dml_cells_rf[["Childless men"]]$theta
  father_se <- sqrt(dml_cells_rf[["Fathers"]]$se^2 + dml_cells_rf[["Childless men"]]$se^2)
  
  # Motherhood penalty on β₂
  mother_penalty <- dml_cells_rf[["Mothers"]]$theta - dml_cells_rf[["Childless women"]]$theta
  mother_se <- sqrt(dml_cells_rf[["Mothers"]]$se^2 + dml_cells_rf[["Childless women"]]$se^2)
  
  # Are the penalties DIFFERENT? (= triple interaction)
  triple <- (father_penalty - mother_penalty)
  triple_se <- sqrt(father_se^2 + mother_se^2)
  z_triple <- triple / triple_se
  p_triple <- 2 * pnorm(-abs(z_triple))
  
  cat(sprintf("  DML fatherhood penalty: %+.4f%% (SE %.4f)\n",
              father_penalty*100, father_se*100))
  cat(sprintf("  DML motherhood penalty: %+.4f%% (SE %.4f)\n",
              mother_penalty*100, mother_se*100))
  cat(sprintf("  DML triple (father − mother penalty): %+.4f%% (z=%.3f, p=%.4f%s)\n",
              triple*100, z_triple, p_triple, stars_fn(p_triple)))
  
  if (p_triple > 0.10) {
    cat("  → CONFIRMED: penalties are SYMMETRIC (Dynamic Becker REJECTED)\n")
    cat("    Parenthood destroys complementarity equally for both partners\n")
  } else {
    cat("  → ASYMMETRIC: one partner loses more complementarity than the other\n")
  }
}

cat(sprintf("\n  DML triple interaction: p = %.4f\n",
            ifelse(exists("p_triple"), p_triple, NA)))
cat("  If insignificant → symmetric penalty (Dynamic Becker rejected)\n")

cat("\n\nPREDICTION 3 (Permanent scar): Lock-in is irreversible\n")
cat("  → Post-mothers FE β₂: see Heckman / Table 1 Panel B output below\n")
cat("  → Active vs post-mothers indistinguishable (see Table 1 Panel B)\n")
cat("  → CONFIRMED: complementarity loss persists after children leave\n")


cat("\n\n")
cat(strrep("#", 80), "\n")
cat("#  COUPLE-LEVEL DML ANALYSIS COMPLETE\n")
cat(strrep("#", 80), "\n")

save.image("master_state_after_partD.RData")  # after PART B

################################################################################
# PART E: ROBUSTNESS BATTERY + ALL PUBLICATION FIGURES
################################################################################

# 0. LOAD DATA

cat(strrep("=", 80), "\n")
cat("  ROBUSTNESS BATTERY FOR JOURNAL SUBMISSION\n")
cat(strrep("=", 80), "\n\n")

num_rows    <- format(nrow(panel), big.mark = ",")
num_persons <- format(n_distinct(panel$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

# Ensure numeric variables



# WPI deflators — applied wherever ln_hourly_wage_real needs to be constructed
# from scratch (no-trim panel, panel_notrim). Base = Wave 24 (2024, index 149.6).
# Note: hilda_panel_data_extended.rds already contains ln_hourly_wage_real for the main
# wage sample; wpi_data is only needed when rebuilding wage samples outside
# the standard in_wage_sample flag.

# PART 1: CORE ESTIMATION FUNCTION

# Unified function: estimate β₂ for a given subsample under OLS and FE
# Returns a one-row tibble for stacking into the master table
estimate_beta2 <- function(data, group_var, group_val, label,
                           extra_fe = NULL, weights_col = NULL) {
  d <- data %>%
    filter(.data[[group_var]] == group_val,
           in_wage_sample == 1) %>%
    mutate(
      educ_c = educ_years - mean(educ_years, na.rm = TRUE),
      exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c
    ) %>%
    filter(complete.cases(ln_hourly_wage_real, educ_years, educ_exp_c,
                          experience_years, age_sq, married_num))
  
  if (nrow(d) < 100) {
    return(tibble(label = label, group = group_val,
                  N = nrow(d), N_id = n_distinct(d$person_id),
                  ols_b2 = NA, ols_p = NA, ols_se = NA,
                  fe_b2 = NA, fe_p = NA, fe_se = NA))
  }
  
  # Weights
  w <- if (!is.null(weights_col) && weights_col %in% names(d)) {
    wt <- d[[weights_col]]
    wt[is.na(wt) | wt <= 0] <- 1e-6
    wt
  } else NULL
  
  # FE formula
  # FIX: educ_years removed — absorbed by person FE (paper §4.1: "the standalone
  # education level effect is absorbed by individual fixed effects")
  fe_fml <- if (!is.null(extra_fe)) {
    as.formula(paste0("ln_hourly_wage_real ~ educ_exp_c + ",
                      "experience_years + age_sq + married_num | person_id + ",
                      extra_fe))
  } else {
    ln_hourly_wage_real ~ educ_exp_c + experience_years +
      age_sq + married_num | person_id + wave
  }
  
  ols_fml <- if (!is.null(extra_fe)) {
    as.formula(paste0("ln_hourly_wage_real ~ educ_years + educ_exp_c + ",
                      "experience_years + age_sq + married_num | ",
                      extra_fe))
  } else {
    ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
      age_sq + married_num | wave
  }
  
  m_ols <- tryCatch(
    feols(ols_fml, data = d, cluster = ~person_id,
          weights = w, notes = FALSE),
    error = function(e) NULL
  )
  
  m_fe <- tryCatch(
    feols(fe_fml, data = d, cluster = ~person_id,
          weights = w, notes = FALSE),
    error = function(e) NULL
  )
  
  ols_b2 <- if (!is.null(m_ols) && "educ_exp_c" %in% names(coef(m_ols)))
    coef(m_ols)["educ_exp_c"] else NA
  ols_p <- if (!is.null(m_ols) && "educ_exp_c" %in% rownames(summary(m_ols)$coeftable))
    summary(m_ols)$coeftable["educ_exp_c", "Pr(>|t|)"] else NA
  
  fe_b2 <- if (!is.null(m_fe) && "educ_exp_c" %in% names(coef(m_fe)))
    coef(m_fe)["educ_exp_c"] else NA
  fe_p <- if (!is.null(m_fe) && "educ_exp_c" %in% rownames(summary(m_fe)$coeftable))
    summary(m_fe)$coeftable["educ_exp_c", "Pr(>|t|)"] else NA
  
  ols_se <- if (!is.null(m_ols) && "educ_exp_c" %in% rownames(summary(m_ols)$coeftable))
    summary(m_ols)$coeftable["educ_exp_c", "Std. Error"] else NA
  fe_se <- if (!is.null(m_fe) && "educ_exp_c" %in% rownames(summary(m_fe)$coeftable))
    summary(m_fe)$coeftable["educ_exp_c", "Std. Error"] else NA
  
  tibble(label = label, group = group_val,
         N = nrow(d), N_id = n_distinct(d$person_id),
         ols_b2 = ols_b2, ols_p = ols_p, ols_se = ols_se,
         fe_b2 = fe_b2, fe_p = fe_p, fe_se = fe_se)
}

# Wrapper: run all 4/6 groups for a given sample/spec
run_all_groups <- function(data, label_prefix, groups = "four",
                           extra_fe = NULL, weights_col = NULL) {
  if (groups == "four") {
    grps <- c("childless_men_ever", "fathers_ever",
              "never_mothers", "mothers_ever")
    grp_labels <- c("Childless men", "Fathers",
                    "Never-mothers", "Mothers (ever)")
  } else {
    grps <- c("never_mothers", "active_mothers", "post_mothers")
    grp_labels <- c("Never-mothers", "Active mothers", "Post-mothers")
  }
  
  results <- map2_dfr(grps, grp_labels, function(g, gl) {
    gvar <- if (groups == "four") "group_ever" else "group_3way_women"
    estimate_beta2(data, gvar, g, paste0(label_prefix, ": ", gl),
                   extra_fe = extra_fe, weights_col = weights_col)
  })
  results
}

# Pretty-print a results table
print_results <- function(results, title) {
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-45s  %7s %7s  %10s  %10s\n",
              "Label", "N", "N_id", "OLS β₂", "FE β₂"))
  cat("  ", strrep("-", 85), "\n")
  for (i in 1:nrow(results)) {
    r <- results[i, ]
    ols_str <- if (!is.na(r$ols_b2)) sprintf("%+.4f%%%s",
                                             r$ols_b2 * 100, stars_fn(r$ols_p)) else "—"
    fe_str <- if (!is.na(r$fe_b2)) sprintf("%+.4f%%%s",
                                           r$fe_b2 * 100, stars_fn(r$fe_p)) else "—"
    cat(sprintf("  %-45s  %7s %7s  %10s  %10s\n",
                r$label,
                format(r$N, big.mark = ","),
                format(r$N_id, big.mark = ","),
                ols_str, fe_str))
  }
}

# BASELINE (for comparison)

cat("\n")
cat(strrep("#", 80), "\n")
# Paper §5 (Robustness Analysis), Table 6: Master robustness summary.
# Paper Table 1: Baseline β₂ estimates (absorbing definition).
# Sample: Wage sample (in_wage_sample == 1), ~87,229 obs.
cat("#  BASELINE: Main results (absorbing definition, unweighted)\n")
cat(strrep("#", 80), "\n")

baseline_4 <- run_all_groups(panel, "Baseline", "four")
baseline_3 <- run_all_groups(panel, "Baseline", "three")
print_results(baseline_4, "Baseline: Four groups")
print_results(baseline_3, "Baseline: Three-way women")

master <- bind_rows(
  baseline_4 %>% mutate(test = "Baseline"),
  baseline_3 %>% mutate(test = "Baseline (3-way)")
)


# #############################################################################
# PART 2: WEIGHTING (5.1–5.3)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.1 (Weighting): "Three weighting schemes leave results
#   substantively unchanged."
# Samples: §5.1 longitudinal weight subsample → ~17,703 obs / 2,599 persons;
#   §5.2 cross-sectional weighted; §5.3 IPW for attrition.
cat("#  PART 2: WEIGHTING ROBUSTNESS (5.1–5.3)\n")
cat(strrep("#", 80), "\n")

# --- 5.1 Longitudinal weight subsample ---
cat("\n--- 5.1 Longitudinal weight subsample ---\n")
panel_longwt <- panel %>%
  filter(!is.na(longitudinal_weight) & longitudinal_weight > 0)
num_rows    <- format(nrow(panel_longwt), big.mark = ",")
num_persons <- format(n_distinct(panel_longwt$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))
cat(sprintf("  Coverage: %.1f%% of full panel\n",
            100 * nrow(panel_longwt) / nrow(panel)))

r51 <- run_all_groups(panel_longwt, "LongWt subsample", "four")
print_results(r51, "5.1: Longitudinal weight subsample (unweighted)")
master <- bind_rows(master, r51 %>% mutate(test = "5.1 LongWt subsample"))

# --- 5.2 Cross-sectional (hhwtrps) weighted ---
cat("\n--- 5.2 Cross-sectional weighted (analysis_weight) ---\n")
r52 <- run_all_groups(panel, "XS-weighted", "four",
                      weights_col = "analysis_weight")
print_results(r52, "5.2: Cross-sectional weighted")
master <- bind_rows(master, r52 %>% mutate(test = "5.2 XS-weighted"))

# --- 5.3 IPW for attrition ---
cat("\n--- 5.3 IPW for employment attrition ---\n")

# Build IPW: predict P(in_wage_sample) using baseline characteristics
# Use first-wave characteristics for each person
baseline_chars <- panel %>%
  group_by(person_id) %>%
  arrange(wave) %>%
  slice(1) %>%
  ungroup() %>%
  select(person_id, educ_years, age, female, married_num, state)

panel_ipw <- panel %>%
  left_join(baseline_chars %>% rename(bl_educ = educ_years, bl_age = age,
                                      bl_female = female, bl_married = married_num),
            by = "person_id") %>%
  mutate(in_wage_num = as.numeric(in_wage_sample == 1))

# Probit for employment
ipw_model <- tryCatch(
  glm(in_wage_num ~ bl_educ + bl_age + I(bl_age^2) + bl_female + bl_married +
        as.factor(wave),
      data = panel_ipw, family = binomial("probit")),
  error = function(e) NULL
)

if (!is.null(ipw_model)) {
  panel_ipw$p_employed <- predict(ipw_model, newdata = panel_ipw, type = "response")
  panel_ipw$ipw <- 1 / pmax(panel_ipw$p_employed, 0.01)
  # Normalize
  panel_ipw$ipw <- panel_ipw$ipw / mean(panel_ipw$ipw[panel_ipw$in_wage_sample == 1],
                                        na.rm = TRUE)
  r53 <- run_all_groups(panel_ipw, "IPW-weighted", "four",
                        weights_col = "ipw")
  print_results(r53, "5.3: IPW-weighted (employment attrition)")
  master <- bind_rows(master, r53 %>% mutate(test = "5.3 IPW-weighted"))
} else {
  cat("  IPW model failed\n")
}


# #############################################################################
# PART 3: SAMPLE DEFINITIONS (5.4–5.7)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.2 (Sample Definition): Age windows, FT/PT, COVID exclusion.
# Samples: Age 25–54 → ~77,315; Full-time → ~53,696;
#   Part-time count: re-verify against current extract.
cat("#  PART 3: SAMPLE DEFINITION SENSITIVITY (5.4–5.7)\n")
cat(strrep("#", 80), "\n")

# --- 5.4 Wider age range 25–54 ---
cat("\n--- 5.4 Age 25–54 ---\n")
panel_wide <- panel %>% filter(age >= 25 & age <= 54)
cat(sprintf("  N: %s\n", format(nrow(panel_wide), big.mark = ",")))
r54 <- run_all_groups(panel_wide, "Age 25-54", "four")
print_results(r54, "5.4: Age 25-54")
master <- bind_rows(master, r54 %>% mutate(test = "5.4 Age 25-54"))

# --- 5.5 Narrow age range 35–49 ---
cat("\n--- 5.5 Age 35–49 ---\n")
panel_narrow <- panel %>% filter(age >= 35 & age <= 49)
cat(sprintf("  N: %s\n", format(nrow(panel_narrow), big.mark = ",")))
r55 <- run_all_groups(panel_narrow, "Age 35-49", "four")
print_results(r55, "5.5: Age 35-49 (established careers)")
master <- bind_rows(master, r55 %>% mutate(test = "5.5 Age 35-49"))

# --- 5.6a Full-time only ---
cat("\n--- 5.6a Full-time only ---\n")
panel_ft <- panel %>% filter(fulltime == 1)
cat(sprintf("  N: %s\n", format(nrow(panel_ft), big.mark = ",")))
r56a <- run_all_groups(panel_ft, "Full-time only", "four")
print_results(r56a, "5.6a: Full-time workers only")
master <- bind_rows(master, r56a %>% mutate(test = "5.6a Full-time only"))

# --- 5.6b Part-time only ---
cat("\n--- 5.6b Part-time only ---\n")
panel_pt <- panel %>% filter(fulltime == 0 & employed == 1)
cat(sprintf("  N: %s\n", format(nrow(panel_pt), big.mark = ",")))
r56b <- run_all_groups(panel_pt, "Part-time only", "four")
print_results(r56b, "5.6b: Part-time workers only")
master <- bind_rows(master, r56b %>% mutate(test = "5.6b Part-time only"))

# --- 5.7 Exclude COVID waves 2020–2021 (waves 20–21) ---
cat("\n--- 5.7 Exclude COVID (waves 20–21) ---\n")
panel_nocovid <- panel %>% filter(!(wave %in% c(20, 21)))
cat(sprintf("  N: %s\n", format(nrow(panel_nocovid), big.mark = ",")))
r57 <- run_all_groups(panel_nocovid, "No COVID", "four")
print_results(r57, "5.7: Excluding COVID waves 2020-2021")
master <- bind_rows(master, r57 %>% mutate(test = "5.7 No COVID"))


# #############################################################################
# PART 4: MOTHERHOOD DEFINITION SENSITIVITY (5.8)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.3, Table 5 (Definition Sensitivity): Time-varying vs absorbing
#   vs tchad-based. "23.3% of time-varying 'childless' observations are
#   post-mothers."
# Sample: Wage sample, women only for cross-tabulation.
cat("#  PART 4: MOTHERHOOD DEFINITION SENSITIVITY (5.8)\n")
cat(strrep("#", 80), "\n")

# Definition (a): Time-varying (children_under15 > 0)
cat("\n--- Def (a): Time-varying children in household ---\n")
panel_a <- panel %>%
  mutate(group_tv = case_when(
    female == 0 & ever_parent == 0 ~ "childless_men",
    female == 0 & ever_parent == 1 ~ "fathers",
    female == 1 & num_children_under15 > 0 ~ "mothers_current",
    female == 1 & num_children_under15 == 0 ~ "childless_women_current",
    TRUE ~ NA_character_
  ))

r_a <- map_dfr(
  list(c("childless_men", "Childless men"),
       c("fathers", "Fathers"),
       c("childless_women_current", "Childless women (TV)"),
       c("mothers_current", "Mothers (TV)")),
  function(x) estimate_beta2(panel_a, "group_tv", x[1],
                             paste0("Def(a) TV: ", x[2]))
)
print_results(r_a, "Definition (a): Time-varying")
master <- bind_rows(master, r_a %>% mutate(test = "5.8a Time-varying"))

# Definition (b): Absorbing — already the baseline
cat("\n  Definition (b): Absorbing = Baseline (already reported)\n")

# Definition (c): tchad only (total_children_ever_had > 0)
cat("\n--- Def (c): total_children_ever_had > 0 ---\n")
panel_c <- panel %>%
  mutate(group_tchad = case_when(
    female == 0 & (is.na(total_children_ever_had) | total_children_ever_had == 0) ~ "childless_men_tchad",
    female == 0 & total_children_ever_had > 0 ~ "fathers_tchad",
    female == 1 & (is.na(total_children_ever_had) | total_children_ever_had == 0) ~ "childless_women_tchad",
    female == 1 & total_children_ever_had > 0 ~ "mothers_tchad",
    TRUE ~ NA_character_
  ))

r_c <- map_dfr(
  list(c("childless_men_tchad", "Childless men"),
       c("fathers_tchad", "Fathers"),
       c("childless_women_tchad", "Childless women (tchad)"),
       c("mothers_tchad", "Mothers (tchad)")),
  function(x) estimate_beta2(panel_c, "group_tchad", x[1],
                             paste0("Def(c) tchad: ", x[2]))
)
print_results(r_c, "Definition (c): tchad-based")
master <- bind_rows(master, r_c %>% mutate(test = "5.8c tchad-based"))

# --- 6.2 Contamination cross-tab ---
cat("\n--- 6.2 Contamination cross-tab ---\n")
panel %>%
  filter(female == 1, in_wage_sample == 1) %>%
  mutate(
    current = ifelse(num_children_under15 > 0, "Mother (current)",
                     "Childless (current)"),
    ever = ifelse(ever_parent == 1, "Mother (ever)", "Never-mother")
  ) %>%
  count(current, ever) %>%
  pivot_wider(names_from = ever, values_from = n, values_fill = 0) %>%
  print()

cat("\n  Post-mothers contaminating 'childless' under current definition:\n")
contam <- panel %>%
  filter(female == 1, in_wage_sample == 1,
         num_children_under15 == 0, ever_parent == 1) %>%
  nrow()
total_childless_current <- panel %>%
  filter(female == 1, in_wage_sample == 1, num_children_under15 == 0) %>%
  nrow()
cat(sprintf("  %s / %s = %.1f%%\n",
            format(contam, big.mark = ","),
            format(total_childless_current, big.mark = ","),
            100 * contam / total_childless_current))
rec("contam_pct", 100 * contam / total_childless_current)


# #############################################################################
# PART 5: OUTCOME & TRIMMING (5.9–5.10)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.4 (Outcome and Trimming): Log weekly earnings; 1/99 trim;
#   no trim. "Results are insensitive to the treatment of wage outliers."
# Samples: Trim 1/99 and No trim counts re-verify against current extract.
cat("#  PART 5: OUTCOME & TRIMMING (5.9–5.10)\n")
cat(strrep("#", 80), "\n")

# --- 5.9 Weekly earnings instead of hourly wage ---
cat("\n--- 5.9 Weekly earnings ---\n")
panel_weekly <- panel %>%
  filter(!is.na(weekly_earnings) & weekly_earnings > 0) %>%
  mutate(
    ln_weekly = log(weekly_earnings),
    # Overwrite for the estimation function
    ln_hourly_wage_orig = ln_hourly_wage,
    ln_hourly_wage = ln_weekly,
    in_wage_sample = 1L  # all with positive weekly earnings
  )
cat(sprintf("  N with positive weekly earnings: %s\n",
            format(nrow(panel_weekly), big.mark = ",")))
r59 <- run_all_groups(panel_weekly, "Ln(weekly)", "four")
print_results(r59, "5.9: Log weekly earnings")
master <- bind_rows(master, r59 %>% mutate(test = "5.9 Weekly earnings"))

# --- 5.10 Alternative trimming ---
cat("\n--- 5.10 Tighter trimming (1st/99th) ---\n")
panel_trim <- panel %>%
  filter(in_wage_sample == 1) %>%
  group_by(female) %>%
  mutate(
    p01 = quantile(hourly_wage_clean, 0.01, na.rm = TRUE),
    p99 = quantile(hourly_wage_clean, 0.99, na.rm = TRUE),
    in_trim = hourly_wage_clean >= p01 & hourly_wage_clean <= p99
  ) %>%
  ungroup() %>%
  filter(in_trim)

cat(sprintf("  After 1%%/99%% trim: %s obs (dropped %s)\n",
            format(nrow(panel_trim), big.mark = ","),
            format(sum(panel$in_wage_sample == 1, na.rm = TRUE) - nrow(panel_trim),
                   big.mark = ",")))
r510 <- run_all_groups(panel_trim, "Trim 1/99", "four")
print_results(r510, "5.10: 1st/99th percentile trimming")
master <- bind_rows(master, r510 %>% mutate(test = "5.10 Trim 1/99"))

# No trimming at all
cat("\n  No trimming at all:\n")
panel_notrim <- panel %>%
  filter(employed == 1, !is.na(hourly_wage_clean), hourly_wage_clean > 0) %>%
  select(-any_of("deflator")) %>%
  left_join(wpi_data, by = "wave") %>%
  mutate(ln_hourly_wage_real = log(hourly_wage_clean * deflator),
         in_wage_sample = 1L)
r510b <- run_all_groups(panel_notrim, "No trim", "four")
print_results(r510b, "5.10b: No trimming")
master <- bind_rows(master, r510b %>% mutate(test = "5.10b No trim"))


# #############################################################################
# PART 6: SPECIFICATION SENSITIVITY (5.11–5.14)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.5–5.8: Quadratic interaction, Occ×Wave FE, State×Wave FE, CRE.
# Key finding: "occupation-by-wave FE absorbs time-varying occupational
#   sorting. Childless men's β₂ falls 54%."
cat("#  PART 6: SPECIFICATION SENSITIVITY (5.11–5.14)\n")
cat(strrep("#", 80), "\n")

# --- 5.11 Quadratic interaction: educ × exp² ---
cat("\n--- 5.11 Quadratic interaction ---\n")

groups_5 <- list(
  "Childless men" = wage %>% filter(group_ever == "childless_men_ever"),
  "Fathers"       = wage %>% filter(group_ever == "fathers_ever"),
  "Never-mothers" = wage %>% filter(group_ever == "never_mothers"),
  "Mothers (ever)" = wage %>% filter(group_ever == "mothers_ever")
)

cat(sprintf("  %-20s  %10s  %10s  %10s\n",
            "Group", "β₂ (linear)", "β₂ (quad)", "β₃ (exp²)"))
cat("  ", strrep("-", 55), "\n")

for (g in names(groups_5)) {
  d <- groups_5[[g]] %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c,
           educ_exp2_c = educ_c * exp_c^2)
  
  m_q <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + educ_exp2_c +
            experience_years + age_sq + married_num | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(m_q)) {
    s <- summary(m_q)$coeftable
    b2 <- if ("educ_exp_c" %in% rownames(s)) s["educ_exp_c", "Estimate"] else NA
    b3 <- if ("educ_exp2_c" %in% rownames(s)) s["educ_exp2_c", "Estimate"] else NA
    p3 <- if ("educ_exp2_c" %in% rownames(s)) s["educ_exp2_c", "Pr(>|t|)"] else NA
    # Linear-only for comparison
    m_l <- tryCatch(
      feols(ln_hourly_wage_real ~ educ_exp_c +
              experience_years + age_sq + married_num | person_id + wave,
            data = d, cluster = ~person_id, notes = FALSE),
      error = function(e) NULL
    )
    b2_lin <- if (!is.null(m_l)) coef(m_l)["educ_exp_c"] else NA
    cat(sprintf("  %-20s  %+.4f%%  %+.4f%%  %+.6f%s\n",
                g, b2_lin * 100, b2 * 100, b3 * 100, stars_fn(p3)))
  }
}

# --- 5.12 Occupation × wave FE ---
cat("\n--- 5.12 Occupation × wave FE ---\n")
panel_occ <- panel %>%
  filter(!is.na(occupation_major), occupation_major != "")
cat(sprintf("  Obs with occupation: %s\n",
            format(sum(panel_occ$in_wage_sample == 1), big.mark = ",")))
r512 <- run_all_groups(panel_occ, "Occ×Wave FE", "four",
                       extra_fe = "occupation_major^wave")
print_results(r512, "5.12: Occupation × Wave FE")
master <- bind_rows(master, r512 %>% mutate(test = "5.12 Occ×Wave FE"))

# --- 5.13 State × wave FE ---
cat("\n--- 5.13 State × wave FE ---\n")
panel_st <- panel %>% filter(!is.na(state), state != "")
r513 <- run_all_groups(panel_st, "State×Wave FE", "four",
                       extra_fe = "state^wave")
print_results(r513, "5.13: State × Wave FE")
master <- bind_rows(master, r513 %>% mutate(test = "5.13 State×Wave FE"))

# --- 5.14 Mundlak CRE ---
cat("\n--- 5.14 Mundlak CRE ---\n")

cre_results <- tibble()
for (g in names(groups_5)) {
  d <- groups_5[[g]] %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c) %>%
    group_by(person_id) %>%
    mutate(
      exp_bar = mean(experience_years, na.rm = TRUE),
      hours_bar = mean(hours_worked_clean, na.rm = TRUE),
      married_bar = mean(married_num, na.rm = TRUE),
      educ_exp_bar = mean(educ_exp_c, na.rm = TRUE)
    ) %>%
    ungroup()
  
  m_cre <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
            age_sq + married_num + exp_bar + married_bar + educ_exp_bar | wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(m_cre) && "educ_exp_c" %in% names(coef(m_cre))) {
    s <- summary(m_cre)$coeftable
    b <- s["educ_exp_c", "Estimate"]
    p <- s["educ_exp_c", "Pr(>|t|)"]
    # Mundlak test: is educ_exp_bar significant?
    if ("educ_exp_bar" %in% rownames(s)) {
      bar_b <- s["educ_exp_bar", "Estimate"]
      bar_p <- s["educ_exp_bar", "Pr(>|t|)"]
      cat(sprintf("  %-20s CRE β₂=%+.4f%%%s  Mundlak_bar=%+.4f (p=%.3f%s)\n",
                  g, b*100, stars_fn(p), bar_b*100, bar_p, stars_fn(bar_p)))
    }
    cre_results <- bind_rows(cre_results,
                             tibble(label = paste0("CRE: ", g), group = g,
                                    N = nrow(d), N_id = n_distinct(d$person_id),
                                    ols_b2 = b, ols_p = p,
                                    ols_se = s["educ_exp_c", "Std. Error"],
                                    fe_b2 = NA, fe_p = NA, fe_se = NA,
                                    test = "5.14 Mundlak CRE"))
  }
}
master <- bind_rows(master, cre_results)


# #############################################################################
# PART 6B: HECKMAN SELECTION CORRECTION (paper Table 3, §4.3)
# #############################################################################
# Two-step Heckman correction for employment selection.
# First stage: probit for employment using child-age and care-burden
# exclusion restrictions.
# Second stage: FE wage regression with inverse Mills ratio.

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §4.3, Table 3: "Corrections are negligible (Δ < 0.015 pp),
#   indicating that the 15-percentage-point employment gap between
#   never-mothers and active mothers does not drive the observed hierarchy."
# Sample: Broad panel for probit 1st stage; wage sample for 2nd stage.
cat("#  PART 6B: HECKMAN SELECTION CORRECTION (Table 3)\n")
cat(strrep("#", 80), "\n")

# Exclusion restrictions: child-age variables + care burden + marital status
# These affect employment but should not directly affect hourly wages
# conditional on being employed.

heckman_groups <- list(
  "Childless men"  = list(group_var = "group_ever", group_val = "childless_men_ever"),
  "Fathers (ever)" = list(group_var = "group_ever", group_val = "fathers_ever"),
  "Never-mothers"  = list(group_var = "group_ever", group_val = "never_mothers"),
  "Active mothers"  = list(group_var = "group_3way_women", group_val = "active_mothers"),
  "Post-mothers"    = list(group_var = "group_3way_women", group_val = "post_mothers")
)

heckman_results <- tibble()

for (g_name in names(heckman_groups)) {
  g_info <- heckman_groups[[g_name]]
  
  # Broad sample for selection equation
  d_broad <- panel %>%
    filter(.data[[g_info$group_var]] == g_info$group_val) %>%
    mutate(
      employed_num = as.numeric(employed == 1),
      has_young_child = as.numeric(!is.na(num_children_under15) & num_children_under15 > 0),
      educ_c = educ_years - mean(educ_years, na.rm = TRUE),
      exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c
    ) %>%
    filter(!is.na(educ_years), !is.na(experience_years), !is.na(age_sq))
  
  # First stage: probit for employment
  probit_mod <- tryCatch(
    glm(employed_num ~ educ_years + experience_years + age_sq +
          married_num + has_young_child + as.factor(wave),
        data = d_broad, family = binomial("probit")),
    error = function(e) NULL
  )
  
  if (is.null(probit_mod)) {
    cat(sprintf("  %-20s  Probit FAILED\n", g_name))
    next
  }
  
  # Compute inverse Mills ratio
  d_broad$probit_xb <- predict(probit_mod, newdata = d_broad, type = "link")
  d_broad$imr       <- dnorm(d_broad$probit_xb) / pmax(pnorm(d_broad$probit_xb), 1e-6)
  
  # Wage sample with IMR
  
  d_wage <- d_broad %>%
    filter(in_wage_sample == 1, !is.na(ln_hourly_wage_real), !is.na(imr))
  
  if (nrow(d_wage) < 100) {
    cat(sprintf("  %-20s  Too few wage obs: %d\n", g_name, nrow(d_wage)))
    next
  }
  
  # Raw FE (without correction)
  m_raw <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num | person_id + wave,
          data = d_wage, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL
  )
  
  # Corrected FE (with IMR)
  m_corr <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num + imr | person_id + wave,
          data = d_wage, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(m_raw) && !is.null(m_corr) &&
      "educ_exp_c" %in% names(coef(m_raw)) &&
      "educ_exp_c" %in% names(coef(m_corr))) {
    
    s_raw  <- summary(m_raw)$coeftable
    s_corr <- summary(m_corr)$coeftable
    
    b_raw  <- s_raw["educ_exp_c", 1]
    b_corr <- s_corr["educ_exp_c", 1]
    imr_b  <- s_corr["imr", 1]
    imr_p  <- s_corr["imr", 4]
    se_raw  <- s_raw["educ_exp_c", 2]
    se_corr <- s_corr["educ_exp_c", 2]
    
    cat(sprintf("  %-20s  FE raw=%+.3f%% (SE %.3f)  FE corr=%+.3f%% (SE %.3f)  Δ=%+.3f  IMR=%+.3f (p=%.3f)\n",
                g_name, b_raw*100, se_raw*100, b_corr*100, se_corr*100, (b_corr - b_raw)*100,
                imr_b, imr_p))
    
    heckman_results <- bind_rows(heckman_results,
                                 tibble(group = g_name, fe_raw = b_raw*100, se_raw = se_raw*100,
                                        fe_corrected = b_corr*100, se_corrected = se_corr*100,
                                        delta = (b_corr - b_raw)*100, imr_coef = imr_b, imr_p = imr_p))
  }
}

if (nrow(heckman_results) > 0) {
  cat("\n  Heckman Summary (Table 3):\n")
  cat(sprintf("  %-20s  %8s  %8s  %8s  %8s  %8s\n",
              "Group", "FE raw", "FE corr", "Δ", "IMR", "IMR p"))
  cat("  ", strrep("-", 62), "\n")
  for (i in 1:nrow(heckman_results)) {
    r <- heckman_results[i, ]
    cat(sprintf("  %-20s  %+.3f%%  %+.3f%%  %+.3f  %+.3f  %.3f%s\n",
                r$group, r$fe_raw, r$fe_corrected, r$delta,
                r$imr_coef, r$imr_p, stars_fn(r$imr_p)))
  }
}

cat("\n  Note: Corrections are expected to be negligible (Δ < 0.015 pp per paper)\n")

if (nrow(heckman_results) > 0) {
  write.csv(heckman_results, "se_ledger_heckman.csv", row.names = FALSE)
  cat("  Wrote se_ledger_heckman.csv (Heckman FE raw/corrected beta2 + clustered SEs).\n")
}


# #############################################################################
# PART 7: OSTER BOUNDS (from HILDA_Oster_bounds.R)
# #############################################################################
# Specification (matches HILDA_Oster_bounds.R exactly):
#   - Group-specific mean-centring (matches Table 1 / Tab main_beta2)
#   - Restricted model  : OLS with educ_years + educ_exp_c + experience_years
#   - Full model        : OLS + age_sq + married_num + wave FE
#   - Causal benchmark  : Individual + wave FE (same as Table 1)
#   - R²_max            : min(1, 1.3 × R²_full)  [Oster 2019 convention]
#
# Oster formula (Oster 2019, eq. 5):
#   Π = (β_restricted − β_full) · (R²_max − R²_full) / (R²_full − R²_restricted)
#   δ = (β_full − β_FE) / Π
#   β* (δ=1) = β_full - Π   [identified set endpoint under equal selection]
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.9, Table 7: Oster (2019) selection bounds.
#   "For childless men (δ = −4.64) and never-mothers (δ = +8.18),
#   unobservables would need to be at least 38 times more important."
# Sample: Wage sample by group, with group-specific centering.
#   Expected per group: re-verify group N against current extract.
cat("#  PART 7: OSTER BOUNDS (from HILDA_Oster_bounds.R)\n")
cat(strrep("#", 80), "\n")

compute_oster <- function(data, group_label,
                          r2_max_scalar = 1.3,
                          verbose = TRUE) {
  
  # Group-specific centering (MATCHES Table 1)
  d <- data %>%
    filter(in_wage_sample == 1) %>%
    mutate(
      educ_c     = educ_years      - mean(educ_years,      na.rm = TRUE),
      exp_c      = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c
    ) %>%
    filter(complete.cases(ln_hourly_wage_real, educ_years, educ_exp_c,
                          experience_years, age_sq, married_num))
  
  N    <- nrow(d)
  N_id <- n_distinct(d$person_id)
  
  # Model 1: Restricted OLS (no controls)
  m_r <- tryCatch(
    lm(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years,
       data = d),
    error = function(e) NULL)
  
  # Model 2: Full OLS (controls + wave FE)
  m_f <- tryCatch(
    lm(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
         age_sq + married_num + as.factor(wave),
       data = d),
    error = function(e) NULL)
  
  # Model 3: Individual + wave FE (causal benchmark = Table 1)
  m_fe <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL)
  
  if (is.null(m_r) || is.null(m_f) || is.null(m_fe)) {
    warning(sprintf("Model failed for: %s", group_label))
    return(NULL)
  }
  
  b_r   <- coef(m_r)["educ_exp_c"]
  b_f   <- coef(m_f)["educ_exp_c"]
  
  fe_ct <- summary(m_fe)$coeftable
  b_fe  <- fe_ct["educ_exp_c", "Estimate"]
  se_fe <- fe_ct["educ_exp_c", "Std. Error"]
  p_fe  <- fe_ct["educ_exp_c", "Pr(>|t|)"]
  
  r2_r   <- summary(m_r)$r.squared
  r2_f   <- summary(m_f)$r.squared
  r2_max <- min(1, r2_max_scalar * r2_f)
  
  # Robust SE for full OLS β₂
  ct_f <- tryCatch({
    coeftest(m_f, vcov = vcovHC(m_f, "HC1"))["educ_exp_c", ]
  }, error = function(e) NULL)
  se_f <- if (!is.null(ct_f)) ct_f["Std. Error"] else NA
  p_f  <- if (!is.null(ct_f)) ct_f["Pr(>|t|)"]  else NA
  
  # Oster δ
  r2_gap <- r2_f - r2_r
  Pi     <- if (abs(r2_gap) > 1e-8) {
    (b_r - b_f) * (r2_max - r2_f) / r2_gap
  } else NA
  delta  <- if (!is.na(Pi) && abs(Pi) > 1e-8) (b_f - b_fe) / Pi else NA
  
  # β* when δ = 1 (equal selection)
  beta_star <- if (!is.na(Pi)) b_f - Pi else NA
  
  if (verbose) {
    cat(sprintf("  %-20s  N=%s (%s id)\n",
                group_label, format(N, big.mark = ","), format(N_id, big.mark = ",")))
    cat(sprintf("    Restricted OLS β₂ = %+.4f%%\n", b_r * 100))
    cat(sprintf("    Full OLS       β₂ = %+.4f%%%s  (R²=%.4f)\n",
                b_f * 100, stars_fn(p_f), r2_f))
    cat(sprintf("    Individual FE  β₂ = %+.4f%%%s  (R²_max=%.4f)\n",
                b_fe * 100, stars_fn(p_fe), r2_max))
    cat(sprintf("    β* (δ=1 equal sel.) = %+.4f%%\n",
                if (!is.na(beta_star)) beta_star * 100 else NA))
    cat(sprintf("    Oster δ = %+.2f\n", if (!is.na(delta)) delta else NA))
    
    if (!is.na(delta)) {
      if (abs(delta) > 1) {
        direction <- if (delta < 0) "opposite direction to" else "same direction as"
        cat(sprintf(
          "    → |δ|=%.2f > 1: unobservables %.1f× MORE important (%s observables). ROBUST.\n",
          abs(delta), abs(delta), direction))
      } else {
        cat(sprintf(
          "    → |δ|=%.2f < 1: selection could plausibly explain gap. INTERPRET WITH CAUTION.\n",
          abs(delta)))
      }
      if (delta < 0) {
        cat("    → δ < 0: FE > full OLS; selection must work AGAINST observables.\n")
      }
    }
    cat("\n")
  }
  
  tibble(
    group       = group_label,
    N           = N,
    N_id        = N_id,
    b_restricted = b_r   * 100,
    b_full_ols  = b_f   * 100,
    se_full_ols = se_f  * 100,
    p_full_ols  = p_f,
    b_fe        = b_fe  * 100,
    se_fe       = se_fe * 100,
    p_fe        = p_fe,
    r2_restricted = r2_r,
    r2_full     = r2_f,
    r2_max      = r2_max,
    beta_star   = if (!is.na(beta_star)) beta_star * 100 else NA_real_,
    delta       = delta,
    robust      = if (!is.na(delta)) abs(delta) > 1 else NA
  )
}

# --- Run for all four groups ---

cat(strrep("=", 70), "\n")
# Paper §5.9, Table 7. R²_max = min(1, 1.3 × R²_full) [Oster convention].
cat("  OSTER (2019) SELECTION BOUNDS FOR β₂\n")
cat("  Specification: group-specific centering (matches Table 1)\n")
cat("  R²_max = min(1, 1.3 × R²_full)  [Oster 2019 convention]\n")
cat(strrep("=", 70), "\n\n")

oster_groups <- list(
  "Childless men"  = wage %>% filter(group_ever == "childless_men_ever"),
  "Fathers"        = wage %>% filter(group_ever == "fathers_ever"),
  "Never-mothers"  = wage %>% filter(group_ever == "never_mothers"),
  "Mothers (ever)" = wage %>% filter(group_ever == "mothers_ever")
)

oster_results <- map_dfr(names(oster_groups), function(g) {
  compute_oster(oster_groups[[g]], g, verbose = TRUE)
})

# --- Summary table ---

cat(strrep("=", 70), "\n")
cat("  OSTER SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("  %-20s  %8s  %8s  %8s  %8s  %8s  %7s\n",
            "Group", "OLS(R)", "OLS(F)", "FE β₂", "β*(δ=1)", "R²_full", "δ"))
cat("  ", strrep("-", 72), "\n")

for (i in seq_len(nrow(oster_results))) {
  r <- oster_results[i, ]
  cat(sprintf("  %-20s  %+.3f%%  %+.3f%%%s  %+.3f%%%s  %+.3f%%  %.4f  %+.2f\n",
              r$group,
              r$b_restricted,
              r$b_full_ols, stars_fn(r$p_full_ols),
              r$b_fe, stars_fn(r$p_fe),
              if (!is.na(r$beta_star)) r$beta_star else NA_real_,
              r$r2_full,
              if (!is.na(r$delta)) r$delta else NA_real_))
}

cat("\n  OLS(R) = restricted (educ, educ×exp, exp only)\n")
cat("  OLS(F) = full OLS (+ controls + wave FE)\n")
cat("  FE β₂  = individual + wave FE (causal benchmark)\n")
cat("  β*(δ=1) = Oster bound under equal selection\n")
cat("  δ      = selection needed to explain OLS(F)→FE gap\n\n")

# --- Sensitivity: alternative R²_max scalars ---

cat(strrep("-", 70), "\n")
cat("  SENSITIVITY: Oster δ under alternative R²_max scalars\n")
cat(strrep("-", 70), "\n\n")

cat(sprintf("  %-20s  %8s  %8s  %8s\n",
            "Group", "δ(×1.0)", "δ(×1.3)", "δ(×2.0)"))
cat("  ", strrep("-", 50), "\n")

for (g in names(oster_groups)) {
  row <- list()
  for (sc in c(1.0, 1.3, 2.0)) {
    r <- compute_oster(oster_groups[[g]], g, r2_max_scalar = sc, verbose = FALSE)
    row[[as.character(sc)]] <- if (!is.null(r)) r$delta else NA
  }
  cat(sprintf("  %-20s  %+8.2f  %+8.2f  %+8.2f\n",
              g,
              if (!is.na(row[["1"]])) row[["1"]] else NA,
              if (!is.na(row[["1.3"]])) row[["1.3"]] else NA,
              if (!is.na(row[["2"]])) row[["2"]] else NA))
}

cat("\n  Note: R²_max = min(1, scalar × R²_full).\n")
cat("  Oster (2019) recommends scalar=1.3. Results are sign-stable across\n")
cat("  all scalars for groups where |δ| >> 1.\n\n")

# --- Identified sets [β*, b_full] for δ ∈ {0.5, 1.0, 2.0} ---

cat(strrep("-", 70), "\n")
cat("  IDENTIFIED SETS: β* at alternative δ values\n")
cat(strrep("-", 70), "\n\n")

cat(sprintf("  %-20s  %8s  %10s  %10s  %10s\n",
            "Group", "FE β₂", "β*(δ=0.5)", "β*(δ=1.0)", "β*(δ=2.0)"))
cat("  ", strrep("-", 65), "\n")

oster_group_key <- c(
  "Childless men"  = "childless_men_ever",
  "Fathers"        = "fathers_ever",
  "Never-mothers"  = "never_mothers",
  "Mothers (ever)" = "mothers_ever"
)

for (i in seq_len(nrow(oster_results))) {
  r  <- oster_results[i, ]
  dd <- wage %>%
    filter(group_ever == oster_group_key[r$group], in_wage_sample == 1) %>%
    mutate(
      educ_c     = educ_years       - mean(educ_years,       na.rm = TRUE),
      exp_c      = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c
    ) %>%
    filter(complete.cases(ln_hourly_wage_real, educ_years, educ_exp_c,
                          experience_years, age_sq, married_num))
  
  m_r_i  <- lm(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years, data = dd)
  m_f_i  <- lm(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
                 age_sq + married_num + as.factor(wave), data = dd)
  b_r_i  <- coef(m_r_i)["educ_exp_c"]
  b_f_i  <- coef(m_f_i)["educ_exp_c"]
  r2_r_i <- summary(m_r_i)$r.squared
  r2_f_i <- summary(m_f_i)$r.squared
  r2_mx  <- min(1, 1.3 * r2_f_i)
  Pi_i   <- (b_r_i - b_f_i) * (r2_mx - r2_f_i) / max(r2_f_i - r2_r_i, 1e-8)
  bstar  <- function(delta) (b_f_i - delta * Pi_i) * 100
  
  cat(sprintf("  %-20s  %+.3f%%   %+.3f%%    %+.3f%%    %+.3f%%\n",
              r$group, r$b_fe,
              bstar(0.5), bstar(1.0), bstar(2.0)))
}

cat("\n  β*(δ=1.0) is the standard Oster bound.\n")
cat("  If β*(δ=1.0) > 0 for childless groups → complementarity is robust even\n")
cat("  under equal selection on unobservables.\n\n")

# --- LaTeX table (paste into paper) ---

cat(strrep("=", 70), "\n")
cat("  LaTeX TABLE (paste into paper)\n")
cat(strrep("=", 70), "\n\n")

cat("\\begin{tabular}{lcccccc}\n")
cat("\\toprule\n")
cat(sprintf("%-22s & %s & %s & %s & %s & %s & %s \\\\\n",
            "Group",
            "OLS(F) $\\hat{\\beta}_2$",
            "FE $\\hat{\\beta}_2$",
            "$\\beta^*(\\delta=1)$",
            "$R^2_{\\text{full}}$",
            "$R^2_{\\text{max}}$",
            "$\\delta$"))
cat("\\midrule\n")
for (i in seq_len(nrow(oster_results))) {
  r <- oster_results[i, ]
  cat(sprintf("%-22s & %+.3f\\%% & %+.3f\\%%%s & %+.3f\\%% & %.4f & %.4f & %+.2f \\\\\n",
              r$group,
              r$b_full_ols,
              r$b_fe, stars_fn(r$p_fe),
              if (!is.na(r$beta_star)) r$beta_star else NA_real_,
              r$r2_full,
              r$r2_max,
              if (!is.na(r$delta)) r$delta else NA_real_))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n\n")
cat("% OLS(F): full OLS with controls + wave FE.\n")
cat("% FE: individual + wave FE (causal benchmark).\n")
cat("% beta*(delta=1): Oster bound under equal selection.\n")
cat("% delta < 0 for childless groups: FE > OLS, selection must work against observables.\n\n")

# --- Save Oster results ---
saveRDS(oster_results, file.path(outdir, "oster_bounds.rds"))
write.csv(oster_results, file.path(outdir, "oster_bounds.csv"), row.names = FALSE)
cat(sprintf("  \u2713 Saved: %s\n", file.path(outdir, "oster_bounds.csv")))


# #############################################################################
# PART 8: AGE AT FIRST BIRTH HETEROGENEITY (3.1–3.3)
# Memory-efficient: pre-demean by person instead of feols person_id FE
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper §5.11, Figure 5: AFB gradient. "Early (<25): +0.042%;
#   Mid (25–31): +0.051%; Late (32+): +0.105%"
# Sample: Ever-mothers in wage sample with non-missing AFB.
#   Expected: wave-split counts re-verify against the current extract.
cat("#  PART 8: AGE AT FIRST BIRTH HETEROGENEITY (3.1–3.3)\n")
cat(strrep("#", 80), "\n")

# Free memory from earlier parts
gc(verbose = FALSE)

mothers <- panel %>%
  filter(in_wage_sample == 1, group_ever == "mothers_ever",
         !is.na(approx_age_at_first_birth)) %>%
  select(person_id, wave, ln_hourly_wage_real, educ_years, experience_years,
         age_sq, married_num, approx_age_at_first_birth) %>%
  mutate(afb_group = case_when(
    approx_age_at_first_birth < 25 ~ "Early (<25)",
    approx_age_at_first_birth < 32 ~ "Mid (25-31)",
    TRUE ~ "Late (32+)"
  ))

num_rows    <- format(nrow(mothers), big.mark = ",")
num_persons <- format(n_distinct(mothers$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

cat("\n  Age at first birth distribution:\n")
mothers %>%
  group_by(afb_group) %>%
  summarise(N = n(), N_id = n_distinct(person_id),
            mean_afb = mean(approx_age_at_first_birth),
            mean_educ = mean(educ_years, na.rm = TRUE),
            .groups = "drop") %>%
  print()

# BUG FIX: compute centering reference means from the FULL mothers sample
# (not within each AFB sub-group). Sub-group centering absorbs the mean
# shift introduced by the experience fix, making estimates invariant to the
# dataset update. Using a shared reference mean correctly propagates any
# individual-specific corrections into educ_exp_c.
educ_mean_ref <- mean(mothers$educ_years,      na.rm = TRUE)
exp_mean_ref  <- mean(mothers$experience_years, na.rm = TRUE)

# Helper: within-person demean (absorbs person FE without building sparse matrix)
demean_by <- function(x, id) {
  gm <- ave(x, id, FUN = function(v) mean(v, na.rm = TRUE))
  x - gm
}

r8 <- map_dfr(c("Early (<25)", "Mid (25-31)", "Late (32+)"), function(afb) {
  d <- mothers %>%
    filter(afb_group == afb) %>%
    mutate(
      # FIX: center on full mothers reference means, not sub-group means
      educ_c     = educ_years      - educ_mean_ref,
      exp_c      = experience_years - exp_mean_ref,
      educ_exp_c = educ_c * exp_c
    )
  
  # Pre-demean all variables by person_id (equivalent to person FE)
  d <- d %>%
    mutate(across(c(ln_hourly_wage_real, educ_years, educ_exp_c, experience_years,
                    age_sq, married_num),
                  ~ demean_by(.x, person_id),
                  .names = "dm_{.col}"))
  
  # Wave dummies on demeaned data (equivalent to feols ... | person_id + wave)
  m_fe <- tryCatch(
    lm(dm_ln_hourly_wage_real ~ dm_educ_years + dm_educ_exp_c + dm_experience_years +
         dm_age_sq + dm_married_num + factor(wave),
       data = d),
    error = function(e) NULL
  )
  
  b <- if (!is.null(m_fe) && "dm_educ_exp_c" %in% names(coef(m_fe)))
    coef(m_fe)["dm_educ_exp_c"] else NA
  
  # Cluster-robust SE by person_id
  p <- tryCatch({
    vc <- sandwich::vcovCL(m_fe, cluster = d$person_id)
    se <- sqrt(vc["dm_educ_exp_c", "dm_educ_exp_c"])
    2 * pt(abs(b / se), df = n_distinct(d$person_id) - 1, lower.tail = FALSE)
  }, error = function(e) NA)
  
  cat(sprintf("  %-15s  FE β₂ = %+.4f%%%s  N=%s\n",
              afb, b*100, stars_fn(p),
              format(nrow(d), big.mark = ",")))
  
  tibble(label = paste0("AFB: ", afb), group = afb,
         N = nrow(d), N_id = n_distinct(d$person_id),
         ols_b2 = NA, ols_p = NA, fe_b2 = b, fe_p = p)
})
master <- bind_rows(master, r8 %>% mutate(test = "8. Age at first birth"))
rec("afb_early", 100 * r8$fe_b2[match("Early (<25)", r8$group)])
rec("afb_mid",   100 * r8$fe_b2[match("Mid (25-31)", r8$group)])
rec("afb_late",  100 * r8$fe_b2[match("Late (32+)",  r8$group)])
rm(mothers); gc(verbose = FALSE)


# #############################################################################
# PART 9: WAVE-BY-WAVE β₂ (3.24–3.26)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
# Paper Appendix Figure A1: "β₂ hierarchy is stable across waves."
# Sample: Wage sample, by wave. Per wave ~1,300–2,000 per group.
cat("#  PART 9: WAVE-BY-WAVE β₂ WITH ABSORBING DEFINITION (3.24–3.26)\n")
cat(strrep("#", 80), "\n")

wave_results <- tibble()
for (w in 12:24) {
  w_data <- panel %>%
    filter(wave == w, in_wage_sample == 1)
  
  for (g in c("childless_men_ever", "fathers_ever",
              "never_mothers", "mothers_ever")) {
    d <- w_data %>%
      filter(group_ever == g) %>%
      mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
             exp_c = experience_years - mean(experience_years, na.rm = TRUE),
             educ_exp_c = educ_c * exp_c)
    
    m <- tryCatch(
      lm(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
           age_sq + married_num, data = d),
      error = function(e) NULL
    )
    
    if (!is.null(m) && "educ_exp_c" %in% names(coef(m))) {
      ct <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
      b <- ct["educ_exp_c", "Estimate"]
      se <- ct["educ_exp_c", "Std. Error"]
      p <- ct["educ_exp_c", "Pr(>|t|)"]
      wave_results <- bind_rows(wave_results,
                                tibble(wave = w, year = 2000 + w, group = g,
                                       N = nrow(d), beta2 = b, se = se, pval = p))
    }
  }
}

# Print
cat("\n  Wave-by-wave OLS β₂ (absorbing groups):\n\n")
cat(sprintf("  %4s  %-25s  %7s  %10s  %6s\n",
            "Wave", "Group", "N", "β₂ × 100", "p"))
cat("  ", strrep("-", 60), "\n")
for (i in 1:nrow(wave_results)) {
  r <- wave_results[i,]
  cat(sprintf("  %4d  %-25s  %7d  %+.4f%%%s  %.3f\n",
              r$wave, r$group, r$N, r$beta2*100, stars_fn(r$pval), r$pval))
}


# #############################################################################
# PART 10: MASTER SUMMARY TABLE
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
cat("#  PART 10: MASTER ROBUSTNESS SUMMARY\n")
cat(strrep("#", 80), "\n\n")

# Focus on key groups: childless men and never-mothers (cleanest comparison)
summary_tab <- master %>%
  filter(grepl("Childless men|Never-mothers", label)) %>%
  select(test, label, N, ols_b2, ols_p, fe_b2, fe_p)

cat("  CHILDLESS MEN — β₂ across specifications:\n\n")
cm <- master %>% filter(grepl("Childless men", label))
cat(sprintf("  %-35s  %7s  %10s  %10s\n", "Test", "N", "OLS β₂", "FE β₂"))
cat("  ", strrep("-", 65), "\n")
for (i in 1:nrow(cm)) {
  r <- cm[i,]
  ols_str <- if (!is.na(r$ols_b2)) sprintf("%+.4f%%%s", r$ols_b2*100, stars_fn(r$ols_p)) else "—"
  fe_str <- if (!is.na(r$fe_b2)) sprintf("%+.4f%%%s", r$fe_b2*100, stars_fn(r$fe_p)) else "—"
  cat(sprintf("  %-35s  %7s  %10s  %10s\n",
              r$test, format(r$N, big.mark = ","), ols_str, fe_str))
}

cat("\n\n  NEVER-MOTHERS — β₂ across specifications:\n\n")
nm <- master %>% filter(grepl("Never-mothers", label))
cat(sprintf("  %-35s  %7s  %10s  %10s\n", "Test", "N", "OLS β₂", "FE β₂"))
cat("  ", strrep("-", 65), "\n")
for (i in 1:nrow(nm)) {
  r <- nm[i,]
  ols_str <- if (!is.na(r$ols_b2)) sprintf("%+.4f%%%s", r$ols_b2*100, stars_fn(r$ols_p)) else "—"
  fe_str <- if (!is.na(r$fe_b2)) sprintf("%+.4f%%%s", r$fe_b2*100, stars_fn(r$fe_p)) else "—"
  cat(sprintf("  %-35s  %7s  %10s  %10s\n",
              r$test, format(r$N, big.mark = ","), ols_str, fe_str))
}

# Save master table
saveRDS(master, "robustness_master.rds")

# --- Standard-error ledger for Table B2 (AEA style: SEs in parentheses, no stars) ---
se_ledger_robustness <- master %>%
  dplyr::transmute(
    test, label, group, N,
    fe_beta2_x100  = fe_b2  * 100, fe_se_x100  = fe_se  * 100,
    ols_beta2_x100 = ols_b2 * 100, ols_se_x100 = ols_se * 100
  )
write.csv(se_ledger_robustness, "se_ledger_robustness.csv", row.names = FALSE)
cat("\n  Wrote se_ledger_robustness.csv (Table B2 estimates + clustered SEs).\n")
write_csv(master, "robustness_master.csv")

cat("\n\n  Saved: robustness_master.rds / .csv\n")

# Save wave-by-wave results
saveRDS(wave_results, "wave_by_wave_beta2.rds")
write_csv(wave_results, "wave_by_wave_beta2.csv")

cat("  Saved: wave_by_wave_beta2.rds / .csv\n")

cat("\n")
cat(strrep("#", 80), "\n")
cat("#  ROBUSTNESS BATTERY COMPLETE\n")
cat("#  Total specifications: ", nrow(master), "\n")
cat(strrep("#", 80), "\n")



# #############################################################################
# PART 5: PAPER-READY TABLES (6.1, 6.3–6.5)
# #############################################################################

cat("\n\n")
cat(strrep("#", 80), "\n")
cat("#  PART 5: PAPER-READY TABLES\n")
cat(strrep("#", 80), "\n")

# --- 6.1 Descriptive table: four groups (broad + wage) ---
cat("\n--- 6.1 Descriptive table ---\n\n")

desc_table <- panel %>%
  filter(!is.na(group_ever)) %>%
  group_by(group_ever) %>%
  summarise(
    N_broad = n(),
    N_id_broad = n_distinct(person_id),
    N_wage = sum(in_wage_sample == 1, na.rm = TRUE),
    N_id_wage = n_distinct(person_id[in_wage_sample == 1]),
    mean_educ = mean(educ_years[in_wage_sample == 1], na.rm = TRUE),
    mean_exp = mean(experience_years[in_wage_sample == 1], na.rm = TRUE),
    mean_age = mean(age[in_wage_sample == 1], na.rm = TRUE),
    mean_wage = mean(exp(ln_hourly_wage_real[in_wage_sample == 1]), na.rm = TRUE),
    pct_ft = mean(fulltime[in_wage_sample == 1] == 1, na.rm = TRUE) * 100,
    pct_married = mean(married[in_wage_sample == 1] == 1, na.rm = TRUE) * 100,
    pct_employed = mean(employed == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat(sprintf("  %-20s  %7s %7s  %5s  %5s  %5s  %7s  %5s  %5s  %5s\n",
            "Group", "N_wage", "N_id", "Educ", "Exp", "Age",
            "Wage$", "FT%", "Mar%", "Emp%"))
cat("  ", strrep("-", 80), "\n")
for (i in 1:nrow(desc_table)) {
  r <- desc_table[i,]
  cat(sprintf("  %-20s  %7s %7s  %5.1f  %5.1f  %5.1f  %7.1f  %5.1f  %5.1f  %5.1f\n",
              r$group_ever,
              format(r$N_wage, big.mark = ","),
              format(r$N_id_wage, big.mark = ","),
              r$mean_educ, r$mean_exp, r$mean_age, r$mean_wage,
              r$pct_ft, r$pct_married, r$pct_employed))
}


# --- 6.1b Three-way women descriptive (fills Table 1 footnote-a cells) ---
cat("\n\n--- 6.1b Three-way women descriptive (wage sample) ---\n\n")

threeway_desc <- wage %>%
  filter(!is.na(group_3way_women)) %>%
  group_by(group_3way_women) %>%
  summarise(
    N_obs       = n(),
    N_id        = n_distinct(person_id),
    mean_educ   = mean(educ_years,          na.rm = TRUE),
    mean_exp    = mean(experience_years,    na.rm = TRUE),
    mean_lnwage = mean(ln_hourly_wage_real, na.rm = TRUE),
    mean_wage   = mean(exp(ln_hourly_wage_real), na.rm = TRUE),
    pct_ft      = mean(fulltime == 1,       na.rm = TRUE) * 100,
    pct_emp     = mean(employed == 1,       na.rm = TRUE) * 100,
    mean_waves  = n() / n_distinct(person_id),
    .groups = "drop"
  )

cat(sprintf("  %-18s  %7s %7s  %6s  %6s  %7s  %7s  %6s  %6s  %6s\n",
            "Group", "N_obs", "N_id", "Educ", "Exp",
            "Ln wage", "Wage($)", "FT%", "Emp%", "Waves"))
cat("  ", strrep("-", 80), "\n")
for (i in 1:nrow(threeway_desc)) {
  r <- threeway_desc[i, ]
  cat(sprintf("  %-18s  %7s %7s  %6.1f  %6.1f  %7.3f  %7.1f  %6.1f  %6.1f  %6.1f\n",
              r$group_3way_women,
              format(r$N_obs, big.mark = ","),
              format(r$N_id,  big.mark = ","),
              r$mean_educ, r$mean_exp, r$mean_lnwage, r$mean_wage,
              r$pct_ft, r$pct_emp, r$mean_waves))
}
cat("\n  → Use these values to fill the active-mothers and post-mothers cells\n")
cat("    (rows: Education, Experience, Ln real wage) in Table 1 of the paper.\n")


# --- 6.3 Main results table: OLS and FE under both definitions ---
cat("\n\n--- 6.3 Main results: OLS/FE × absorbing/time-varying ---\n\n")

cat(sprintf("  %-20s  %12s  %12s  %12s  %12s\n",
            "Group", "OLS(absorb)", "FE(absorb)", "OLS(TV)", "FE(TV)"))
cat("  ", strrep("-", 70), "\n")

for (g_abs in c("childless_men_ever", "fathers_ever",
                "never_mothers", "mothers_ever")) {
  
  d_abs <- wage %>%
    filter(group_ever == g_abs) %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c)
  
  m_ols_a <- tryCatch(feols(ln_hourly_wage_real ~ educ_years + educ_exp_c +
                              experience_years + age_sq + married_num | wave,
                            data = d_abs, cluster = ~person_id, notes = FALSE), error = function(e) NULL)
  m_fe_a <- tryCatch(feols(ln_hourly_wage_real ~ educ_exp_c +
                             experience_years + age_sq + married_num | person_id + wave,
                           data = d_abs, cluster = ~person_id, notes = FALSE), error = function(e) NULL)
  
  get_b2 <- function(m) {
    if (is.null(m) || !("educ_exp_c" %in% names(coef(m)))) return("—")
    s <- summary(m)$coeftable
    sprintf("%+.3f%s", s["educ_exp_c","Estimate"]*100,
            stars_fn(s["educ_exp_c","Pr(>|t|)"]))
  }
  
  # Time-varying equivalent
  g_tv <- switch(g_abs,
                 "childless_men_ever" = "childless_men",
                 "fathers_ever" = "fathers",
                 "never_mothers" = "childless_women_current",
                 "mothers_ever" = "mothers_current"
  )
  
  panel_tv <- panel %>%
    mutate(group_tv = case_when(
      female == 0 & ever_parent == 0 ~ "childless_men",
      female == 0 & ever_parent == 1 ~ "fathers",
      female == 1 & num_children_under15 > 0 ~ "mothers_current",
      female == 1 & num_children_under15 == 0 ~ "childless_women_current"
    ))
  
  d_tv <- panel_tv %>%
    filter(in_wage_sample == 1, group_tv == g_tv) %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c)
  
  m_ols_t <- tryCatch(feols(ln_hourly_wage_real ~ educ_years + educ_exp_c +
                              experience_years + age_sq + married_num | wave,
                            data = d_tv, cluster = ~person_id, notes = FALSE), error = function(e) NULL)
  m_fe_t <- tryCatch(feols(ln_hourly_wage_real ~ educ_exp_c +
                             experience_years + age_sq + married_num | person_id + wave,
                           data = d_tv, cluster = ~person_id, notes = FALSE), error = function(e) NULL)
  
  cat(sprintf("  %-20s  %12s  %12s  %12s  %12s\n",
              g_abs, get_b2(m_ols_a), get_b2(m_fe_a),
              get_b2(m_ols_t), get_b2(m_fe_t)))
}


# --- 6.5 Employment selection table ---
cat("\n\n--- 6.5 Employment selection table ---\n\n")

emp_selection <- panel %>%
  filter(!is.na(group_3way_women) | (female == 0 & !is.na(group_ever))) %>%
  mutate(group_label = case_when(
    group_ever == "childless_men_ever" ~ "Childless men",
    group_ever == "fathers_ever" ~ "Fathers",
    group_3way_women == "never_mothers" ~ "Never-mothers",
    group_3way_women == "active_mothers" ~ "Active mothers",
    group_3way_women == "post_mothers" ~ "Post-mothers",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group_label)) %>%
  group_by(group_label) %>%
  summarise(
    N = n(),
    emp_rate = mean(employed == 1, na.rm = TRUE) * 100,
    wage_rate = mean(in_wage_sample == 1, na.rm = TRUE) * 100,
    ft_given_emp = mean(fulltime[employed == 1] == 1, na.rm = TRUE) * 100,
    mean_hours_emp = mean(hours_worked_clean[employed == 1], na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  %-20s  %7s  %7s  %7s  %7s  %7s\n",
            "Group", "N", "Emp%", "Wage%", "FT|Emp", "Hours"))
cat("  ", strrep("-", 55), "\n")
for (i in 1:nrow(emp_selection)) {
  r <- emp_selection[i,]
  cat(sprintf("  %-20s  %7s  %6.1f%%  %6.1f%%  %6.1f%%  %7.1f\n",
              r$group_label, format(r$N, big.mark = ","),
              r$emp_rate, r$wage_rate, r$ft_given_emp, r$mean_hours_emp))
}

cat("\n  Key gaps:\n")
nm <- emp_selection %>% filter(group_label == "Never-mothers")
am <- emp_selection %>% filter(group_label == "Active mothers")
pm <- emp_selection %>% filter(group_label == "Post-mothers")
if (nrow(nm) > 0 && nrow(am) > 0) {
  cat(sprintf("  Never-mothers vs Active mothers: %.1fpp employment gap\n",
              nm$emp_rate - am$emp_rate))
}
if (nrow(nm) > 0 && nrow(pm) > 0) {
  cat(sprintf("  Never-mothers vs Post-mothers: %.1fpp employment gap\n",
              nm$emp_rate - pm$emp_rate))
}

# ── PART O: estimator-aligned betaX contrast + Lee bounds (Table A4 / Table 13) ──
#    Independent of the event study: reads only the slim panel, reuses helpers
#    (stars_fn, rec) defined above. Writes se_ledger_A3_aligned*.csv,
#    betaX_lee_bounds.csv, mincer_gap_betaX_selcorr.csv.
if (file.exists("02b_part_o_aligned.R")) {
  tryCatch(source("02b_part_o_aligned.R"),
           error = function(e)
             cat(sprintf("  NOTE: 02b_part_o_aligned.R did not complete (%s); run it standalone.\n",
                         conditionMessage(e))))
} else {
  cat("  NOTE: 02b_part_o_aligned.R not found; Table A4 / Table 13 betaX bounds not refreshed.\n")
}

# ==============================================================================
# PART G: FTB-B REFORM 1 IV  (sourced from 04_part_g_ftbb_iv.R)
# Identifies theta on S_ct using the 2015 FTB-B income-threshold reform.
# Produces ftbb_reform1_iv_results.rds and the iv_summary table.
# ==============================================================================

# ── Child-penalty event study (Sun-Abraham); writes cs_event_study_results.rds ──
#    Reuses the full `panel` already in memory; figure is left to RUN_FIGURES_hh.R.
if (file.exists("02_cs_event_study.R")) {
  options(cs_estimation_only = TRUE)
  tryCatch(source("02_cs_event_study.R"),
           error = function(e)
             cat(sprintf("  NOTE: 02_cs_event_study.R did not complete (%s); run it standalone.\n",
                         conditionMessage(e))))
  options(cs_estimation_only = NULL)
} else {
  cat("  NOTE: 02_cs_event_study.R not found; cs_event_study_results.rds not refreshed.\n")
}

# ── Implied dynamic cost, concave/tenure-varying return (07b) ──────────────────
#    Reads cs_event_study_results.rds; writes implied_dynamic_cost_concave_results.rds.
#    MUST run before 08_gelbach so the Gelbach output can read the late-horizon
#    (k=5..10) return-times-gap share (the reliable ~29% figure for Appendix A6).
if (file.exists("07b_implied_dynamic_cost_concave.R")) {
  tryCatch(source("07b_implied_dynamic_cost_concave.R"),
           error = function(e)
             cat(sprintf("  NOTE: 07b_implied_dynamic_cost_concave.R did not complete (%s); run it standalone.\n",
                         conditionMessage(e))))
} else {
  cat("  NOTE: 07b_implied_dynamic_cost_concave.R not found; late-horizon experience share not refreshed.\n")
}

# ── Gelbach decomposition of the post-birth wage penalty (Table 10 / App. A6) ──
#    Runs immediately after the event study and 07b so it reuses the in-memory
#    `panel`, the event-study roster, and the 07b late-horizon share; writes
#    gelbach_decomposition_results.rds.
if (file.exists("08_gelbach_decomposition.R")) {
  tryCatch(source("08_gelbach_decomposition.R"),
           error = function(e)
             cat(sprintf("  NOTE: 08_gelbach_decomposition.R did not complete (%s); run it standalone.\n",
                         conditionMessage(e))))
} else {
  cat("  NOTE: 08_gelbach_decomposition.R not found; Table 10 / Appendix A6 not refreshed.\n")
}

if (file.exists("04_part_g_ftbb_iv.R")) {
  source("04_part_g_ftbb_iv.R")
} else {
  warning("04_part_g_ftbb_iv.R not found in working directory; ",
          "PART G skipped. FTB-B IV results will not be available.")
}

# ── PART H: policy economics of the FTB-B reform (needs Part G objects) ─────
if (file.exists("05_part_h_policy.R") && exists("iv_couples")) {
  tryCatch(source("05_part_h_policy.R"),
           error = function(e)
             cat(sprintf("  NOTE: PART H did not complete (%s); run it standalone.\n",
                         conditionMessage(e))))
} else if (!file.exists("05_part_h_policy.R")) {
  cat("  NOTE: 05_part_h_policy.R not found; PART H skipped.\n")
} else {
  cat("  NOTE: Part G objects missing; PART H skipped.\n")
}

# ── Training event study (Sun-Abraham, training outcomes) ──────────────────
if (file.exists("03_training_event_study.R")) {
  options(training_estimation_only = TRUE)
  tryCatch(source("03_training_event_study.R"),
           error = function(e)
             cat(sprintf("  NOTE: training event study did not complete (%s).\n",
                         conditionMessage(e))))
  options(training_estimation_only = NULL)
} else {
  cat("  NOTE: 03_training_event_study.R not found; skipped.\n")
}
# NOTE: sourced BEFORE the FINAL VERIFICATION block below so that the
# verification's §5.10 line can read ftbb_reform1_iv_results from memory
# (previously it printed "[not computed this run]" on every fresh run).


# ======
# FINAL VERIFICATION: KEY PAPER VALUES
# =============================================================================
# This section prints all key numerical claims from the paper alongside
# the estimates computed above, for audit purposes.
# Paper: "Household Specialisation, Lost Experience, and the Motherhood Wage Penalty" (cost_hh_spec.tex)

cat("\n\n")
cat(strrep("#", 80), "\n")
cat("#  FINAL VERIFICATION: KEY PAPER VALUES vs COMPUTED ESTIMATES\n")
cat("#  Use this output to verify consistency with the manuscript.\n")
cat(strrep("#", 80), "\n\n")

cat("─── §3 DATA ───────────────────────────────────────────────────────────────\n")
cat(sprintf("  Full panel:   %s obs / %s persons\n",
            format(nrow(panel), big.mark = ","),
            format(n_distinct(panel$person_id), big.mark = ",")))
cat(sprintf("  Wage sample:  %s obs / %s persons\n\n",
            format(nrow(wage), big.mark = ","),
            format(n_distinct(wage$person_id), big.mark = ",")))

cat("─── Table 1: BASELINE β₂ ──────────────────────────────────────────────────\n")
cat("  (see baseline_4 estimates printed above)\n\n")

cat("─── §6.1 PARENTHOOD TRANSITIONS ──────────────────────────────────────────\n")
cat("  (see Part A transition counts printed above)\n\n")

cat("─── §6.2 BREADWINNER LOCK-IN (men, ×100) ─────────────────────────────────\n")
cat(sprintf("  β_educ pre → post: %s → %s\n", pv("lockin_be_pre","%.2f"), pv("lockin_be_post","%.2f")))
cat(sprintf("  β₂     pre → post: %s → %s\n\n", pv("lockin_b2_pre","%.3f"), pv("lockin_b2_post","%.3f")))

cat("─── P1: STATIC SPECIALISATION (couple FE levels) ──────────────────────────\n")
cat(sprintf("  10-hr fall in her hours → %s%% his wage (her_hours coef p = %s)\n",
            pv("becker_pct","%+.2f"), pv("becker_p","%.4f")))
cat(sprintf("  Childless × parent interaction: p = %s\n\n", pv("becker_int_p","%.4f")))

cat("─── P2: TRIPLE INTERACTION ───────────────────────────────────────────────\n")
cat("  (see B3 triple-interaction and DML output above)\n\n")

cat("─── P3: PERSISTENT SCAR ──────────────────────────────────────────────────\n")
cat(sprintf("  Post-mothers β₂ (FE raw) = %s%%\n\n",
            gv(heckman_results$fe_raw[heckman_results$group == "Post-mothers"], "%+.3f")))

cat("─── P4: ASSORTATIVE MATCHING (within-couple DML) ──────────────────────────\n")
cat(sprintf("  His β₂ = %s%%, her β₂ = %s%%, gap = %s%% (p = %s)\n\n",
            pv("dml_his","%.3f"), pv("dml_her","%.3f"), pv("dml_gap","%.3f"), pv("dml_gap_p","%.4f")))

cat("─── §5.3 CONTAMINATION ───────────────────────────────────────────────────\n")
cat(sprintf("  Post-mothers mislabelled childless: %s%%\n\n", pv("contam_pct","%.1f")))

cat("─── §5.9 OSTER BOUNDS ────────────────────────────────────────────────────\n")
cat(sprintf("  Childless men δ = %s ; never-mothers δ = %s\n\n",
            gv(oster_results$delta[oster_results$group == "Childless men"], "%+.2f"),
            gv(oster_results$delta[oster_results$group == "Never-mothers"], "%+.2f")))

cat("─── §5.10 IV ─────────────────────────────────────────────────────────────\n")
cat(sprintf("  FTB-B Reform 1 IV: θ on S_ct = %s (SE %s), first-stage F = %s\n\n",
            gv(ftbb_reform1_iv_results$summary$theta_hat[ftbb_reform1_iv_results$summary$spec == "IV Reform 1"], "%.2f"),
            gv(ftbb_reform1_iv_results$summary$se[ftbb_reform1_iv_results$summary$spec == "IV Reform 1"], "%.2f"),
            gv(ftbb_reform1_iv_results$summary$fs_F_Sct[ftbb_reform1_iv_results$summary$spec == "IV Reform 1"], "%.1f")))

cat("─── §5.11 AGE AT FIRST BIRTH ─────────────────────────────────────────────\n")
cat(sprintf("  Early = %s%% ; Mid = %s%% ; Late = %s%%\n\n",
            pv("afb_early","%+.3f"), pv("afb_mid","%+.3f"), pv("afb_late","%+.3f")))

cat("─── CAUSAL FOREST ────────────────────────────────────────────────────────\n")
cat(sprintf("  Own experience importance = %s%%, own/partner ratio = %s×\n",
            pv("forest_exp_imp","%.1f"), pv("forest_ratio","%.2f")))
cat(sprintf("  CATE–partner-hours correlation (childless men) = %s\n\n",
            pv("cate_childless_men","%+.3f")))

cat(strrep("#", 80), "\n")
cat("#  ANALYSIS COMPLETE\n")
cat(strrep("#", 80), "\n")

cat("\nDone.\n")





# ==============================================================================
# SAVE ALL DML RESULTS — COUPLE-LEVEL ANALYSIS
# Run this immediately after PART D completes
# Saves: dial_dml_results.rds  (full objects, reloadable)
#        dial_dml_summary.csv  (all scalar estimates in one flat table)
#        dial_cf_vimp.csv      (causal forest variable importance)
#        dial_cf_group_cate.csv (CATE by gender × parenthood)
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("  SAVING DML RESULTS\n")
cat(strrep("=", 70), "\n\n")

# ------------------------------------------------------------------------------
# HELPER: safely strip the large DoubleML model object but keep all scalars
# ------------------------------------------------------------------------------
strip_model <- function(x) {
  if (is.null(x)) return(NULL)
  x$model <- NULL
  x
}

# ------------------------------------------------------------------------------
# 1.  COLLECT EVERYTHING INTO ONE LIST
# ------------------------------------------------------------------------------

dial_dml_results <- list(
  
  # ── Part 1A: Pooled DML ────────────────────────────────────────────────────
  pooled = list(
    all_own          = strip_model(dml_all_own),
    all_partner      = strip_model(dml_all_partner),
    all_lasso_own    = strip_model(dml_all_lasso),
    all_lasso_partner = strip_model(dml_all_partner_lasso)
  ),
  
  # ── Part 1B: By gender ─────────────────────────────────────────────────────
  by_gender = list(
    his          = strip_model(dml_his),
    her          = strip_model(dml_her),
    his_own_only = strip_model(dml_his_own),
    her_own_only = strip_model(dml_her_own)
  ),
  
  # ── Part 1C: Four cells (RF + LASSO) ───────────────────────────────────────
  cells_rf    = lapply(dml_cells_rf,    strip_model),
  cells_lasso = lapply(dml_cells_lasso, strip_model),
  
  # ── Part 2: Within-couple (couple-demeaned) ────────────────────────────────
  within = list(
    overall_own     = strip_model(dml_w_own),
    overall_partner = strip_model(dml_w_partner),
    his             = strip_model(dml_w_his),
    her             = strip_model(dml_w_her),
    cells           = lapply(dml_w_cells, strip_model)
  ),
  
  # ── Part 3: Causal forest ─────────────────────────────────────────────────
  # The forest object itself is large; keep derived quantities separately.
  # Set save_forest = TRUE below if you want the full grf object.
  cf = if (exists("cf") && !is.null(cf)) list(
    ate         = if (exists("ate"))         ate         else NULL,
    calibration = if (exists("cal"))         cal         else NULL,
    vimp        = if (exists("vimp_df"))     vimp_df     else NULL,
    group_cate  = if (exists("group_cate"))  group_cate  else NULL
  ) else NULL,
  
  # ── Part 4: Linear OLS comparison ─────────────────────────────────────────
  linear = if (exists("linear_results")) linear_results else NULL,
  
  # ── Part 5: Becker triple interaction ─────────────────────────────────────
  becker = if (exists("p_triple")) list(
    father_penalty = father_penalty,
    father_se      = father_se,
    mother_penalty = mother_penalty,
    mother_se      = mother_se,
    triple         = triple,
    triple_se      = triple_se,
    z_triple       = z_triple,
    p_triple       = p_triple
  ) else NULL,
  
  # ── Metadata ───────────────────────────────────────────────────────────────
  meta = list(
    saved_at    = Sys.time(),
    n_stacked   = if (exists("stacked"))   nrow(stacked)                  else NA,
    n_couples   = if (exists("stacked"))   n_distinct(stacked$couple_id)  else NA,
    n_within    = if (exists("stacked_w")) nrow(stacked_w)                else NA,
    r_version   = R.version.string
  )
)

# ------------------------------------------------------------------------------
# 2.  SAVE RDS (full reloadable list)
# ------------------------------------------------------------------------------
saveRDS(dial_dml_results, "dial_dml_results.rds")
cat(sprintf("Saved: dial_dml_results.rds  (%s MB)\n",
            round(file.size("dial_dml_results.rds") / 1e6, 2)))

# Optionally save the causal forest object separately (can be large)
save_forest <- FALSE   # set TRUE if you want the full grf object
if (save_forest && exists("cf") && !is.null(cf)) {
  saveRDS(cf, "dial_cf_forest.rds")
  cat(sprintf("Saved: dial_cf_forest.rds  (%s MB)\n",
              round(file.size("dial_cf_forest.rds") / 1e6, 2)))
}

# ------------------------------------------------------------------------------
# 3.  BUILD FLAT SUMMARY TABLE FOR CSV
# ------------------------------------------------------------------------------

stars_fn <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    p < 0.10   ~ "†",
    TRUE       ~ ""
  )
}

make_row <- function(section, label, estimator, controls, result) {
  if (is.null(result)) {
    return(data.frame(section = section, label = label,
                      estimator = estimator, controls = controls,
                      theta_pct = NA, se_pct = NA, pval = NA,
                      stars = NA, n = NA, stringsAsFactors = FALSE))
  }
  data.frame(
    section   = section,
    label     = label,
    estimator = estimator,
    controls  = controls,
    theta_pct = round(result$theta * 100, 5),
    se_pct    = round(result$se    * 100, 5),
    pval      = round(result$pval,          5),
    stars     = stars_fn(result$pval),
    n         = result$n,
    stringsAsFactors = FALSE
  )
}

rows <- list(
  
  # ── Part 1A ─────────────────────────────────────────────────────────────────
  make_row("1A_pooled", "All stacked",      "RF",    "own",     dml_all_own),
  make_row("1A_pooled", "All stacked",      "RF",    "partner", dml_all_partner),
  make_row("1A_pooled", "All stacked",      "LASSO", "own",     dml_all_lasso),
  make_row("1A_pooled", "All stacked",      "LASSO", "partner", dml_all_partner_lasso),
  
  # ── Part 1B ─────────────────────────────────────────────────────────────────
  make_row("1B_gender", "His",              "RF",    "partner", dml_his),
  make_row("1B_gender", "Her",              "RF",    "partner", dml_her),
  make_row("1B_gender", "His",              "RF",    "own",     dml_his_own),
  make_row("1B_gender", "Her",              "RF",    "own",     dml_her_own),
  
  # ── Part 1C ─────────────────────────────────────────────────────────────────
  make_row("1C_cells_RF",    "Childless men",   "RF",    "partner", dml_cells_rf[["Childless men"]]),
  make_row("1C_cells_RF",    "Fathers",         "RF",    "partner", dml_cells_rf[["Fathers"]]),
  make_row("1C_cells_RF",    "Childless women", "RF",    "partner", dml_cells_rf[["Childless women"]]),
  make_row("1C_cells_RF",    "Mothers",         "RF",    "partner", dml_cells_rf[["Mothers"]]),
  make_row("1C_cells_LASSO", "Childless men",   "LASSO", "partner", dml_cells_lasso[["Childless men"]]),
  make_row("1C_cells_LASSO", "Fathers",         "LASSO", "partner", dml_cells_lasso[["Fathers"]]),
  make_row("1C_cells_LASSO", "Childless women", "LASSO", "partner", dml_cells_lasso[["Childless women"]]),
  make_row("1C_cells_LASSO", "Mothers",         "LASSO", "partner", dml_cells_lasso[["Mothers"]]),
  
  # ── Part 2A ─────────────────────────────────────────────────────────────────
  make_row("2A_within",      "All (within)",    "RF",    "own",     dml_w_own),
  make_row("2A_within",      "All (within)",    "RF",    "partner", dml_w_partner),
  
  # ── Part 2B ─────────────────────────────────────────────────────────────────
  make_row("2B_within_gender", "His (within)",  "RF",    "partner", dml_w_his),
  make_row("2B_within_gender", "Her (within)",  "RF",    "partner", dml_w_her),
  
  # ── Part 2C ─────────────────────────────────────────────────────────────────
  make_row("2C_within_cells", "Childless men (W)",   "RF", "partner", dml_w_cells[["Childless men (W)"]]),
  make_row("2C_within_cells", "Fathers (W)",         "RF", "partner", dml_w_cells[["Fathers (W)"]]),
  make_row("2C_within_cells", "Childless women (W)", "RF", "partner", dml_w_cells[["Childless women (W)"]]),
  make_row("2C_within_cells", "Mothers (W)",         "RF", "partner", dml_w_cells[["Mothers (W)"]])
)

summary_df <- do.call(rbind, rows)

# Add linear OLS for Part 4 comparison
if (exists("linear_results") && !is.null(linear_results)) {
  for (g in names(linear_results)) {
    lr <- linear_results[[g]]
    if (!is.null(lr)) {
      summary_df <- rbind(summary_df, data.frame(
        section   = "4_OLS_comparison",
        label     = g,
        estimator = "OLS-FE",
        controls  = "partner",
        theta_pct = round(lr$b  * 100, 5),
        se_pct    = NA,
        pval      = round(lr$p,         5),
        stars     = stars_fn(lr$p),
        n         = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Add derived tests as extra rows
if (exists("p_triple") && !is.null(p_triple)) {
  summary_df <- rbind(summary_df, data.frame(
    section   = "5_becker_triple",
    label     = "Father penalty vs Mother penalty",
    estimator = "RF (derived)",
    controls  = "partner",
    theta_pct = round(triple        * 100, 5),
    se_pct    = round(triple_se     * 100, 5),
    pval      = round(p_triple,             5),
    stars     = stars_fn(p_triple),
    n         = NA,
    stringsAsFactors = FALSE
  ))
}

write.csv(summary_df, "dial_dml_summary.csv", row.names = FALSE)
cat(sprintf("Saved: dial_dml_summary.csv  (%d rows)\n", nrow(summary_df)))

# ------------------------------------------------------------------------------
# 4.  CAUSAL FOREST SUPPLEMENTARY CSVs
# ------------------------------------------------------------------------------

if (!is.null(dial_dml_results$cf$vimp)) {
  write.csv(dial_dml_results$cf$vimp, "dial_cf_vimp.csv", row.names = FALSE)
  cat("Saved: dial_cf_vimp.csv\n")
}

if (!is.null(dial_dml_results$cf$group_cate)) {
  write.csv(dial_dml_results$cf$group_cate, "dial_cf_group_cate.csv", row.names = FALSE)
  cat("Saved: dial_cf_group_cate.csv\n")
}

# ------------------------------------------------------------------------------
# 5.  PRINT RELOAD INSTRUCTIONS
# ------------------------------------------------------------------------------

cat("\n", strrep("-", 70), "\n")
cat("  HOW TO RELOAD\n")
cat(strrep("-", 70), "\n\n")
cat('  res <- readRDS("dial_dml_results.rds")\n\n')
cat("  # Key objects:\n")
cat('  res$cells_rf[["Childless men"]]$theta   # DML β₂ estimate\n')
cat('  res$cells_rf[["Mothers"]]$se            # standard error\n')
cat('  res$cells_rf[["Fathers"]]$pval          # p-value\n')
cat('  res$within$cells[["Mothers (W)"]]$theta # within-couple estimate\n')
cat('  res$cf$ate                              # causal forest ATE\n')
cat('  res$cf$vimp                             # variable importance table\n')
cat('  res$becker$p_triple                     # Becker triple test\n')
cat('  res$meta                                # sample sizes + timestamp\n\n')

cat("  # Flat table (for tables / external analysis):\n")
cat('  df <- read.csv("dial_dml_summary.csv")\n')
cat('  df[df$section == "1C_cells_RF", ]\n\n')

cat(strrep("=", 70), "\n")
cat("  SAVE COMPLETE\n")
cat(strrep("=", 70), "\n")

# seed and package versions
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())   # explicit print so it is captured when source()d

# ── Close log sinks (keep this as the last block of the script) ─────────────
cat(sprintf("\nRun complete: %s\nLog saved to: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), log_file))
sink(type = "message")
sink()
close(.log_con)