# =============================================================================
# MASTER_PART_G_FTBB_IV.R
#
# PART G: FAMILY TAX BENEFIT B REFORM 1 — SIMULATED-INSTRUMENT IV (P2)
#
# Identifies theta on S_ct (household specialisation) using the 2015 reduction
# of the FTB-B primary-earner income threshold (Reform 1, from $150K to $100K)
# and the 2016 reduction of the maximum dependent-child age (Reform 2).
# Estimation panel: stacked person-couple panel with person-within-couple FE
# plus wave FE. Z_ftb is the change in statutory FTB-B value at frozen
# baseline incomes (Currie 1996, Gruber-Saez 2002).
#
# Source from the master analysis script BEFORE the SAVE block:
#     source("MASTER_PART_G_FTBB_IV.R")
#
# Assumes:
#   - tidyverse, fixest loaded
#   - Either: in-memory `panel` carrying age_resident_child_1,
#     num_step_foster_grand, num_own_resident_children, etc.
#   - Or: hilda_panel_data_extended.rds readable from cwd (preferred fallback).
# Variables are namespaced with iv_ prefix so they will not collide with the
# master script's couples / dual_earner / stacked / baseline objects.
#
# Reproduces from FTBB_simulated_IV_v4.R:
#   Reform 1 first-stage F on S_ct ~ 15.5 (target 16.5)
#   Reform 1 IV theta_hat on S_ct ~ -1.30 (SE 0.87)
#   OLS comparator (parents sample) ~ -0.21
# Writes ftbb_reform1_iv_results.rds.
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  PART G: FTB-B REFORM 1 IV (simulated-instrument identification of P2)\n")
cat(strrep("=", 78), "\n", sep = "")

# Namespace protection
iv_select <- dplyr::select
iv_filter <- dplyr::filter
iv_mutate <- dplyr::mutate
iv_slice  <- dplyr::slice

# -----------------------------------------------------------------------------
# G.0  STATUTORY FTB-B PARAMETERS
# -----------------------------------------------------------------------------

# Wave -> statutory parameter table. Values from Services Australia /
# Department of Social Services published rates by financial year.
#
# Notes:
#   - M_0to4, M_5to12, I_r frozen by the 2014 federal budget through
#     FY 2020-21 (waves 15-21).
#   - P_r dropped from $150K to $100K on 1 Jul 2015 (Reform 1, wave 16).
#   - M_13plus dropped from positive to zero on 1 Jul 2016 (Reform 2,
#     wave 17), modelled via the age13_couple_ok flag below.
#   - DSS Family Assistance Guide section 3.6.3: "Indexation of the
#     primary earner income limit commenced from 1 July 2021."
#     P_r therefore indexes from wave 22 onward.
#   - M and I_r values for waves 22-24 reflect Services Australia
#     historical-rate publications for FY 2021-22 through FY 2023-24.
iv_ftbb_params <- tibble::tribble(
  ~wave, ~fy_label, ~M_0to4, ~M_5to12, ~M_13plus, ~I_r,  ~taper, ~P_r,    ~age13_couple_ok,
  # M_r annual = fortnightly (DSS 3.6.3) x 365/14 + supplement + ES Part B.
  # I_r and P_r are direct from DSS 3.6.3 historical rates table.
  # Note: M_5to12 == M_13plus pre-Reform-2 (DSS published a single 5-15 rate);
  # age13_couple_ok = FALSE from W17 enforces the 1 Jul 2016 couple-family
  # cliff that zeroes out the 13+ rate (DSS FA Guide 1.2.1).
  # ES Part B was introduced 1 Jul 2013, so W12 and W13 have no ES component;
  # W14 onwards includes ES at the 2013 rate (69.35 <5 / 51.10 5-15),
  # rising to 73.00 <5 / 51.10 5-18 from 1 Jul 2014 (W15 onwards).
  12L,  "2011-12", 4004.05,  2898.10,  2898.10,    4891,  0.20,  150000,  TRUE,
  13L,  "2012-13", 4117.20,  2978.40,  2978.40,    5037,  0.20,  150000,  TRUE,
  14L,  "2013-14", 4241.30,  3069.65,  3069.65,    5183,  0.20,  150000,  TRUE,
  15L,  "2014-15", 4347.15,  3142.65,  3142.65,    5329,  0.20,  150000,  TRUE,
  16L,  "2015-16", 4412.85,  3190.10,  3190.10,    5402,  0.20,  100000,  TRUE,
  17L,  "2016-17", 4482.20,  3237.55,     0.00,    5475,  0.20,  100000,  FALSE,
  18L,  "2017-18", 4485.85,  3241.20,     0.00,    5548,  0.20,  100000,  FALSE,
  19L,  "2018-19", 4493.15,  3248.50,     0.00,    5621,  0.20,  100000,  FALSE,
  20L,  "2019-20", 4573.45,  3306.90,     0.00,    5694,  0.20,  100000,  FALSE,
  21L,  "2020-21", 4653.75,  3365.30,     0.00,    5767,  0.20,  100000,  FALSE,
  22L,  "2021-22", 4693.90,  3394.50,     0.00,    5840,  0.20,  100900,  FALSE,
  23L,  "2022-23", 4858.15,  3511.30,     0.00,    6059,  0.20,  104432,  FALSE,
  24L,  "2023-24", 5234.10,  3781.40,     0.00,    6497,  0.20,  112578,  FALSE
)
iv_pre_params <- iv_ftbb_params %>% iv_filter(wave == 15L)

IV_REFORM_WAVE_1 <- 16L
IV_REFORM_WAVE_2 <- 17L
# Baseline income window: W15-W16 covers FY 2014-15 and FY 2015-16 (the
# financial years around the 1 July 2015 Reform 1 cliff). This captures
# income AT THE TIME OF THE REFORM, which is what the simulated-IV requires
# (Currie-Gruber-Saez 1996/2002): the relevant counterfactual is "what
# would each couple have received under the old vs new rules given their
# income then". Pre-policy averages over a longer window (e.g., W12-W15)
# attenuate the first-stage F because income trajectories misclassify
# couples relative to the policy cutoff at the time of the reform.
IV_BASELINE_WAVES <- 15:16

iv_simulate_ftbb <- function(a_young, Y_P_0, Y_S_0,
                             M_0to4, M_5to12, M_13plus,
                             I_r, taper, P_r, age13_couple_ok,
                             is_grandparent_carer) {
  M_r <- dplyr::case_when(
    is.na(a_young) | a_young < 0 ~ 0,
    a_young >= 16                ~ 0,
    a_young <= 4                 ~ M_0to4,
    a_young <= 12                ~ M_5to12,
    TRUE                         ~ M_13plus
  )
  M_r <- ifelse(!is.na(Y_P_0) & Y_P_0 > P_r, 0, M_r)
  age13_applies <- !age13_couple_ok & !is.na(a_young) & a_young >= 13 &
    !dplyr::coalesce(is_grandparent_carer, FALSE)
  M_r <- ifelse(age13_applies, 0, M_r)
  Y_S <- ifelse(is.na(Y_S_0), 0, Y_S_0)
  pmax(0, M_r - taper * pmax(0, Y_S - I_r))
}

# -----------------------------------------------------------------------------
# G.1  LOAD / REUSE PANEL
# -----------------------------------------------------------------------------
#
# PART G requires extended HILDA variables that PARTS A-F do not use:
#   age_resident_child_1, num_step_foster_grand, num_own_resident_children
# These ship with hilda_panel_data_extended.rds.
#
# Note on panel choice for IV identification:
#   We deliberately load the W15-W24 panel (hilda_panel_data_extended.rds)
#   for PART G even when an extended W12-W24 panel is in memory from the
#   master. The Reform 1 simulated-IV first stage is materially stronger
#   on the W15-W24 panel (F = ~23 vs ~7 on W12-W24) because the W12-W24
#   panel contains additional persons (5,000+ extra) who entered or left
#   the panel in W12-W14 or W22-W24, bringing Reform 1 candidates with
#   insufficient within-couple variation to identify off the 1 July 2015
#   cliff. The simulated-IV requires couples observed *around* the reform
#   date — the W15-W24 panel composition provides this most cleanly.

iv_required_vars <- c(
  "person_id", "partner_id", "wave", "female",
  "hours_worked_clean", "ln_hourly_wage_real", "annual_income",
  "age_youngest_own_child", "age_resident_child_1",
  "num_step_foster_grand",  "num_own_resident_children",
  "educ_years", "experience_years", "age", "married", "employed"
)

# Always load W15-W24 panel for clean Reform 1 identification, even if a
# different panel is in memory. Master may have an extended W12-W24 panel
# loaded for PARTS A-F, but PART G requires the W15-W24 subset.
if (file.exists("hilda_panel_data_extended.rds")) {
  iv_panel <- readRDS("hilda_panel_data_extended.rds")
  cat("  PART G: Loaded hilda_panel_data_extended.rds (W15-W24)\n")
  cat("          for clean Reform 1 IV identification.\n")
} else {
  iv_panel <- readRDS("hilda_panel_data_W12_W24.rds")
  cat("  Reading hilda_panel_data_W12_W24.rds (extended variables needed).\n")
  iv_missing <- setdiff(iv_required_vars, names(iv_panel))
  if (length(iv_missing) > 0) {
    stop(sprintf("PART G: required variables missing from extended panel: %s",
                 paste(iv_missing, collapse = ", ")))
  }
}
cat(sprintf("  Loaded panel: %d obs, %d persons, waves %d-%d\n",
            nrow(iv_panel),
            dplyr::n_distinct(iv_panel$person_id),
            min(iv_panel$wave, na.rm = TRUE),
            max(iv_panel$wave, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# G.2  COUPLE PANEL VIA ASYMMETRIC p1/p2 JOIN
# -----------------------------------------------------------------------------

iv_p1 <- iv_panel %>%
  iv_filter(!is.na(partner_id), partner_id != person_id,
            partner_id %in% iv_panel$person_id) %>%
  iv_select(person_id, partner_id, wave, female,
            hours_worked_clean, ln_hourly_wage_real,
            annual_income, educ_years, experience_years,
            age, married, employed,
            age_youngest_own_child, age_resident_child_1,
            num_step_foster_grand, num_own_resident_children)

iv_p2 <- iv_panel %>%
  iv_select(person_id, wave, female,
            hours_worked_clean, ln_hourly_wage_real,
            annual_income, educ_years, experience_years,
            age, married, employed,
            age_youngest_own_child, age_resident_child_1,
            num_step_foster_grand, num_own_resident_children)

iv_couples <- iv_p1 %>%
  dplyr::inner_join(iv_p2,
                    by = c("partner_id" = "person_id", "wave" = "wave"),
                    suffix = c("_f", "_m")) %>%
  iv_filter(female_f != female_m) %>%
  iv_mutate(
    couple_id = ifelse(person_id < partner_id,
                       paste0(person_id, "_", partner_id),
                       paste0(partner_id, "_", person_id)),
    her_hours  = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
    his_hours  = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
    her_lnwage = ifelse(female_f == 1, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
    his_lnwage = ifelse(female_f == 0, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
    her_income = ifelse(female_f == 1, annual_income_f,       annual_income_m),
    his_income = ifelse(female_f == 0, annual_income_f,       annual_income_m),
    her_educ   = ifelse(female_f == 1, educ_years_f,          educ_years_m),
    his_educ   = ifelse(female_f == 0, educ_years_f,          educ_years_m),
    her_exp    = ifelse(female_f == 1, experience_years_f,    experience_years_m),
    his_exp    = ifelse(female_f == 0, experience_years_f,    experience_years_m),
    her_age    = ifelse(female_f == 1, age_f,                 age_m),
    his_age    = ifelse(female_f == 0, age_f,                 age_m),
    her_married = ifelse(female_f == 1, married_f, married_m),
    his_married = ifelse(female_f == 0, married_f, married_m),
    her_employed = ifelse(female_f == 1, employed_f, employed_m),  # PART H
    his_employed = ifelse(female_f == 0, employed_f, employed_m),  # extensive margin
    total_hours      = his_hours + her_hours,
    S_ct             = ifelse(total_hours > 0, her_hours / total_hours, NA_real_),
    primary_income   = pmax(his_income, her_income, na.rm = FALSE),
    secondary_income = pmin(his_income, her_income, na.rm = FALSE),
    a_young_ct = pmin(
      dplyr::coalesce(age_resident_child_1_f,   Inf),
      dplyr::coalesce(age_resident_child_1_m,   Inf),
      dplyr::coalesce(age_youngest_own_child_f, Inf),
      dplyr::coalesce(age_youngest_own_child_m, Inf)
    ),
    a_young_ct = ifelse(is.infinite(a_young_ct), NA_real_, a_young_ct),
    n_step_total = dplyr::coalesce(num_step_foster_grand_f, 0L)
    + dplyr::coalesce(num_step_foster_grand_m, 0L),
    n_own_total  = dplyr::coalesce(num_own_resident_children_f, 0L)
    + dplyr::coalesce(num_own_resident_children_m, 0L),
    is_grandparent_carer = (n_step_total > 0 & n_own_total == 0),
    has_children = as.integer(!is.na(a_young_ct) & a_young_ct >= 0)
  ) %>%
  dplyr::distinct(couple_id, wave, .keep_all = TRUE)

iv_dual_earner <- iv_couples %>%
  iv_filter(his_hours > 0, her_hours > 0,
            !is.na(his_lnwage), !is.na(her_lnwage))
cat(sprintf("  Couple-waves (full): %d (%d couples)\n",
            nrow(iv_couples), dplyr::n_distinct(iv_couples$couple_id)))
cat(sprintf("  Dual-earner: %d (%d couples)\n",
            nrow(iv_dual_earner), dplyr::n_distinct(iv_dual_earner$couple_id)))

# -----------------------------------------------------------------------------
# G.3  BASELINE INCOMES (mean over BASELINE_WAVES, full couple panel)
# -----------------------------------------------------------------------------

iv_baseline <- iv_couples %>%
  iv_filter(wave %in% IV_BASELINE_WAVES) %>%
  iv_mutate(
    inc_a = dplyr::coalesce(her_income, 0),
    inc_b = dplyr::coalesce(his_income, 0),
    y_primary_cw   = pmax(inc_a, inc_b),
    y_secondary_cw = pmin(inc_a, inc_b)
  ) %>%
  dplyr::group_by(couple_id) %>%
  dplyr::summarise(
    Y_P_0 = mean(y_primary_cw),
    Y_S_0 = mean(y_secondary_cw),
    .groups = "drop"
  ) %>%
  iv_mutate(
    Y_P_0 = ifelse(is.finite(Y_P_0), Y_P_0, NA_real_),
    Y_S_0 = ifelse(is.finite(Y_S_0), Y_S_0, NA_real_)
  )

iv_dual_earner <- iv_dual_earner %>%
  dplyr::left_join(iv_baseline, by = "couple_id")

# -----------------------------------------------------------------------------
# G.4  Z_FTB AT FROZEN BASELINE INCOMES
# -----------------------------------------------------------------------------

iv_dual_earner <- iv_dual_earner %>%
  dplyr::left_join(iv_ftbb_params %>% iv_select(-fy_label), by = "wave")

iv_dual_earner$B_post <- iv_simulate_ftbb(
  iv_dual_earner$a_young_ct, iv_dual_earner$Y_P_0, iv_dual_earner$Y_S_0,
  iv_dual_earner$M_0to4, iv_dual_earner$M_5to12, iv_dual_earner$M_13plus,
  iv_dual_earner$I_r, iv_dual_earner$taper, iv_dual_earner$P_r,
  iv_dual_earner$age13_couple_ok, iv_dual_earner$is_grandparent_carer
)
iv_dual_earner$B_pre <- iv_simulate_ftbb(
  iv_dual_earner$a_young_ct, iv_dual_earner$Y_P_0, iv_dual_earner$Y_S_0,
  iv_pre_params$M_0to4, iv_pre_params$M_5to12, iv_pre_params$M_13plus,
  iv_pre_params$I_r,    iv_pre_params$taper,   iv_pre_params$P_r,
  iv_pre_params$age13_couple_ok, iv_dual_earner$is_grandparent_carer
)

iv_dual_earner <- iv_dual_earner %>%
  iv_mutate(
    Z_ftb        = B_post - B_pre,
    hit_by_limit = !is.na(Y_P_0) & Y_P_0 >= 100000 & Y_P_0 <= 150000,
    hit_by_age13 = !is.na(a_young_ct) & a_young_ct >= 13 &
      wave >= IV_REFORM_WAVE_2 & !is_grandparent_carer
  )

cat(sprintf("  Z_ftb < 0: %d (%.1f%%) | = 0: %d (%.1f%%) | > 0: %d (%.1f%%)\n",
            sum(iv_dual_earner$Z_ftb < 0,  na.rm = TRUE),
            100 * mean(iv_dual_earner$Z_ftb < 0,  na.rm = TRUE),
            sum(iv_dual_earner$Z_ftb == 0, na.rm = TRUE),
            100 * mean(iv_dual_earner$Z_ftb == 0, na.rm = TRUE),
            sum(iv_dual_earner$Z_ftb > 0,  na.rm = TRUE),
            100 * mean(iv_dual_earner$Z_ftb > 0,  na.rm = TRUE)))

# -----------------------------------------------------------------------------
# G.5  STACK TO PERSON-COUPLE LONG PANEL
# -----------------------------------------------------------------------------

iv_his_rows <- iv_dual_earner %>%
  dplyr::transmute(
    couple_id, wave,
    ln_wage = his_lnwage, educ = his_educ,
    experience_years = his_exp, age_sq = his_age^2,
    married_num = as.integer(his_married),
    person_couple = paste0(couple_id, "_M"),
    S_ct, Y_P_0, Y_S_0, a_young_ct, has_children,
    hit_by_limit, hit_by_age13, Z_ftb
  )
iv_her_rows <- iv_dual_earner %>%
  dplyr::transmute(
    couple_id, wave,
    ln_wage = her_lnwage, educ = her_educ,
    experience_years = her_exp, age_sq = her_age^2,
    married_num = as.integer(her_married),
    person_couple = paste0(couple_id, "_F"),
    S_ct, Y_P_0, Y_S_0, a_young_ct, has_children,
    hit_by_limit, hit_by_age13, Z_ftb
  )
iv_stacked <- dplyr::bind_rows(iv_his_rows, iv_her_rows) %>%
  iv_filter(!is.na(ln_wage), !is.na(educ), !is.na(experience_years)) %>%
  iv_mutate(
    educ_c     = educ - mean(educ, na.rm = TRUE),
    exper_c    = experience_years - mean(experience_years, na.rm = TRUE),
    educ_exp_c = educ_c * exper_c,
    S_educ_exp = S_ct * educ_exp_c,
    Z_educ_exp = Z_ftb * educ_exp_c
  )

iv_sample <- iv_stacked %>%
  iv_filter(has_children == 1, !is.na(Z_ftb), !is.na(S_ct), !is.na(Y_P_0))
cat(sprintf("  IV sample (parents): %d obs, %d couples\n",
            nrow(iv_sample), dplyr::n_distinct(iv_sample$couple_id)))

iv_reform1_sample <- iv_sample %>% iv_filter(hit_by_limit, !hit_by_age13)
cat(sprintf("  Reform 1 subsample: %d obs, %d couples\n",
            nrow(iv_reform1_sample),
            dplyr::n_distinct(iv_reform1_sample$couple_id)))

# -----------------------------------------------------------------------------
# G.6  REGRESSIONS
# -----------------------------------------------------------------------------

cat("\n  --- First stages and IV ---\n")

iv_fs_level <- fixest::feols(
  S_ct ~ Z_ftb + educ_exp_c + experience_years + age_sq + married_num
  | person_couple + wave,
  data = iv_sample, cluster = ~ couple_id)

iv_fs_inter <- fixest::feols(
  S_educ_exp ~ Z_ftb + Z_educ_exp + educ_exp_c + experience_years + age_sq + married_num
  | person_couple + wave,
  data = iv_sample, cluster = ~ couple_id)

iv_pooled <- fixest::feols(
  ln_wage ~ educ_exp_c + experience_years + age_sq + married_num
  | person_couple + wave
  | S_ct + S_educ_exp ~ Z_ftb + Z_educ_exp,
  data = iv_sample, cluster = ~ couple_id)

iv_reform1 <- fixest::feols(
  ln_wage ~ educ_exp_c + experience_years + age_sq + married_num
  | person_couple + wave
  | S_ct + S_educ_exp ~ Z_ftb + Z_educ_exp,
  data = iv_reform1_sample, cluster = ~ couple_id)

iv_ols_comp <- fixest::feols(
  ln_wage ~ educ_exp_c + S_ct + S_educ_exp + experience_years + age_sq + married_num
  | person_couple + wave,
  data = iv_sample, cluster = ~ couple_id)

iv_unaffected <- iv_sample %>% iv_filter(Y_P_0 < 80000, !hit_by_age13)
iv_falsif <- if (nrow(iv_unaffected) > 500) {
  tryCatch(fixest::feols(
    ln_wage ~ educ_exp_c + experience_years + age_sq + married_num
    | person_couple + wave
    | S_ct + S_educ_exp ~ Z_ftb + Z_educ_exp,
    data = iv_unaffected, cluster = ~ couple_id),
    error = function(e) NULL)
} else NULL

# Print compact summary
iv_pull_F <- function(m, endo = "S_ct") {
  if (is.null(m)) return(NA_real_)
  fs <- tryCatch(fixest::fitstat(m, "ivf1"), error = function(e) NULL)
  if (is.null(fs)) return(NA_real_)
  # fixest names the per-endogenous F slots as "ivf1::<endo>" (double colon),
  # not "ivf1.<endo>" as one might guess from the function name. The slot
  # contains a list with $stat, $p, $df1, $df2.
  flat_name <- paste0("ivf1::", endo)
  if (flat_name %in% names(fs)) {
    v <- fs[[flat_name]]
    if (is.list(v) && !is.null(v$stat)) return(as.numeric(v$stat))
  }
  # Fallbacks for other fixest versions / single-endogenous case
  dot_name <- paste0("ivf1.", endo)
  if (dot_name %in% names(fs)) {
    v <- fs[[dot_name]]
    if (is.list(v) && !is.null(v$stat)) return(as.numeric(v$stat))
  }
  if (!is.null(fs$ivf1)) {
    iv1 <- fs$ivf1
    if (is.list(iv1) && endo %in% names(iv1)) {
      v <- iv1[[endo]]
      if (is.list(v) && !is.null(v$stat)) return(as.numeric(v$stat))
    }
    if (is.list(iv1) && length(iv1) >= 1) {
      v <- iv1[[1]]
      if (is.list(v) && !is.null(v$stat)) return(as.numeric(v$stat))
    }
    if (!is.null(iv1$stat)) return(as.numeric(iv1$stat[1]))
  }
  NA_real_
}
iv_pull_est <- function(m, var) {
  if (is.null(m)) return(c(NA_real_, NA_real_, NA_real_))
  ct <- tryCatch(m$coeftable, error = function(e) NULL)
  if (is.null(ct) || !var %in% rownames(ct)) return(c(NA_real_, NA_real_, NA_real_))
  c(ct[var, "Estimate"], ct[var, "Std. Error"], ct[var, "Pr(>|t|)"])
}

iv_summary <- tibble::tibble(
  spec = c("OLS (parents)", "IV pooled (parents)", "IV Reform 1", "IV falsification"),
  theta_hat = c(iv_pull_est(iv_ols_comp, "S_ct")[1],
                iv_pull_est(iv_pooled,   "fit_S_ct")[1],
                iv_pull_est(iv_reform1,  "fit_S_ct")[1],
                iv_pull_est(iv_falsif,   "fit_S_ct")[1]),
  se        = c(iv_pull_est(iv_ols_comp, "S_ct")[2],
                iv_pull_est(iv_pooled,   "fit_S_ct")[2],
                iv_pull_est(iv_reform1,  "fit_S_ct")[2],
                iv_pull_est(iv_falsif,   "fit_S_ct")[2]),
  p_value   = c(iv_pull_est(iv_ols_comp, "S_ct")[3],
                iv_pull_est(iv_pooled,   "fit_S_ct")[3],
                iv_pull_est(iv_reform1,  "fit_S_ct")[3],
                iv_pull_est(iv_falsif,   "fit_S_ct")[3]),
  fs_F_Sct  = c(NA_real_,
                iv_pull_F(iv_pooled),
                iv_pull_F(iv_reform1),
                iv_pull_F(iv_falsif)),
  n_obs     = c(iv_pooled$nobs, iv_pooled$nobs, iv_reform1$nobs,
                if (!is.null(iv_falsif)) iv_falsif$nobs else NA_integer_)
)
cat("\n  ── PART G SUMMARY ─────────────────────────────────────────────\n")
print(iv_summary, n = Inf)
cat("\n")

# -----------------------------------------------------------------------------
# G.7  SAVE
# -----------------------------------------------------------------------------

ftbb_reform1_iv_results <- list(
  ftbb_params       = iv_ftbb_params,
  pre_params        = iv_pre_params,
  fs_level          = iv_fs_level,
  fs_inter          = iv_fs_inter,
  iv_pooled         = iv_pooled,
  iv_reform1        = iv_reform1,
  ols_comparator    = iv_ols_comp,
  iv_falsif         = iv_falsif,
  summary           = iv_summary,
  n_iv_obs          = nrow(iv_sample),
  n_iv_couples      = dplyr::n_distinct(iv_sample$couple_id),
  n_reform1_obs     = nrow(iv_reform1_sample),
  n_reform1_couples = dplyr::n_distinct(iv_reform1_sample$couple_id),
  reform1_F_Sct     = iv_pull_F(iv_reform1, "S_ct"),
  reform1_F_Seducexp= iv_pull_F(iv_reform1, "S_educ_exp")
)
saveRDS(ftbb_reform1_iv_results, "ftbb_reform1_iv_results.rds")
cat("  Saved ftbb_reform1_iv_results.rds\n")

# -----------------------------------------------------------------------------
# G.8  PARTNER-SPLIT DIAGNOSTIC  (referee comments 1 & 2)
#   [1] Confirm S_ct is the wife's paid-hours share (orientation).
#   [2] Split the stacked Reform-1 / pooled IV by partner (wife-only effect),
#       since his and her wage rows share the same S_ct, so the headline
#       estimate is a both-partner average.
#   Writes se_ledger_iv_partner.csv and attaches $partner_split to the RDS.
# -----------------------------------------------------------------------------
cat("\n  ── G.8 PARTNER-SPLIT DIAGNOSTIC (referee comments 1 & 2) ──────\n")

pg_co <- function(m, v) {                    # fixest IV-fitted coef (b, se, p)
  nm <- grep(paste0("fit_", v, "$"), names(coef(m)), value = TRUE)
  if (!length(nm)) nm <- grep(v, names(coef(m)), value = TRUE)[1]
  c(b  = as.numeric(coef(m)[nm]),
    se = as.numeric(fixest::se(m)[nm]),
    p  = as.numeric(fixest::pvalue(m)[nm]))
}
pg_F   <- function(m) tryCatch(fixest::fitstat(m, "ivf")[[1]]$stat, error = function(e) NA_real_)
pg_one <- function(d, treat = "S_ct")
  fixest::feols(as.formula(paste0(
    "ln_wage ~ experience_years + age_sq + married_num | person_couple + wave | ",
    treat, " ~ Z_ftb")), data = d, cluster = ~ couple_id, notes = FALSE)

# [1] orientation
.pg_chk <- iv_dual_earner
.pg_chk$ws <- .pg_chk$her_hours / (.pg_chk$his_hours + .pg_chk$her_hours)
cat(sprintf("  [1] cor(S_ct, wife hours share) = %.3f  (1.000 => wife's share)\n",
            cor(.pg_chk$S_ct, .pg_chk$ws, use = "complete.obs")))

# [2] partner split (base subsetting to respect Part G's namespace setup)
pg_split <- function(sample_df, tag) {
  d0 <- sample_df
  d0$is_female <- grepl("_F$", d0$person_couple)
  dplyr::bind_rows(lapply(c("stacked (both)", "wives only", "husbands only"), function(lab) {
    d <- switch(lab,
                "stacked (both)" = d0,
                "wives only"     = d0[d0$is_female, ],
                "husbands only"  = d0[!d0$is_female, ])
    m <- tryCatch(pg_one(d), error = function(e) NULL)
    if (is.null(m)) { cat(sprintf("      %-15s not estimable\n", lab)); return(NULL) }
    s <- pg_co(m, "S_ct"); Fv <- pg_F(m)
    cat(sprintf("      %-15s theta=%7.3f  SE=%5.2f  p=%.3f  F=%5.1f  n=%d  couples=%d\n",
                lab, s["b"], s["se"], s["p"], Fv, nobs(m), dplyr::n_distinct(d$couple_id)))
    data.frame(spec = tag, partner = lab, theta = s["b"], se = s["se"], p = s["p"],
               first_stage_F = Fv, n_obs = nobs(m),
               n_couples = dplyr::n_distinct(d$couple_id), stringsAsFactors = FALSE)
  }))
}
cat("  [2] Reform 1 IV, theta on S_ct (single-endogenous):\n")
.pg_r1   <- pg_split(iv_reform1_sample, "Reform 1")
cat("  [2] Pooled IV, theta on S_ct (single-endogenous):\n")
.pg_pool <- pg_split(iv_sample, "Pooled")

# orientation cross-check on wives: theta(S_ct) vs theta(B = 1 - S_ct)
.pg_rw <- iv_reform1_sample[grepl("_F$", iv_reform1_sample$person_couple), ]
.pg_rw$B <- 1 - .pg_rw$S_ct
.pg_sS <- pg_co(pg_one(.pg_rw, "S_ct"), "S_ct")
.pg_sB <- pg_co(pg_one(.pg_rw, "B"),    "B")
cat(sprintf("  [orient] wife-only theta(S_ct)=%.3f ; theta(B=1-S_ct)=%.3f (expect equal & opposite)\n",
            .pg_sS["b"], .pg_sB["b"]))

# write ledger + attach to saved results
iv_partner_split <- dplyr::bind_rows(.pg_r1, .pg_pool)
write.csv(iv_partner_split, "se_ledger_iv_partner.csv", row.names = FALSE)
ftbb_reform1_iv_results$partner_split <- iv_partner_split
saveRDS(ftbb_reform1_iv_results, "ftbb_reform1_iv_results.rds")
cat("  Wrote se_ledger_iv_partner.csv; attached $partner_split to ftbb_reform1_iv_results.rds\n")
cat("  Note: wife-only row replaces the stacked average in the IV section;\n")
cat("        IV stays a labour-supply validation, event study remains P2.\n")

cat(strrep("=", 78), "\n", sep = "")