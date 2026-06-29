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

# =============================================================================
# PART G ROBUSTNESS BATTERY  (Appendix C: G.9-G.17)
# Merged in from the former 04b_PART_G_IV_ROBUSTNESS.R. It reuses the objects
# built in G.0-G.8 above (iv_dual_earner, iv_baseline, iv_ftbb_params,
# iv_pre_params, iv_simulate_ftbb, iv_pull_F), which are live in this session,
# so the stand-alone bootstrap that reloaded them has been removed. The bare
# dplyr verbs below need the namespace guard the estimator's iv_-prefixed
# verbs do not provide, so it is (re)set here.
# =============================================================================

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  PART G ROBUSTNESS BATTERY (Appendix C: G.9-G.17)\n")
cat(strrep("=", 78), "\n", sep = "")

# Guard: G.0-G.8 above must have populated Part G's frames and helpers.
rb_needed <- c("iv_dual_earner", "iv_baseline", "iv_ftbb_params",
               "iv_pre_params", "iv_simulate_ftbb", "iv_pull_F")
if (!all(vapply(rb_needed, exists, logical(1), envir = globalenv())))
  stop("Part G robustness battery: estimator objects missing; G.0-G.8 must run first.")

iv_rb <- list()  # results container

# Controls and FE shared across specs (match Part G G.6)
RB_CTRL <- c("educ_exp_c", "experience_years", "age_sq", "married_num")
RB_FE   <- "person_couple + wave"

# -----------------------------------------------------------------------------
# G.9  SEPARATE-REFORM INSTRUMENTS  Z1 (threshold), Z2 (age-13 cliff)
#
#   Z1 = b(R_R1) - b(R_pre)     1 Jul 2015 threshold drop $150k -> $100k
#   Z2 = b(R_R2) - b(R_R1)      1 Jul 2016 couple-family age-13 exclusion
#
# Both are evaluated at frozen baseline incomes (Y_P_0, Y_S_0) and the
# wave-varying youngest-child age a_young_ct, identically to Part G's bundled
# Z_ftb. R_pre = W15 rules (iv_pre_params); R_R1 = W16 rules (post-threshold,
# pre-age13); R_R2 = W17 rules (post-age13). Within-couple time variation comes
# from child ageing; cross-couple variation from position in the
# (Y_P_0, Y_S_0, a_young_ct) distribution.
# -----------------------------------------------------------------------------
cat("\n  -- G.9 separate-reform instruments (Z1, Z2) --------------------\n")

rb_R1_params <- iv_ftbb_params %>% filter(wave == 16L)   # post-Reform-1
rb_R2_params <- iv_ftbb_params %>% filter(wave == 17L)   # post-Reform-2

rb_sim_with <- function(df, p) {
  iv_simulate_ftbb(
    df$a_young_ct, df$Y_P_0, df$Y_S_0,
    p$M_0to4, p$M_5to12, p$M_13plus,
    p$I_r, p$taper, p$P_r, p$age13_couple_ok,
    df$is_grandparent_carer
  )
}

iv_dual_earner <- iv_dual_earner %>%
  mutate(
    B_pre_rb = rb_sim_with(., iv_pre_params),   # R_pre (should equal Part G B_pre)
    B_R1     = rb_sim_with(., rb_R1_params),    # R_R1
    B_R2     = rb_sim_with(., rb_R2_params),    # R_R2
    Z1       = B_R1 - B_pre_rb,                 # threshold-only
    Z2       = B_R2 - B_R1                      # age-13-only
  )

cat(sprintf("    Z1 (threshold)  nonzero: %d (%.1f%%)\n",
            sum(iv_dual_earner$Z1 != 0, na.rm = TRUE),
            100 * mean(iv_dual_earner$Z1 != 0, na.rm = TRUE)))
cat(sprintf("    Z2 (age-13)     nonzero: %d (%.1f%%)\n",
            sum(iv_dual_earner$Z2 != 0, na.rm = TRUE),
            100 * mean(iv_dual_earner$Z2 != 0, na.rm = TRUE)))

# Restack to person-couple long panel carrying Z1, Z2 and their interactions,
# mirroring Part G's iv_his_rows / iv_her_rows / centring exactly.
rb_stack_rows <- function(de, who = c("M", "F")) {
  who <- match.arg(who)
  if (who == "M") {
    de %>% transmute(
      couple_id, wave, ln_wage = his_lnwage, educ = his_educ,
      experience_years = his_exp, age_sq = his_age^2,
      married_num = as.integer(his_married),
      person_couple = paste0(couple_id, "_M"),
      S_ct, Y_P_0, Y_S_0, a_young_ct, has_children,
      hit_by_limit, hit_by_age13, Z_ftb, Z1, Z2, B_pre = B_pre_rb)
  } else {
    de %>% transmute(
      couple_id, wave, ln_wage = her_lnwage, educ = her_educ,
      experience_years = her_exp, age_sq = her_age^2,
      married_num = as.integer(her_married),
      person_couple = paste0(couple_id, "_F"),
      S_ct, Y_P_0, Y_S_0, a_young_ct, has_children,
      hit_by_limit, hit_by_age13, Z_ftb, Z1, Z2, B_pre = B_pre_rb)
  }
}

rb_stacked <- dplyr::bind_rows(rb_stack_rows(iv_dual_earner, "M"),
                               rb_stack_rows(iv_dual_earner, "F")) %>%
  filter(!is.na(ln_wage), !is.na(educ), !is.na(experience_years)) %>%
  mutate(
    educ_c     = educ - mean(educ, na.rm = TRUE),
    exper_c    = experience_years - mean(experience_years, na.rm = TRUE),
    educ_exp_c = educ_c * exper_c,
    S_educ_exp = S_ct * educ_exp_c,
    Z1_educ_exp = Z1 * educ_exp_c,
    Z2_educ_exp = Z2 * educ_exp_c,
    is_female   = grepl("_F$", person_couple)
  )

rb_sample   <- rb_stacked %>%
  filter(has_children == 1, !is.na(S_ct), !is.na(Y_P_0),
         !is.na(Z1), !is.na(Z2))
rb_reform1  <- rb_sample %>% filter(hit_by_limit, !hit_by_age13)

cat(sprintf("    rb_sample (parents): %d obs, %d couples\n",
            nrow(rb_sample), dplyr::n_distinct(rb_sample$couple_id)))
cat(sprintf("    rb_reform1 band:     %d obs, %d couples\n",
            nrow(rb_reform1), dplyr::n_distinct(rb_reform1$couple_id)))

iv_rb$instruments <- list(
  Z1_nonzero_frac = mean(iv_dual_earner$Z1 != 0, na.rm = TRUE),
  Z2_nonzero_frac = mean(iv_dual_earner$Z2 != 0, na.rm = TRUE),
  cor_Z1_Z2 = suppressWarnings(cor(iv_dual_earner$Z1, iv_dual_earner$Z2,
                                   use = "complete.obs"))
)

# -----------------------------------------------------------------------------
# G.10  OVERIDENTIFIED IV + OVERID J TEST (Sargan AND cluster-robust Hansen J)
#
# Two endogenous {S_ct, S_educ_exp}; four excluded instruments
# {Z1, Z2, Z1_educ_exp, Z2_educ_exp}; overid df = 4 - 2 = 2.
# -----------------------------------------------------------------------------
cat("\n  -- G.10 overidentified IV + overid J test ----------------------\n")

rb_iv_formula <- as.formula(paste0(
  "ln_wage ~ ", paste(RB_CTRL, collapse = " + "),
  " | ", RB_FE,
  " | S_ct + S_educ_exp ~ Z1 + Z2 + Z1_educ_exp + Z2_educ_exp"))

rb_overid_pooled  <- fixest::feols(rb_iv_formula, data = rb_sample,
                                   cluster = ~ couple_id, notes = FALSE)
rb_overid_reform1 <- tryCatch(
  fixest::feols(rb_iv_formula, data = rb_reform1,
                cluster = ~ couple_id, notes = FALSE),
  error = function(e) NULL)

# Homoskedastic Sargan straight from fixest
rb_sargan <- function(m) {
  if (is.null(m)) return(c(stat = NA, p = NA, df = NA))
  s <- tryCatch(fixest::fitstat(m, "sargan"), error = function(e) NULL)
  if (is.null(s)) return(c(stat = NA, p = NA, df = NA))
  c(stat = as.numeric(s$stat), p = as.numeric(s$p), df = as.numeric(s$df))
}

# Cluster-robust two-step-GMM Hansen J on FE-demeaned data.
# Y: n-vector outcome; Wn: endogenous+exogenous regressors; Hn: instruments+
# exogenous; cl: cluster id. Returns J ~ chi^2(L-k), df = #excluded - #endog.
rb_hansenJ_cluster <- function(df, endog, excl_iv, ctrl, fe, outcome, cl) {
  vars <- unique(c(outcome, endog, ctrl, excl_iv))
  mm   <- as.matrix(df[, vars, drop = FALSE])
  dm   <- tryCatch(fixest::demean(X = mm, f = df[, strsplit(fe, " \\+ ")[[1]]]),
                   error = function(e) NULL)
  if (is.null(dm)) return(c(J = NA, df = NA, p = NA))
  Y <- dm[, outcome]
  W <- dm[, c(endog, ctrl), drop = FALSE]            # k columns
  H <- dm[, c(excl_iv, ctrl), drop = FALSE]          # L columns
  n <- nrow(W); k <- ncol(W); L <- ncol(H)
  if (L <= k) return(c(J = NA, df = 0, p = NA))
  clid <- df[[cl]]
  HHi  <- tryCatch(solve(crossprod(H)), error = function(e) NULL)
  if (is.null(HHi)) return(c(J = NA, df = NA, p = NA))
  PW <- H %*% (HHi %*% crossprod(H, W))
  b1 <- tryCatch(solve(crossprod(PW, W), crossprod(PW, Y)),
                 error = function(e) NULL)                 # 2SLS step 1
  if (is.null(b1)) return(c(J = NA, df = NA, p = NA))
  u1 <- as.vector(Y - W %*% b1)
  # Homoskedastic Sargan from the same demeaned 2SLS residuals:
  # u' P_H u / (u'u / n), with P_H = H (H'H)^{-1} H'; df = L - k.
  uPu     <- as.numeric(crossprod(crossprod(H, u1), HHi %*% crossprod(H, u1)))
  sigma2  <- sum(u1^2) / n
  sargan  <- if (sigma2 > 0) uPu / sigma2 else NA_real_
  sargan_p <- if (is.finite(sargan)) pchisq(sargan, L - k, lower.tail = FALSE) else NA_real_
  G  <- H * u1                                             # n x L moments
  # cluster-summed moment covariance
  Gc <- rowsum(G, group = clid)                            # (#clusters) x L
  S  <- crossprod(Gc)                                      # L x L
  Si <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Si)) return(c(J = NA, df = NA, p = NA, sargan = sargan, sargan_p = sargan_p))
  A  <- crossprod(H, W)                                    # L x k  (= H'W)
  cY <- crossprod(H, Y)                                    # L x 1  (= H'Y)
  b2 <- tryCatch(solve(t(A) %*% Si %*% A, t(A) %*% Si %*% cY),
                 error = function(e) NULL)                 # efficient GMM
  if (is.null(b2)) return(c(J = NA, df = NA, p = NA, sargan = sargan, sargan_p = sargan_p))
  gsum <- crossprod(H, as.vector(Y - W %*% b2))            # L x 1
  Jval <- as.numeric(t(gsum) %*% Si %*% gsum)
  dfJ  <- L - k
  c(J = Jval, df = dfJ, p = pchisq(Jval, dfJ, lower.tail = FALSE),
    sargan = sargan, sargan_p = sargan_p)
}

rb_J_pooled  <- rb_hansenJ_cluster(
  rb_sample, endog = c("S_ct", "S_educ_exp"),
  excl_iv = c("Z1", "Z2", "Z1_educ_exp", "Z2_educ_exp"),
  ctrl = RB_CTRL, fe = RB_FE, outcome = "ln_wage", cl = "couple_id")
rb_J_reform1 <- tryCatch(rb_hansenJ_cluster(
  rb_reform1, endog = c("S_ct", "S_educ_exp"),
  excl_iv = c("Z1", "Z2", "Z1_educ_exp", "Z2_educ_exp"),
  ctrl = RB_CTRL, fe = RB_FE, outcome = "ln_wage", cl = "couple_id"),
  error = function(e) c(J = NA, df = NA, p = NA))

rb_sg  <- function(jvec, m) { v <- jvec["sargan"];   if (length(v) && is.finite(v)) as.numeric(v) else rb_sargan(m)["stat"] }
rb_sgp <- function(jvec, m) { v <- jvec["sargan_p"]; if (length(v) && is.finite(v)) as.numeric(v) else rb_sargan(m)["p"] }

rb_overid_tab <- data.frame(
  sample          = c("Pooled (parents)", "Reform 1 band"),
  sargan_stat     = c(rb_sg(rb_J_pooled,  rb_overid_pooled),
                      rb_sg(rb_J_reform1, rb_overid_reform1)),
  sargan_p        = c(rb_sgp(rb_J_pooled,  rb_overid_pooled),
                      rb_sgp(rb_J_reform1, rb_overid_reform1)),
  hansenJ_cluster = c(rb_J_pooled["J"],  rb_J_reform1["J"]),
  hansenJ_df      = c(rb_J_pooled["df"], rb_J_reform1["df"]),
  hansenJ_p       = c(rb_J_pooled["p"],  rb_J_reform1["p"]),
  row.names = NULL
)
cat("    Overidentification tests (H0: instruments valid):\n")
print(rb_overid_tab, digits = 3)
iv_rb$overid <- rb_overid_tab

# -----------------------------------------------------------------------------
# G.11  PER-ENDOGENOUS FIRST-STAGE F (Sanderson-Windmeijer-style) + EFFECTIVE F
# -----------------------------------------------------------------------------
cat("\n  -- G.11 per-endogenous first-stage F + effective F -------------\n")

# fixest exposes the per-endogenous first-stage F as ivf1::<endo>; with multiple
# endogenous regressors this is the Sanderson-Windmeijer conditional F.
rb_F_tab <- data.frame(
  sample        = c("Pooled", "Reform 1"),
  F_S_ct        = c(iv_pull_F(rb_overid_pooled,  "S_ct"),
                    iv_pull_F(rb_overid_reform1, "S_ct")),
  F_S_educ_exp  = c(iv_pull_F(rb_overid_pooled,  "S_educ_exp"),
                    iv_pull_F(rb_overid_reform1, "S_educ_exp")),
  row.names = NULL
)
cat("    Sanderson-Windmeijer per-endogenous first-stage F:\n")
print(rb_F_tab, digits = 3)
iv_rb$sw_F <- rb_F_tab

# Olea-Pflueger effective F for the single-endogenous wife-only spec via ivDiag
rb_effF <- NA_real_
if (requireNamespace("ivDiag", quietly = TRUE)) {
  rb_wife <- rb_reform1 %>% filter(is_female)
  # FE-partial-out then feed numeric matrices to ivDiag for effective F + AR
  rb_dm <- tryCatch(fixest::demean(
    X = as.matrix(rb_wife[, c("ln_wage", "S_ct", RB_CTRL, "Z1", "Z2")]),
    f = rb_wife[, c("person_couple", "wave")]), error = function(e) NULL)
  if (!is.null(rb_dm)) {
    eff <- tryCatch(ivDiag::ivDiag(
      data = as.data.frame(rb_dm),
      Y = "ln_wage", D = "S_ct", Z = c("Z1", "Z2"),
      controls = RB_CTRL, cl = NULL, run_AR = TRUE),
      error = function(e) NULL)
    if (!is.null(eff)) {
      rb_effF <- tryCatch(eff$F_stat[["F.effective"]], error = function(e) NA_real_)
      iv_rb$ivDiag_wife <- eff
      cat(sprintf("    Olea-Pflueger effective F (wife-only, Reform 1): %.2f\n",
                  rb_effF))
    }
  }
} else {
  cat("    [skip] package 'ivDiag' not installed; effective F not computed.\n")
}
iv_rb$effective_F_wife <- rb_effF

# -----------------------------------------------------------------------------
# G.12  REDUCED FORM  (ln wage on the instruments directly)
# -----------------------------------------------------------------------------
cat("\n  -- G.12 reduced form -------------------------------------------\n")

rb_rf_formula <- as.formula(paste0(
  "ln_wage ~ Z1 + Z2 + ", paste(RB_CTRL, collapse = " + "), " | ", RB_FE))
rb_rf_pooled  <- fixest::feols(rb_rf_formula, data = rb_sample,
                               cluster = ~ couple_id, notes = FALSE)
rb_rf_reform1 <- tryCatch(fixest::feols(rb_rf_formula, data = rb_reform1,
                                        cluster = ~ couple_id, notes = FALSE),
                          error = function(e) NULL)

# first stage of S_ct on the same instruments (for the RF/FS ratio check)
rb_fs_formula <- as.formula(paste0(
  "S_ct ~ Z1 + Z2 + ", paste(RB_CTRL, collapse = " + "), " | ", RB_FE))
rb_fs_reform1 <- tryCatch(fixest::feols(rb_fs_formula, data = rb_reform1,
                                        cluster = ~ couple_id, notes = FALSE),
                          error = function(e) NULL)

rb_rf_get <- function(m, v) if (is.null(m) || !v %in% names(coef(m))) NA_real_ else coef(m)[v]
cat("    Reduced-form coefficients (Reform 1):\n")
cat(sprintf("      Z1: %.5f   Z2: %.5f\n",
            rb_rf_get(rb_rf_reform1, "Z1"), rb_rf_get(rb_rf_reform1, "Z2")))
if (!is.null(rb_rf_reform1) && !is.null(rb_fs_reform1)) {
  cat(sprintf("      RF/FS ratio on Z1 (indirect-least-squares check): %.3f\n",
              rb_rf_get(rb_rf_reform1, "Z1") / rb_rf_get(rb_fs_reform1, "Z1")))
}
iv_rb$reduced_form <- list(pooled = rb_rf_pooled, reform1 = rb_rf_reform1,
                           first_stage_reform1 = rb_fs_reform1)

# -----------------------------------------------------------------------------
# G.13  WEAK-IV-ROBUST INFERENCE: ANDERSON-RUBIN + CLR  (wife-only, single endog)
#
# Single endogenous S_ct, instruments {Z1, Z2}; AR and CLR sets are valid
# regardless of instrument strength. FE partialled out, then ivmodel.
# -----------------------------------------------------------------------------
cat("\n  -- G.13 Anderson-Rubin / CLR confidence sets -------------------\n")

iv_rb$weak_iv_robust <- NA
if (requireNamespace("ivmodel", quietly = TRUE)) {
  rb_wife <- rb_reform1 %>% filter(is_female)
  rb_dm <- tryCatch(fixest::demean(
    X = as.matrix(rb_wife[, c("ln_wage", "S_ct", RB_CTRL, "Z1", "Z2")]),
    f = rb_wife[, c("person_couple", "wave")]), error = function(e) NULL)
  if (!is.null(rb_dm)) {
    rb_dm <- as.data.frame(rb_dm)
    m_iv <- tryCatch(ivmodel::ivmodel(
      Y = rb_dm$ln_wage, D = rb_dm$S_ct,
      Z = rb_dm[, c("Z1", "Z2")],
      X = rb_dm[, RB_CTRL]), error = function(e) NULL)
    if (!is.null(m_iv)) {
      ar  <- tryCatch(ivmodel::AR.test(m_iv), error = function(e) NULL)
      clr <- tryCatch(ivmodel::CLR(m_iv),     error = function(e) NULL)
      if (!is.null(ar))
        cat(sprintf("    AR test:  p = %.3f ; CI = [%s]\n",
                    ar$p.value,
                    paste(round(ar$ci, 3), collapse = ", ")))
      if (!is.null(clr))
        cat(sprintf("    CLR test: p = %.3f\n", clr$p.value))
      iv_rb$weak_iv_robust <- list(AR = ar, CLR = clr)
    }
  }
} else {
  cat("    [skip] package 'ivmodel' not installed; AR/CLR not computed.\n")
}

# -----------------------------------------------------------------------------
# G.14  FALSIFICATIONS
#   (a) husband's ABSOLUTE hours on the instrument  (should be ~null)
#   (b) pre-reform placebo first stage on W12-W14    (should be ~null)
#   (c) donut-hole around the $100k threshold        (first stage should hold)
# -----------------------------------------------------------------------------
cat("\n  -- G.14 falsifications -----------------------------------------\n")

# (a) Husband's hours: FTB-B targets the secondary earner; the instrument
#     should move the wife's share, not the husband's absolute hours.
rb_husb <- iv_dual_earner %>%
  filter(!is.na(his_hours), !is.na(Z1), !is.na(Z2),
         !is.na(Y_P_0), Y_P_0 >= 100000, Y_P_0 <= 150000)
rb_husb_test <- tryCatch(fixest::feols(
  his_hours ~ Z1 + Z2 | couple_id + wave,
  data = rb_husb, cluster = ~ couple_id, notes = FALSE),
  error = function(e) NULL)
if (!is.null(rb_husb_test)) {
  cat("    (a) husband absolute hours on Z1,Z2 (expect null):\n")
  .ct_husb <- fixest::coeftable(rb_husb_test)
  .rows_husb <- intersect(c("Z1", "Z2"), rownames(.ct_husb))
  if (length(.rows_husb)) {
    print(.ct_husb[.rows_husb, , drop = FALSE], digits = 3)
  } else {
    cat("        [note] Z1/Z2 dropped (collinear after FE); nothing to report ",
        "(can occur on synthetic data).\n", sep = "")
  }
}
iv_rb$falsif_husband_hours <- rb_husb_test

# (b) Pre-reform placebo first stage: does the (counterfactual) instrument
#     predict the wife's share BEFORE the reform (W12-W14)? Needs a W12-W24
#     panel carrying person_id, partner_id, female, hours_worked_clean, wave,
#     and youngest-child age. Z1 is built with the ACTUAL per-wave child age so
#     it carries within-couple variation (via child ageing), exactly as the
#     real first stage does; freezing child age would make Z1 constant within
#     couple and collinear with the couple FE.
rb_placebo <- NULL
rb_pre_panel_candidates <- c("hilda_panel_data_W12_W24_slim.rds",
                             "hilda_panel_data_W12_W24.rds")
rb_pre_panel_file <- rb_pre_panel_candidates[file.exists(rb_pre_panel_candidates)][1]
if (is.na(rb_pre_panel_file)) {
  cat("    (b) [skip] no W12-W24 panel found in working directory.\n")
} else if (!exists("iv_baseline")) {
  cat("    (b) [skip] iv_baseline not available.\n")
} else {
  cat(sprintf("    (b) pre-reform placebo using %s\n", rb_pre_panel_file))
  rb_full <- tryCatch(readRDS(rb_pre_panel_file), error = function(e) NULL)
  rb_age_cols <- intersect(c("age_resident_child_1", "age_youngest_own_child"),
                           names(rb_full))
  rb_pre_couples <- tryCatch({
    p1 <- rb_full %>%
      filter(!is.na(partner_id), partner_id != person_id,
             partner_id %in% rb_full$person_id, wave %in% 12:14) %>%
      select(person_id, partner_id, wave, female, hours_worked_clean,
             dplyr::all_of(rb_age_cols))
    p2 <- rb_full %>% filter(wave %in% 12:14) %>%
      select(person_id, wave, female, hours_worked_clean,
             dplyr::all_of(rb_age_cols))
    j <- p1 %>% dplyr::inner_join(
      p2, by = c("partner_id" = "person_id", "wave" = "wave"),
      suffix = c("_f", "_m")) %>%
      filter(female_f != female_m) %>%
      mutate(
        couple_id = ifelse(person_id < partner_id,
                           paste0(person_id, "_", partner_id),
                           paste0(partner_id, "_", person_id)),
        her_hours = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
        his_hours = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
        S_ct = ifelse((his_hours + her_hours) > 0,
                      her_hours / (his_hours + her_hours), NA_real_))
    # youngest-child age per couple-wave: min of available child-age fields
    age_f <- if ("age_resident_child_1_f"   %in% names(j)) j$age_resident_child_1_f   else NA_real_
    age_m <- if ("age_resident_child_1_m"   %in% names(j)) j$age_resident_child_1_m   else NA_real_
    ay_f  <- if ("age_youngest_own_child_f" %in% names(j)) j$age_youngest_own_child_f else NA_real_
    ay_m  <- if ("age_youngest_own_child_m" %in% names(j)) j$age_youngest_own_child_m else NA_real_
    j$a_young_ct <- pmin(dplyr::coalesce(age_f, Inf), dplyr::coalesce(age_m, Inf),
                         dplyr::coalesce(ay_f, Inf), dplyr::coalesce(ay_m, Inf))
    j$a_young_ct <- ifelse(is.infinite(j$a_young_ct), NA_real_, j$a_young_ct)
    j %>% dplyr::distinct(couple_id, wave, .keep_all = TRUE) %>%
      dplyr::left_join(iv_baseline, by = "couple_id")   # frozen Y_P_0, Y_S_0
  }, error = function(e) { cat(sprintf("    (b) build error: %s\n", conditionMessage(e))); NULL })
  
  if (!is.null(rb_pre_couples)) {
    rb_pre_couples$is_grandparent_carer <- FALSE
    # Counterfactual Reform-1 threshold instrument at ACTUAL per-wave child age.
    rb_pre_couples$Z1 <-
      rb_sim_with(rb_pre_couples, rb_R1_params) -
      rb_sim_with(rb_pre_couples, iv_pre_params)
    rb_pre_couples <- rb_pre_couples %>%
      filter(!is.na(S_ct), !is.na(Z1), !is.na(Y_P_0), !is.na(a_young_ct),
             Y_P_0 >= 100000, Y_P_0 <= 150000)
    # Differential pre-trend test: a couple-level reform-withdrawal intensity
    # (Zint, scalar per couple) interacted with a wave trend. The interaction
    # always carries within-couple variation, so it is identified where a level
    # first stage on the near-constant Z1 would not be. A null interaction means
    # high- and low-exposure couples were on parallel S_ct paths before reform.
    rb_pre_couples <- rb_pre_couples %>%
      dplyr::group_by(couple_id) %>%
      dplyr::mutate(Zint = mean(Z1, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      mutate(wave_c = wave - mean(wave, na.rm = TRUE))
    n_cpl <- dplyr::n_distinct(rb_pre_couples$couple_id)
    cat(sprintf("    (b) pre-reform affected band: %d obs, %d couples\n",
                nrow(rb_pre_couples), n_cpl))
    rb_placebo <- tryCatch(fixest::feols(
      S_ct ~ Zint:wave_c | couple_id + wave,
      data = rb_pre_couples, cluster = ~ couple_id, notes = FALSE),
      error = function(e) { cat(sprintf("    (b) feols error: %s\n", conditionMessage(e))); NULL })
  }
}
rb_placebo_term <- "Zint:wave_c"
if (!is.null(rb_placebo) && rb_placebo_term %in% names(coef(rb_placebo))) {
  cat("    (b) pre-reform differential trend (S_ct on intensity x wave; expect null):\n")
  print(fixest::coeftable(rb_placebo)[rb_placebo_term, , drop = FALSE], digits = 3)
} else if (!is.null(rb_placebo)) {
  cat("    (b) [warn] pre-trend term absorbed by fixed effects.\n")
}
iv_rb$placebo_pretrend <- rb_placebo

# (c) Donut-hole: drop couples whose baseline primary income sits right at the
#     $100k cliff (noisiest / most manipulable), confirm first stage survives.
rb_donut_fs <- function(lo, hi) {
  d <- rb_reform1 %>% filter(!(Y_P_0 > lo & Y_P_0 < hi))
  m <- tryCatch(fixest::feols(
    ln_wage ~ educ_exp_c + experience_years + age_sq + married_num
    | person_couple + wave
    | S_ct + S_educ_exp ~ Z1 + Z2 + Z1_educ_exp + Z2_educ_exp,
    data = d, cluster = ~ couple_id, notes = FALSE), error = function(e) NULL)
  c(F_S_ct = iv_pull_F(m, "S_ct"),
    n = if (!is.null(m)) m$nobs else NA_integer_,
    couples = dplyr::n_distinct(d$couple_id))
}
rb_donut_tab <- rbind(
  `+/- $5k`  = rb_donut_fs(95000, 105000),
  `+/- $10k` = rb_donut_fs(90000, 110000)
)
cat("    (c) donut-hole first-stage F on S_ct (expect first stage to survive):\n")
print(rb_donut_tab, digits = 3)
iv_rb$donut <- rb_donut_tab

# -----------------------------------------------------------------------------
# G.15  INCOME-EFFECT MITIGATION
#   The reform removes thousands of dollars of benefit; the level loss could
#   touch wages through non-labour channels. Add the baseline benefit LEVEL
#   (B_pre) as a control so identification is the CHANGE net of the level, and
#   run a placebo on a non-labour outcome if one is available in the panel.
# -----------------------------------------------------------------------------
cat("\n  -- G.15 income-effect mitigation -------------------------------\n")

rb_iv_bctrl <- tryCatch(fixest::feols(as.formula(paste0(
  "ln_wage ~ ", paste(c(RB_CTRL, "B_pre"), collapse = " + "), " | ", RB_FE,
  " | S_ct + S_educ_exp ~ Z1 + Z2 + Z1_educ_exp + Z2_educ_exp")),
  data = rb_reform1, cluster = ~ couple_id, notes = FALSE),
  error = function(e) NULL)
if (!is.null(rb_iv_bctrl)) {
  cat("    IV with baseline-benefit-level control (theta on S_ct):\n")
  .ct_bctrl <- fixest::coeftable(rb_iv_bctrl)
  if ("fit_S_ct" %in% rownames(.ct_bctrl)) {
    print(.ct_bctrl["fit_S_ct", , drop = FALSE], digits = 3)
  } else {
    cat("        [note] fit_S_ct not in the coefficient table ",
        "(can occur on synthetic data).\n", sep = "")
  }
}
iv_rb$income_effect_bctrl <- rb_iv_bctrl

# Non-labour wellbeing placebo (income-effect exclusion check).
#   Reads an external LONG-format satisfaction extract with columns
#   (xwaveid OR person_id), wave, and one or more losat* variables. Save it as
#   hilda_losat_long.rds or .csv in the working directory. If the extract keys
#   on xwaveid, it is crosswalked to the panel's person_id. Wife-side outcomes
#   are joined onto iv_dual_earner and tested on the instruments over the
#   affected band. losat and the off-channel items are placebos (expect null);
#   losatfs is reported separately as an income-channel diagnostic (it may move
#   if the benefit cut is subjectively felt, which does not by itself violate
#   exclusion). losatft / losateo are intentionally excluded as labour-channel.
#
# Source for the satisfaction columns, in priority order:
#   (1) an external long extract (hilda_losat_long.rds/.csv), if present;
#   (2) iv_panel in memory, if the rebuilt panel carries losat* columns;
#   (3) the slim panel file on disk, if it carries losat* columns.
# Paths (2) and (3) require the panel to have been rebuilt with the loader
# that adds the losat* domains, so no separate extract step is needed.
rb_sat_file <- c("hilda_losat_long.rds", "hilda_losat_long.csv",
                 "losat_long.rds", "losat_long.csv")
rb_sat_file <- rb_sat_file[file.exists(rb_sat_file)][1]
rb_sat <- if (!is.na(rb_sat_file)) tryCatch({
  if (grepl("\\.csv$", rb_sat_file)) utils::read.csv(rb_sat_file, stringsAsFactors = FALSE)
  else readRDS(rb_sat_file)
}, error = function(e) NULL) else NULL

# Fallback (2): pull losat* straight from the in-memory extended panel.
if (is.null(rb_sat) && exists("iv_panel")) {
  iv_panel_sat <- grep("^losat", names(iv_panel), value = TRUE)
  if (length(iv_panel_sat) && all(c("person_id", "wave") %in% names(iv_panel))) {
    rb_sat <- iv_panel[, c("person_id", "wave", iv_panel_sat)]
    cat("    G.15 satisfaction source: in-memory iv_panel (extended).\n")
  }
}
# Fallback (3): read the slim panel and take the losat* columns.
if (is.null(rb_sat) && file.exists("hilda_panel_data_W12_W24_slim.rds")) {
  .sp <- tryCatch(readRDS("hilda_panel_data_W12_W24_slim.rds"), error = function(e) NULL)
  if (!is.null(.sp)) {
    sp_sat <- grep("^losat", names(.sp), value = TRUE)
    if (length(sp_sat) && all(c("person_id", "wave") %in% names(.sp))) {
      rb_sat <- .sp[, c("person_id", "wave", sp_sat)]
      cat("    G.15 satisfaction source: slim panel on disk.\n")
    }
  }
}

iv_rb$wellbeing_placebo <- NULL
iv_rb$income_channel_check <- NULL
if (is.null(rb_sat)) {
  cat("    [skip] no losat* source found: rebuild the panel with the loader that\n")
  cat("           adds the satisfaction domains, or place hilda_losat_long.rds in wd.\n")
} else {
  rb_sat <- as.data.frame(rb_sat)
  sat_vars <- grep("^losat", names(rb_sat), value = TRUE)
  # resolve a person_id key (.pid) for the extract
  if ("person_id" %in% names(rb_sat)) {
    rb_sat$.pid <- rb_sat$person_id
  } else if ("xwaveid" %in% names(rb_sat)) {
    xwalk <- NULL
    if (exists("iv_panel") && all(c("person_id", "xwaveid") %in% names(iv_panel))) {
      xwalk <- unique(iv_panel[, c("person_id", "xwaveid")])
    } else if (file.exists("hilda_panel_data_W12_W24_slim.rds")) {
      .sp <- tryCatch(readRDS("hilda_panel_data_W12_W24_slim.rds"), error = function(e) NULL)
      if (!is.null(.sp) && all(c("person_id", "xwaveid") %in% names(.sp)))
        xwalk <- unique(.sp[, c("person_id", "xwaveid")])
    }
    if (!is.null(xwalk)) rb_sat <- dplyr::left_join(rb_sat, xwalk, by = "xwaveid")
    rb_sat$.pid <- if ("person_id" %in% names(rb_sat)) rb_sat$person_id else NA
  } else rb_sat$.pid <- NA
  
  if (!length(sat_vars) || all(is.na(rb_sat$.pid))) {
    cat("    [skip] extract has no losat* variables or no usable id crosswalk.\n")
  } else {
    satL <- rb_sat[, c(".pid", "wave", sat_vars)]
    satL <- satL[!is.na(satL$.pid) & !duplicated(satL[, c(".pid", "wave")]), ]
    de <- iv_dual_earner
    de$her_pid <- ifelse(de$female_f == 1, de$person_id, de$partner_id)
    de$his_pid <- ifelse(de$female_f == 0, de$person_id, de$partner_id)
    her <- satL; names(her) <- c("her_pid", "wave", paste0("her_", sat_vars))
    de  <- dplyr::left_join(de, her, by = c("her_pid", "wave"))
    band <- de %>% filter(hit_by_limit, !hit_by_age13)
    
    rb_sat_run <- function(v) {
      yv <- paste0("her_", v)
      if (!yv %in% names(band) || all(is.na(band[[yv]]))) return(NULL)
      m <- tryCatch(fixest::feols(
        as.formula(paste0(yv, " ~ Z1 + Z2 | couple_id + wave")),
        data = band, cluster = ~ couple_id, notes = FALSE),
        error = function(e) NULL)
      if (is.null(m) || !all(c("Z1", "Z2") %in% names(coef(m)))) return(NULL)
      ct <- fixest::coeftable(m)
      data.frame(outcome = v,
                 Z1 = ct["Z1", "Estimate"], Z1_p = ct["Z1", "Pr(>|t|)"],
                 Z2 = ct["Z2", "Estimate"], Z2_p = ct["Z2", "Pr(>|t|)"],
                 n = m$nobs, row.names = NULL)
    }
    placebo_items <- intersect(c("losat", "losathl", "losatnl", "losatsf"), sat_vars)
    channel_items <- intersect(c("losatfs"), sat_vars)
    rb_sat_placebo <- do.call(rbind, lapply(placebo_items, rb_sat_run))
    rb_sat_channel <- do.call(rbind, lapply(channel_items, rb_sat_run))
    if (!is.null(rb_sat_placebo)) {
      cat("    wellbeing placebos (wife outcome on Z1,Z2; expect null):\n")
      print(rb_sat_placebo, digits = 3)
    } else cat("    [note] placebo outcomes produced no rows ",
               "(all-NA after join or instruments collinear; can occur on synthetic data).\n",
               sep = "")
    if (!is.null(rb_sat_channel)) {
      cat("    income-channel diagnostic losatfs (may move if the cut is felt):\n")
      print(rb_sat_channel, digits = 3)
    }
    iv_rb$wellbeing_placebo <- rb_sat_placebo
    iv_rb$income_channel_check <- rb_sat_channel
  }
}

# -----------------------------------------------------------------------------
# G.16  COMPLIER CHARACTERISATION  (Abadie 2003 kappa)
#
# Abadie's kappa is defined for a BINARY instrument and binary treatment, so we
# profile compliers descriptively on a binarised version: instrument
# Zb = hit_by_limit (in affected band, child <13), treatment Db = 1{S_ct below
# the sample median} (i.e. "specialised toward him"). kappa-weighted covariate
# means describe who the compliers are; this is a descriptive complier profile,
# NOT the structural continuous-treatment LATE.
# -----------------------------------------------------------------------------
cat("\n  -- G.16 complier characterisation (Abadie kappa) ---------------\n")

rb_kappa_df <- iv_dual_earner %>%
  filter(has_children == 1, !is.na(S_ct), !is.na(Y_P_0)) %>%
  mutate(
    Zb = as.integer(hit_by_limit & !hit_by_age13),
    Db = as.integer(S_ct < median(S_ct, na.rm = TRUE))
  )
rb_kappa_covs <- intersect(
  c("her_educ", "his_educ", "her_exp", "his_exp", "a_young_ct",
    "her_age", "his_age", "Y_P_0", "Y_S_0"),
  names(rb_kappa_df))
rb_complier <- tryCatch({
  ps_fit <- glm(reformulate(rb_kappa_covs, response = "Zb"),
                data = rb_kappa_df, family = binomial())
  pZ <- predict(ps_fit, type = "response")
  pZ <- pmin(pmax(pZ, 0.02), 0.98)
  kap <- 1 - rb_kappa_df$Db * (1 - rb_kappa_df$Zb) / (1 - pZ) -
    (1 - rb_kappa_df$Db) * rb_kappa_df$Zb / pZ
  data.frame(
    covariate     = rb_kappa_covs,
    sample_mean   = sapply(rb_kappa_covs, function(v) mean(rb_kappa_df[[v]], na.rm = TRUE)),
    complier_mean = sapply(rb_kappa_covs, function(v)
      sum(kap * rb_kappa_df[[v]], na.rm = TRUE) / sum(kap, na.rm = TRUE)),
    row.names = NULL)
}, error = function(e) NULL)
if (!is.null(rb_complier)) {
  cat("    Complier vs full-sample covariate means (kappa-weighted):\n")
  print(rb_complier, digits = 3)
}
iv_rb$complier_profile <- rb_complier

# -----------------------------------------------------------------------------
# G.17  LEE (2009) BOUNDS on the employment-selection margin
#
# Wages are observed only for the employed; the reform shifts employment, so
# the wage effect is bounded by trimming the more-selected group. Treatment is
# binarised (Db = affected by Reform-1 band). Sharp Lee bounds with no
# covariates; monotonicity direction inferred from the selection gap. This
# documents how severe the wage underpowering is in bound form.
#
# Frame note. This bound MUST start from the pre-selection couple frame
# (iv_couples), NOT iv_dual_earner. iv_dual_earner is filtered on
# !is.na(her_lnwage) (and her_hours > 0), so wife employment there is constant
# at 1 and the selection wedge is mechanically zero -- the bound is undefined by
# construction. iv_couples retains couples where the wife is not employed, so
# selection (her_employed) genuinely varies across the reform band and the trim
# fraction is identified. The band flags (hit_by_limit / hit_by_age13) live only
# on iv_dual_earner, so we reconstruct them here on iv_couples by re-attaching
# baseline incomes (iv_baseline) and re-deriving the same two conditions used at
# G.4; the age-13 condition reuses iv_pre/post wave constants exactly as above.
# -----------------------------------------------------------------------------
cat("\n  -- G.17 Lee (2009) bounds (employment selection) ---------------\n")

rb_lee <- tryCatch({
  # Pre-selection couple-waves with children; selection = wife employment,
  # outcome = wife's ln wage (observed iff employed).
  d <- iv_couples %>%
    dplyr::filter(has_children == 1) %>%
    dplyr::left_join(iv_baseline, by = "couple_id") %>%
    dplyr::filter(!is.na(Y_P_0)) %>%
    iv_mutate(
      hit_by_limit = !is.na(Y_P_0) & Y_P_0 >= 100000 & Y_P_0 <= 150000,
      hit_by_age13 = !is.na(a_young_ct) & a_young_ct >= 13 &
        wave >= IV_REFORM_WAVE_2 & !is_grandparent_carer,
      Db  = as.integer(hit_by_limit & !hit_by_age13),
      # wife employed this couple-wave: prefer the explicit flag, fall back to
      # positive hours, and treat an observed wage as employed (it implies work).
      sel = as.integer(
        dplyr::coalesce(her_employed == 1, FALSE) |
          dplyr::coalesce(her_hours > 0,    FALSE) |
          !is.na(her_lnwage)),
      y   = her_lnwage)
  d <- d %>% dplyr::filter(!is.na(Db))
  s0 <- mean(d$sel[d$Db == 0]); s1 <- mean(d$sel[d$Db == 1])
  n0 <- sum(d$Db == 0); n1 <- sum(d$Db == 1)
  # trim the group with the higher selection rate
  if (is.na(s0) || is.na(s1) || s0 == s1 || min(n0, n1) < 30) {
    NULL
  } else {
    if (s1 > s0) { p <- (s1 - s0) / s1; hi_grp <- 1 } else { p <- (s0 - s1) / s0; hi_grp <- 0 }
    yk  <- d$y[d$Db == hi_grp & d$sel == 1]; yk <- sort(yk[!is.na(yk)])
    yo  <- d$y[d$Db == (1 - hi_grp) & d$sel == 1]; mo <- mean(yo, na.rm = TRUE)
    if (length(yk) < 30 || !is.finite(mo)) {
      NULL
    } else {
      nlo <- floor(length(yk) * p); nhi <- length(yk)
      trimmed_low  <- mean(yk[(nlo + 1):nhi])               # drop lowest p
      trimmed_high <- mean(yk[1:(nhi - nlo)])               # drop highest p
      # bounds on (treated - control) difference in observed-wage means
      if (hi_grp == 1) {
        lower <- trimmed_low  - mo
        upper <- trimmed_high - mo
      } else {
        lower <- mo - trimmed_high
        upper <- mo - trimmed_low
      }
      list(sel_rate_D0 = s0, sel_rate_D1 = s1, trim_fraction = p,
           n_D0 = n0, n_D1 = n1, n_trimmed_grp = length(yk),
           lower = lower, upper = upper)
    }
  }
}, error = function(e) NULL)
if (!is.null(rb_lee)) {
  cat(sprintf("    selection (wife employed): D0=%.3f (n=%d) D1=%.3f (n=%d) ; trim p=%.3f\n",
              rb_lee$sel_rate_D0, rb_lee$n_D0, rb_lee$sel_rate_D1, rb_lee$n_D1,
              rb_lee$trim_fraction))
  cat(sprintf("    Lee bounds on reform wage effect: [%.4f, %.4f]\n",
              rb_lee$lower, rb_lee$upper))
} else {
  cat("    [skip] Lee bounds not computable (no selection gap / too few obs).\n")
}
iv_rb$lee_bounds <- rb_lee

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
saveRDS(iv_rb, "iv_robustness_results.rds")
cat("\n  Saved iv_robustness_results.rds\n")
cat(strrep("=", 78), "\n", sep = "")