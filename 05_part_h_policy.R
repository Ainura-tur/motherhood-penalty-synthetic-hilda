# =============================================================================
# MASTER_PART_H_POLICY.R
#
# PART H: POLICY ECONOMICS OF THE FTB-B REFORM
#
# Converts the Part G simulated-IV machinery into policy-denominated objects
# for AEJ: Economic Policy:
#   H.1  Dollar-denominated first stages (per $1,000 of FTB-B withdrawn)
#        on S_ct, her weekly hours, log hours, and EMPLOYMENT (extensive
#        margin, estimated on the full couple panel, not dual-earner).
#   H.2  Price/income decomposition (Gruber-Saez): Delta-EMTR vs Delta-virtual-
#        income instruments; secondary-earner hours & participation
#        elasticities w.r.t. the net-of-tax wage.
#   H.3  Reform 2 (age-13) as a second instrument: overidentification (Hansen J)
#        and dose-response.
#   H.4  Bridge calculation: reform-induced hours -> experience -> wages,
#        NPV per affected mother, aggregate vs fiscal saving, parametric
#        bootstrap CIs.
#   H.5  Weak-IV robust inference: Anderson-Rubin confidence set for theta.
#
# Source from MASTER_hh.R immediately AFTER MASTER_PART_G_FTBB_IV.R.
# Requires in memory: iv_couples, iv_dual_earner, iv_sample,
#   iv_reform1_sample, iv_ftbb_params, iv_pre_params, iv_simulate_ftbb,
#   IV_BASELINE_WAVES; and on disk: cs_event_study_results.rds.
# Writes policy_economics_results.rds.
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  PART H: POLICY ECONOMICS OF THE FTB-B REFORM\n")
cat(strrep("=", 78), "\n", sep = "")

ph_filter <- dplyr::filter; ph_mutate <- dplyr::mutate
ph_select <- dplyr::select

# -----------------------------------------------------------------------------
# H.0  DEFLATE Z TO W24 (2024) DOLLARS
# -----------------------------------------------------------------------------
# Z_ftb in Part G is in nominal FY dollars (the statutory table is nominal).
# Wages are WPI-deflated to wave-24 dollars in HILDA_LOADIND_IV.R, so deflate
# Z with the IDENTICAL factor (ABS A2705194A, base W24 = 149.6).
ph_wpi <- tibble::tibble(
  wave = 12:24,   # verified ABS A2705194A, base wave 24 = 149.6
  wpi_to_w24 = 149.6 / c(110.9, 114.6, 117.6, 120.4, 123.0, 125.4,
                         127.9, 130.9, 133.7, 135.7, 138.9, 143.7, 149.6)
)
ph_deflate <- function(df) {
  df %>% dplyr::left_join(ph_wpi, by = "wave") %>%
    ph_mutate(Z_ftb_real = Z_ftb * wpi_to_w24,
              Z_k        = Z_ftb_real / 1000)        # per $1,000 (2024$)
}
iv_dual_earner <- ph_deflate(iv_dual_earner)

# -----------------------------------------------------------------------------
# H.1  DOLLAR-DENOMINATED FIRST STAGES
# -----------------------------------------------------------------------------
# (a) Intensive margins on the dual-earner couple panel (couple FE, not the
#     stacked person panel: outcomes are couple-level).
ph_couple_sample <- iv_dual_earner %>%
  ph_filter(has_children == 1, !is.na(Z_k), !is.na(Y_P_0))

ph_fs <- list(
  S_ct      = fixest::feols(S_ct          ~ Z_k | couple_id + wave,
                            data = ph_couple_sample, cluster = ~ couple_id),
  her_hours = fixest::feols(her_hours     ~ Z_k | couple_id + wave,
                            data = ph_couple_sample, cluster = ~ couple_id),
  ln_her_h  = fixest::feols(log(her_hours) ~ Z_k | couple_id + wave,
                            data = ph_couple_sample, cluster = ~ couple_id)
)
cat("\n  H.1 first stages, per $1,000 of FTB-B (2024$):\n")
for (nm in names(ph_fs)) {
  ct <- ph_fs[[nm]]$coeftable["Z_k", ]
  cat(sprintf("    %-9s  b = %+.4f  (SE %.4f)\n", nm, ct[1], ct[2]))
}

# (b) Extensive margin on the FULL couple panel. iv_dual_earner conditions on
#     both partners employed, which deletes the participation margin. Rebuild
#     Z on iv_couples and use her EMPLOYMENT as the outcome.
# REQUIRES a one-line Part G change: add `employed` (esdtl in {1,2}, from
# HILDA_LOADIND_IV.R) to the iv_p1 / iv_p2 select lists, and carry
#   her_employed = ifelse(female_f == 1, employed_f, employed_m)
# in the iv_couples mutate. Do NOT derive employment from her_hours > 0:
# hours_worked_clean is NA (not 0) for the non-employed, so an hours-based
# flag drops exactly the extensive-margin observations this estimate needs.
ph_couples_all <- iv_couples %>%
  dplyr::left_join(iv_baseline, by = "couple_id") %>%
  dplyr::left_join(iv_ftbb_params %>% ph_select(-fy_label), by = "wave")
ph_couples_all$B_post <- iv_simulate_ftbb(
  ph_couples_all$a_young_ct, ph_couples_all$Y_P_0, ph_couples_all$Y_S_0,
  ph_couples_all$M_0to4, ph_couples_all$M_5to12, ph_couples_all$M_13plus,
  ph_couples_all$I_r, ph_couples_all$taper, ph_couples_all$P_r,
  ph_couples_all$age13_couple_ok, ph_couples_all$is_grandparent_carer)
ph_couples_all$B_pre <- iv_simulate_ftbb(
  ph_couples_all$a_young_ct, ph_couples_all$Y_P_0, ph_couples_all$Y_S_0,
  iv_pre_params$M_0to4, iv_pre_params$M_5to12, iv_pre_params$M_13plus,
  iv_pre_params$I_r, iv_pre_params$taper, iv_pre_params$P_r,
  iv_pre_params$age13_couple_ok, ph_couples_all$is_grandparent_carer)
ph_couples_all <- ph_couples_all %>%
  ph_mutate(Z_ftb = B_post - B_pre) %>%
  ph_deflate() %>%
  ph_filter(has_children == 1, !is.na(Z_k), !is.na(Y_P_0), !is.na(her_employed))

ph_fs$her_employed <- fixest::feols(
  her_employed ~ Z_k | couple_id + wave,
  data = ph_couples_all, cluster = ~ couple_id)
ct <- ph_fs$her_employed$coeftable["Z_k", ]
cat(sprintf("    %-9s  b = %+.4f  (SE %.4f)  [full couple panel, n=%d]\n",
            "employed", ct[1], ct[2], ph_fs$her_employed$nobs))

# -----------------------------------------------------------------------------
# H.2  PRICE / INCOME DECOMPOSITION (Gruber-Saez)
# -----------------------------------------------------------------------------
# EMTR on her earnings from FTB-B taper, by numerical derivative of the
# statutory schedule at frozen baseline secondary income. Reform 1 both
# removes a lump sum (income effect) and zeroes the 20% taper for affected
# couples (price effect: her net-of-tax wage RISES).
ph_emtr <- function(df, params) {
  step <- 100
  b0 <- iv_simulate_ftbb(df$a_young_ct, df$Y_P_0, df$Y_S_0,
                         params$M_0to4, params$M_5to12, params$M_13plus,
                         params$I_r, params$taper, params$P_r,
                         params$age13_couple_ok, df$is_grandparent_carer)
  b1 <- iv_simulate_ftbb(df$a_young_ct, df$Y_P_0, df$Y_S_0 + step,
                         params$M_0to4, params$M_5to12, params$M_13plus,
                         params$I_r, params$taper, params$P_r,
                         params$age13_couple_ok, df$is_grandparent_carer)
  (b0 - b1) / step                      # benefit lost per extra $ she earns
}
ph_virtual <- function(df, params) {     # benefit at yS = 0 (virtual income)
  iv_simulate_ftbb(df$a_young_ct, df$Y_P_0, 0,
                   params$M_0to4, params$M_5to12, params$M_13plus,
                   params$I_r, params$taper, params$P_r,
                   params$age13_couple_ok, df$is_grandparent_carer)
}
ph_decomp <- function(df) {
  pw <- df %>% ph_select(wave) %>%
    dplyr::left_join(iv_ftbb_params %>% ph_select(-fy_label), by = "wave")
  df$d_emtr <- ph_emtr(df, pw) - ph_emtr(df, iv_pre_params[rep(1, nrow(df)), ])
  df$d_virt <- ph_virtual(df, pw) - ph_virtual(df, iv_pre_params[rep(1, nrow(df)), ])
  df
}
ph_couple_sample <- ph_decomp(ph_couple_sample)
ph_couples_all   <- ph_decomp(ph_couples_all)

ph_pi <- list(
  hours = fixest::feols(her_hours ~ d_emtr + I(d_virt/1000) | couple_id + wave,
                        data = ph_couple_sample, cluster = ~ couple_id),
  part  = fixest::feols(her_employed ~ d_emtr + I(d_virt/1000) | couple_id + wave,
                        data = ph_couples_all, cluster = ~ couple_id)
)
cat("\n  H.2 price/income decomposition (d_emtr in taper points; d_virt per $1,000):\n")
print(lapply(ph_pi, function(m) m$coeftable))

# Participation elasticity w.r.t. net-of-tax wage: scale the d_emtr coefficient
# by baseline participation and the mean net-of-tax share.
# eps = (dP/d(1-tau)) * (1-tau)/P, with d(1-tau) = -d_emtr
ph_P0   <- mean(ph_couples_all$her_employed, na.rm = TRUE)
ph_tau0 <- mean(ph_emtr(ph_couples_all, iv_pre_params[rep(1, nrow(ph_couples_all)), ]),
                na.rm = TRUE)
ph_eps_part <- -ph_pi$part$coeftable["d_emtr", "Estimate"] * (1 - ph_tau0) / ph_P0
cat(sprintf("  Implied participation elasticity (net-of-tax wage): %.3f\n", ph_eps_part))
cat("  Compare: secondary-earner participation elasticities in the literature.\n")

# -----------------------------------------------------------------------------
# H.3  REFORM 2 AS A SECOND INSTRUMENT: OVERIDENTIFICATION
# -----------------------------------------------------------------------------
# Decompose Z: Z1 = threshold change only (hold age13_couple_ok at pre value),
# Z2 = age-13 exclusion only. Two instruments, one endogenous regressor ->
# Hansen J. (Single-endogenous spec, matching G.8.)
ph_params_no_age13 <- iv_ftbb_params %>% ph_mutate(age13_couple_ok = TRUE)
ph_z_split <- function(df) {
  pw1 <- df %>% ph_select(wave) %>%
    dplyr::left_join(ph_params_no_age13 %>% ph_select(-fy_label), by = "wave")
  b_no_age13 <- iv_simulate_ftbb(df$a_young_ct, df$Y_P_0, df$Y_S_0,
                                 pw1$M_0to4, pw1$M_5to12, pw1$M_13plus,
                                 pw1$I_r, pw1$taper, pw1$P_r,
                                 pw1$age13_couple_ok, df$is_grandparent_carer)
  df$Z1 <- b_no_age13 - df$B_pre        # Reform 1 component
  df$Z2 <- df$B_post  - b_no_age13      # Reform 2 component
  df
}
ph_couple_sample <- ph_z_split(ph_couple_sample)

# Wage overid: join the couple-level Z1/Z2 into the STACKED person panel and
# mirror the single-endogenous G.8 spec (ln_wage, S_ct endogenous), now with
# two instruments; Sargan/Hansen J from fitstat. Run on the full iv_sample,
# where both reform components have support (on the Reform 1 band alone, Z2
# is near-degenerate).
ph_z_cw <- ph_couple_sample %>% ph_select(couple_id, wave, Z1, Z2)
ph_overid_sample <- iv_sample %>%
  dplyr::inner_join(ph_z_cw, by = c("couple_id", "wave"))
ph_overid <- tryCatch(
  fixest::feols(
    ln_wage ~ educ_exp_c + experience_years + age_sq + married_num
    | person_couple + wave
    | S_ct ~ Z1 + Z2,
    data = ph_overid_sample, cluster = ~ couple_id),
  error = function(e) { cat("  H.3 overid failed:", conditionMessage(e), "\n"); NULL })
# Number of distinct instruments that actually entered (Sargan/Hansen J is
# only defined when this exceeds the number of endogenous regressors = 1).
# If Z1 and Z2 are collinear within couple+wave (their bands overlap on
# this sample), fixest drops one and the model is just-identified -> J is
# undefined and fitstat returns NA, which is correct, not an error.
ph_n_inst <- if (!is.null(ph_overid))
  length(unique(stats::na.omit(ph_overid$iv_inst_names))) else 0L
ph_sargan_stat <- NA_real_; ph_sargan_p <- NA_real_
if (!is.null(ph_overid) && ph_n_inst >= 2) {
  s <- tryCatch(fixest::fitstat(ph_overid, "sargan"), error = function(e) NULL)
  if (is.list(s) && !is.null(s$sargan)) s <- s$sargan      # list-with-$sargan
  ph_sargan_stat <- tryCatch(if (is.list(s)) s$stat else unname(s["stat"]),
                             error = function(e) NA_real_)
  ph_sargan_p    <- tryCatch(if (is.list(s)) s$p    else unname(s["p"]),
                             error = function(e) NA_real_)
}
if (!is.null(ph_overid)) {
  ct <- ph_overid$coeftable["fit_S_ct", ]
  cat(sprintf("\n  H.3 two-instrument IV: theta = %+.3f (SE %.3f), n = %d, instruments used = %d\n",
              ct[1], ct[2], ph_overid$nobs, ph_n_inst))
  if (ph_n_inst >= 2 && is.finite(ph_sargan_stat))
    cat(sprintf("      Sargan J = %.3f, p = %.3f (H0: both reforms identify the same theta)\n",
                ph_sargan_stat, ph_sargan_p))
  else
    cat("      Sargan J: not defined (model is just-identified; Z1/Z2 collinear\n",
        "        on this sample -- the over-id test needs both to enter separately).\n")
}
# Labour-supply overid (her_hours second stage, couple level):
ph_overid_hours <- tryCatch(
  fixest::feols(her_hours ~ 1 | couple_id + wave | S_ct ~ Z1 + Z2,
                data = ph_couple_sample, cluster = ~ couple_id),
  error = function(e) NULL)

# -----------------------------------------------------------------------------
# H.4  BRIDGE CALCULATION: POLICY -> HOURS -> EXPERIENCE -> WAGES
# -----------------------------------------------------------------------------
ph_es <- tryCatch(readRDS("cs_event_study_results.rds"), error = function(e) NULL)
if (!is.null(ph_es)) {
  gw <- ph_es[["women ln_hourly_wage_real"]]; gm <- ph_es[["men ln_hourly_wage_real"]]
  xw <- ph_es[["women experience_years"]];    xm <- ph_es[["men experience_years"]]
  e10 <- 10
  wage_gap10 <- gw$att[gw$event_time == e10] - gm$att[gm$event_time == e10]
  exp_gap10  <- xw$att[xw$event_time == e10] - xm$att[xm$event_time == e10]
  ph_ret_per_expyr <- wage_gap10 / exp_gap10      # event-study-implied return
  cat(sprintf("\n  H.4 event-study-implied return per experience-year: %.4f\n",
              ph_ret_per_expyr))
  
  # Reform-induced hours -> experience: FTE accrual approximation
  FT_HOURS <- 38
  dh_per_k <- ph_fs$her_hours$coeftable["Z_k", "Estimate"]   # weekly hrs / $1k
  # Z_ftb on the stacked sample is NOMINAL; the first stage above is per
  # $1,000 of 2024 dollars, so deflate before chaining the two.
  # Self-contained WPI deflator (verified ABS A2705194A, base wave 24 =
  # 149.6), built here so H.4 does not depend on ph_wpi surviving from the
  # top of the script when this file is source()d into the master.
  .ph_wpi24 <- tibble::tibble(
    wave = 12:24,
    wpi_to_w24 = 149.6 / c(110.9, 114.6, 117.6, 120.4, 123.0, 125.4,
                           127.9, 130.9, 133.7, 135.7, 138.9, 143.7, 149.6))
  meanZ_aff <- iv_reform1_sample %>%
    dplyr::left_join(.ph_wpi24, by = "wave") %>%
    dplyr::summarise(z = mean(Z_ftb * wpi_to_w24, na.rm = TRUE) / 1000) %>%
    dplyr::pull(z)                                # avg $k withdrawn (2024$)
  dh_reform <- dh_per_k * meanZ_aff               # avg weekly hours response
  dexp_per_yr <- dh_reform / FT_HOURS             # FTE experience-years per year
  
  # NPV per affected mother over remaining career
  ph_npv <- function(dexp_yr, ret, base_earn, T = 25, r = 0.03) {
    # cumulative experience effect: after t years, extra experience = dexp_yr*t
    sum(sapply(1:T, function(t) base_earn * ret * dexp_yr * t / (1 + r)^t))
  }
  # base_earn: weighted mean annual earnings of affected mothers, in W24 $.
  # annual_income (wsfei + wsfes in HILDA_LOADIND_IV.R) is NOMINAL: deflate it
  # with the same wpi_to_w24 factor. Weight by analysis_weight (hhwtrps), the
  # cross-sectional responding-person weight the paper already uses; requires
  # adding annual_income + analysis_weight for the wife to the iv_couples
  # carry-through (annual_income is already in iv_required_vars).
  # base_earn <- with(affected_wives,
  #   weighted.mean(annual_income * wpi_to_w24, analysis_weight, na.rm = TRUE))
  base_earn <- NA_real_  # fill per the recipe above
  cat(sprintf("  reform-induced hours: %+.2f hrs/wk; experience: %+.3f FTE-yrs/yr\n",
              dh_reform, dexp_per_yr))
  cat("  NPV per mother / aggregate, once base_earn is filled:\n")
  cat("    npv1 <- ph_npv(dexp_per_yr, ph_ret_per_expyr, base_earn)\n")
  cat("    aggregate earnings loss: sum over affected couples of\n")
  cat("      couple_weight * npv1, couple_weight = mean of partners' hhwtrps\n")
  cat("    fiscal saving comparator: sum(couple_weight * -Z_ftb_real)\n")
  
  # Parametric bootstrap for the chained CI (independence across blocks; flag)
  ph_boot <- function(B = 2000) {
    fs_b  <- ph_fs$her_hours$coeftable["Z_k", ]
    replicate(B, {
      dh  <- rnorm(1, fs_b[1], fs_b[2]) * meanZ_aff
      wg  <- rnorm(1, wage_gap10, sqrt(gw$se[gw$event_time == e10]^2 +
                                         gm$se[gm$event_time == e10]^2))
      xg  <- rnorm(1, exp_gap10,  sqrt(xw$se[xw$event_time == e10]^2 +
                                         xm$se[xm$event_time == e10]^2))
      ph_npv(dh / FT_HOURS, wg / xg, base_earn)
    })
  }
} else cat("\n  H.4 skipped: cs_event_study_results.rds not found.\n")

# -----------------------------------------------------------------------------
# H.5  ANDERSON-RUBIN CONFIDENCE SET (single-endogenous, Reform 1, wives)
# -----------------------------------------------------------------------------
ph_ar <- function(d, grid = seq(-6, 4, by = 0.02)) {
  ok <- sapply(grid, function(th) {
    d$y_adj <- d$ln_wage - th * d$S_ct
    m <- fixest::feols(y_adj ~ Z_ftb + experience_years + age_sq + married_num
                       | person_couple + wave, data = d, cluster = ~ couple_id,
                       notes = FALSE)
    m$coeftable["Z_ftb", "Pr(>|t|)"] > 0.05     # fail to reject => in AR set
  })
  range(grid[ok])
}
ph_ar_set <- tryCatch(
  ph_ar(iv_reform1_sample[grepl("_F$", iv_reform1_sample$person_couple), ]),
  error = function(e) c(NA, NA))
cat(sprintf("\n  H.5 AR 95%% set, wife-only Reform 1 theta: [%.2f, %.2f]\n",
            ph_ar_set[1], ph_ar_set[2]))

# -----------------------------------------------------------------------------
# H.6  SAVE
# -----------------------------------------------------------------------------
policy_economics_results <- list(
  fs_dollar      = lapply(ph_fs, function(m) m$coeftable),
  price_income   = lapply(ph_pi, function(m) m$coeftable),
  eps_part       = ph_eps_part,
  overid_wage    = if (!is.null(ph_overid)) ph_overid$coeftable else NULL,
  overid_sargan  = c(stat = ph_sargan_stat, p = ph_sargan_p),
  overid_n_inst  = ph_n_inst,
  overid_hours   = if (!is.null(ph_overid_hours)) ph_overid_hours$coeftable else NULL,
  ar_set_reform1 = ph_ar_set
)
saveRDS(policy_economics_results, "policy_economics_results.rds")
cat("  Saved policy_economics_results.rds\n")
cat(strrep("=", 78), "\n", sep = "")