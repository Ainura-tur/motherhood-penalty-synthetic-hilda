# =============================================================================
# MASTER_PART_O_revisions.R
#
# Implements the two outstanding revision-plan items (revision_plan.md):
#
#   O-1  Align estimators in the gender contrast (Online Appendix Table A3).
#        Re-estimate so men's and women's pre/post-parenthood beta2 both use
#        the SAME estimator on the transitioner subsample, producing a clean
#        like-for-like contrast and removing the "within-person OLS for men vs
#        cross-sectional FE for women" asymmetry a referee could read as
#        confounded by estimator choice. Writes a table-note caveat string too,
#        so either resolution path (re-run OR note) is available from one run.
#
#   O-2  Broaden P3 beyond dual-earner couples (Section V.F / tab:mincer_gap).
#        Apply a selection correction to the within-couple betaX (return to
#        experience) gender-gap comparison itself, not merely the pre-birth
#        check already in the paper. Two complementary corrections:
#          (a) Heckman-type IMR control built from a couple-level dual-earner
#              selection probit, added to the linear within-couple betaX gap;
#          (b) Lee (2009) trimming bounds on the within-couple betaX gap under
#              the differential dual-earner-exit wedge, mirroring the
#              event-study Lee bounds in CS_event_study.R.
#
# RUN MODE
#   Preferred: source() AFTER MASTER_hh.R in the SAME session, so the in-memory
#   objects (panel, couples, stacked, stacked_w, men_trans_data,
#   women_trans_data, run_dml, .dml_gender_gap, rec/stars_fn) are reused.
#   Standalone: if those objects are absent, the script rebuilds the minimal
#   ones it needs from hilda_panel_data_W12_W24_slim.rds, replicating the
#   master's couple-construction and stacking logic.
#
# OUTPUTS (written to getwd())
#   se_ledger_A3_aligned.csv        O-1 like-for-like men-vs-women beta2 contrast
#   tableA3_estimator_note.txt      O-1 ready-to-paste table-note caveat
#   mincer_gap_betaX_selcorr.csv    O-2 selection-corrected betaX gender gap
#   betaX_lee_bounds.csv            O-2 Lee bounds on the within-couple betaX gap
#
# Neither item changes any headline number; both are appendix / robustness.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

# ── Namespace protection (avoid MASS/stats masking dplyr verbs) ───────────────
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

# stars helper (reuse master's if present)
if (!exists("stars_fn") || !is.function(get("stars_fn"))) {
  stars_fn <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                           ifelse(p < 0.05, "*", ifelse(p < 0.10, "\u2020", "")))))
  }
}
# results-registry recorder: no-op shim if the master's rec() is absent
if (!exists("rec") || !is.function(get("rec"))) rec <- function(...) invisible(NULL)

# Data-object existence guard. A bare exists() is unsafe here: several packages
# (e.g. lattice) export a function named `panel`, so exists("panel") can be TRUE
# while `panel` is that function rather than our data frame. Require an actual
# data frame in the calling scope before treating the object as "already built".
.has_df <- function(nm) {
  exists(nm, inherits = TRUE) &&
    is.data.frame(tryCatch(get(nm, inherits = TRUE),
                           error = function(e) NULL))
}

cat("\n"); cat(strrep("#", 80), "\n")
cat("#  PART O: REVISION ITEMS O-1 (Table A3 estimator alignment) and\n")
cat("#          O-2 (selection-corrected within-couple betaX gender gap)\n")
cat(strrep("#", 80), "\n")

# =============================================================================
# 0.  SELF-HEALING DATA BUILD (only runs the pieces that are missing)
# =============================================================================
# Hardened object guards. `.has_df()` only checks data-frame-ness, so an
# in-memory object that EXISTS but is the wrong shape (e.g. `panel` after a
# later module subset it to a handful of columns, or a truncated frame) would
# pass and then error downstream. `.frame_ok()` additionally requires the
# columns and the minimum row count each frame must carry; `.need_build()`
# records, per object, whether it was reused or rebuilt and why, stops instead
# of rebuilding when options(part_o_strict = TRUE), and otherwise warns loudly
# before overwriting an inadequate object with a fresh rebuild.
.REQ <- list(
  panel = c("person_id","wave","female","ever_parent","in_wage_sample",
            "ln_hourly_wage_real","educ_years","experience_years",
            "hours_worked_clean","employed","fulltime","age_sq","married"),
  wage  = c("person_id","wave","female","ln_hourly_wage_real","educ_years",
            "experience_years","age_sq","married"),
  men_trans_data   = c("person_id","wave","period","ln_hourly_wage_real",
                       "educ_years","experience_years","age_sq","married"),
  women_trans_data = c("person_id","wave","period","ln_hourly_wage_real",
                       "educ_years","experience_years","age_sq","married"),
  couples = c("couple_id","wave","has_children","his_in_wage","her_in_wage",
              "his_wage","her_wage","his_educ","her_educ","his_exp","her_exp",
              "his_married","her_married","his_age_sq","her_age_sq",
              "his_hours","her_hours","his_fulltime","her_fulltime"))
.MINROWS <- c(panel = 1000L, wage = 1000L, men_trans_data = 50L,
              women_trans_data = 50L, couples = 500L)

.po_strict <- isTRUE(getOption("part_o_strict", FALSE))
.po_prov   <- list()   # per-object provenance, summarised at end of section 0

.frame_ok <- function(nm) {
  if (!.has_df(nm)) return(list(ok = FALSE, why = "absent or not a data frame"))
  d    <- get(nm, inherits = TRUE)
  miss <- setdiff(.REQ[[nm]], names(d))
  if (length(miss))
    return(list(ok = FALSE,
                why = paste0("missing column(s): ", paste(miss, collapse = ", "))))
  mr <- .MINROWS[[nm]]; if (is.null(mr) || is.na(mr)) mr <- 1L
  if (nrow(d) < mr)
    return(list(ok = FALSE,
                why = sprintf("only %s rows (< %s expected)",
                              format(nrow(d), big.mark = ","),
                              format(mr, big.mark = ","))))
  list(ok = TRUE, why = sprintf("reused from memory (%s rows)",
                                format(nrow(d), big.mark = ",")))
}

# TRUE when `nm` must be (re)built. Records provenance; honours strict mode.
.need_build <- function(nm) {
  chk <- .frame_ok(nm)
  .po_prov[[nm]] <<- c(source = if (chk$ok) "reused" else "rebuilt",
                       detail = chk$why)
  if (!chk$ok && .po_strict)
    stop(sprintf(paste0("PART O strict mode: required object '%s' is unusable ",
                        "(%s).\n  Source PART O immediately after PART D, before ",
                        "any module subsets `panel`."), nm, chk$why), call. = FALSE)
  if (!chk$ok && exists(nm, inherits = TRUE))
    warning(sprintf("PART O: in-memory '%s' is unusable (%s); rebuilding it.",
                    nm, chk$why), call. = FALSE)
  !chk$ok
}

# 0a. panel ------------------------------------------------------------------
if (.need_build("panel")) {
  slim <- "hilda_panel_data_W12_W24_slim.rds"
  if (!file.exists(slim))
    stop("panel not in memory and ", slim, " not found. ",
         "Run MASTER_hh.R first, or place the slim RDS in the working directory.")
  cat("\n  [standalone] loading panel from ", slim, "\n", sep = "")
  panel <- readRDS(slim)
  if (!"married_num" %in% names(panel) && "married" %in% names(panel))
    panel <- panel %>% dplyr::mutate(married_num = as.numeric(married))
  # partner_id_clean (same rule as the master)
  if (!"partner_id_clean" %in% names(panel)) {
    if (!"partner_id" %in% names(panel)) panel$partner_id <- NA_character_
    panel <- panel %>%
      dplyr::mutate(partner_id_clean = dplyr::case_when(
        is.na(partner_id)        ~ NA_character_,
        partner_id == ""         ~ NA_character_,
        nchar(partner_id) < 5    ~ NA_character_,
        grepl("^-", partner_id)  ~ NA_character_,
        TRUE                     ~ partner_id))
  }
}

# Ensure couple-link key and married_num exist whether panel was reused or built.
if (!"partner_id_clean" %in% names(panel)) {
  if (!"partner_id" %in% names(panel)) panel$partner_id <- NA_character_
  panel <- panel %>% dplyr::mutate(partner_id_clean = dplyr::case_when(
    is.na(partner_id)       ~ NA_character_,
    partner_id == ""        ~ NA_character_,
    nchar(partner_id) < 5   ~ NA_character_,
    grepl("^-", partner_id) ~ NA_character_,
    TRUE                    ~ partner_id))
}
if (!"married_num" %in% names(panel) && "married" %in% names(panel))
  panel <- panel %>% dplyr::mutate(married_num = as.numeric(married))

# 0b. wage sample ------------------------------------------------------------
if (.need_build("wage")) wage <- panel %>% dplyr::filter(in_wage_sample == 1)

# 0c. transitioner frames (needed for O-1) -----------------------------------
.build_trans <- function(sex) {
  fem <- if (sex == "men") 0L else 1L
  tr <- panel %>%
    dplyr::filter(female == fem) %>%
    dplyr::group_by(person_id) %>%
    dplyr::arrange(wave) %>%
    dplyr::summarise(
      transition_wave = if (any(diff(ever_parent) > 0))
        wave[which(diff(ever_parent) > 0)[1] + 1] else NA_integer_,
      .groups = "drop") %>%
    dplyr::mutate(transitioner = !is.na(transition_wave))
  ids  <- tr %>% dplyr::filter(transitioner) %>% dplyr::pull(person_id)
  lab_pre  <- if (sex == "men") "Pre-fatherhood"  else "Pre-motherhood"
  lab_post <- if (sex == "men") "Post-fatherhood" else "Post-motherhood"
  wage %>%
    dplyr::filter(female == fem, person_id %in% ids) %>%
    dplyr::left_join(tr %>% dplyr::select(person_id, transition_wave), by = "person_id") %>%
    dplyr::mutate(period = dplyr::case_when(wave <  transition_wave ~ lab_pre,
                                            wave >= transition_wave ~ lab_post),
                  years_from_transition = wave - transition_wave)
}
if (.need_build("men_trans_data"))   men_trans_data   <- .build_trans("men")
if (.need_build("women_trans_data")) women_trans_data <- .build_trans("women")

# 0d. couple panel + dual-earner stack (needed for O-2) ----------------------
if (.need_build("couples")) {
  cat("\n  [standalone] rebuilding couple panel\n")
  p1 <- panel %>%
    dplyr::filter(!is.na(partner_id_clean)) %>%
    dplyr::select(person_id, partner_id = partner_id_clean, wave, female, ever_parent,
                  ln_hourly_wage_real, educ_years, experience_years,
                  hours_worked_clean, employed, fulltime, in_wage_sample,
                  married, age_sq)
  p2 <- panel %>%
    dplyr::select(person_id, wave, female, ever_parent, ln_hourly_wage_real,
                  educ_years, experience_years, hours_worked_clean, employed,
                  fulltime, in_wage_sample, married, age_sq)
  couples <- p1 %>%
    dplyr::inner_join(p2, by = c("partner_id" = "person_id", "wave" = "wave"),
                      suffix = c("_f", "_m")) %>%
    dplyr::filter(female_f != female_m) %>%
    dplyr::mutate(
      her_wage = ifelse(female_f == 1, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
      his_wage = ifelse(female_f == 0, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
      her_educ = ifelse(female_f == 1, educ_years_f, educ_years_m),
      his_educ = ifelse(female_f == 0, educ_years_f, educ_years_m),
      her_exp  = ifelse(female_f == 1, experience_years_f, experience_years_m),
      his_exp  = ifelse(female_f == 0, experience_years_f, experience_years_m),
      her_hours = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
      his_hours = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
      her_fulltime = ifelse(female_f == 1, fulltime_f, fulltime_m),
      his_fulltime = ifelse(female_f == 0, fulltime_f, fulltime_m),
      her_age_sq = ifelse(female_f == 1, age_sq_f, age_sq_m),
      his_age_sq = ifelse(female_f == 0, age_sq_f, age_sq_m),
      her_married = ifelse(female_f == 1, married_f, married_m),
      his_married = ifelse(female_f == 0, married_f, married_m),
      her_in_wage = ifelse(female_f == 1, in_wage_sample_f, in_wage_sample_m),
      his_in_wage = ifelse(female_f == 0, in_wage_sample_f, in_wage_sample_m),
      her_emp = ifelse(female_f == 1, employed_f, employed_m),
      his_emp = ifelse(female_f == 0, employed_f, employed_m),
      has_children = as.integer(ever_parent_f == 1 | ever_parent_m == 1),
      couple_id = ifelse(person_id < partner_id,
                         paste0(person_id, "_", partner_id),
                         paste0(partner_id, "_", person_id)))
}

# ---- data-object provenance summary (reused vs rebuilt, per object) ---------
.po_prov_df <- do.call(rbind, lapply(names(.po_prov), function(nm)
  data.frame(object = nm,
             source = unname(.po_prov[[nm]]["source"]),
             detail = unname(.po_prov[[nm]]["detail"]),
             stringsAsFactors = FALSE)))
cat("\n  [PART O] data-object provenance:\n\n")
print(.po_prov_df, row.names = FALSE, right = FALSE)
if (any(.po_prov_df$source == "rebuilt")) {
  cat("\n  NOTE: one or more frames were rebuilt from the slim panel rather than\n",
      "  reused from PART D. The numbers stay valid but may differ in the last\n",
      "  digits from Tables 10/12 if the rebuild and the PART D stacking diverge.\n",
      "  For an exact match, source PART O immediately after PART D, before any\n",
      "  later module subsets `panel`. To make a rebuild a hard error instead,\n",
      "  set options(part_o_strict = TRUE) before sourcing.\n", sep = "")
} else {
  cat("\n  All five frames reused from PART D; PART O matches Tables 10/12.\n")
}
rec("part_o_provenance", .po_prov_df)

# =============================================================================
# O-1.  ESTIMATOR-ALIGNED GENDER CONTRAST (Online Appendix Table A3)
# =============================================================================
# Issue (revision_plan.md O-1): Table A3 contrasts men estimated one way against
# women estimated another, so the gender contrast is confounded with estimator
# choice. Resolution: estimate BOTH sides with the SAME estimator on the SAME
# transitioner subsample, for two estimators side by side, so the contrast is
# like-for-like under either reading.
#
#   Estimator P  (pooled / within-period OLS, wave FE):
#       ln w ~ educ + (educ_c*exp_c) + exp + age_sq + married | wave
#     identifies beta2 off cross-person variation inside the pre (or post) cell.
#
#   Estimator W  (within-person FE, person + wave FE):
#       ln w ~ (educ_c*exp_c) + exp + age_sq + married | person_id + wave
#     identifies beta2 off within-person variation inside the pre (or post) cell.
#
# Both are run on men and on women, pre and post, on the transitioner wage
# sample. The aligned contrast a referee should read is column-by-column:
# men vs women under the SAME estimator.
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n")
cat("  O-1. Estimator-aligned gender contrast (Table A3)\n")
cat(strrep("=", 78), "\n")

.b2_aligned <- function(data, estimator) {
  if (!"married_num" %in% names(data))
    data <- data %>% dplyr::mutate(married_num = as.numeric(married))
  d <- data %>%
    dplyr::filter(!is.na(ln_hourly_wage_real), !is.na(educ_years),
                  !is.na(experience_years), !is.na(age_sq)) %>%
    dplyr::mutate(
      educ_c     = educ_years - mean(educ_years, na.rm = TRUE),
      exp_c      = experience_years - mean(experience_years, na.rm = TRUE),
      educ_exp_c = educ_c * exp_c)
  if (estimator == "OLS") {
    fml <- ln_hourly_wage_real ~ educ_years + educ_exp_c +
      experience_years + age_sq + married_num | wave
  } else {  # within-person FE
    fml <- ln_hourly_wage_real ~ educ_exp_c +
      experience_years + age_sq + married_num | person_id + wave
  }
  m <- tryCatch(feols(fml, data = d, cluster = ~person_id, notes = FALSE),
                error = function(e) NULL)
  if (is.null(m) || !"educ_exp_c" %in% names(coef(m))) return(NULL)
  ct <- summary(m)$coeftable["educ_exp_c", ]
  list(b2 = unname(ct[1]), se = unname(ct[2]), p = unname(ct[4]),
       n = nrow(d), n_persons = dplyr::n_distinct(d$person_id))
}

# men: pre / post  ×  {OLS, FE};  women: pre / post  ×  {OLS, FE}
o1_cells <- tibble::tribble(
  ~sex,     ~period,           ~data_period,        ~frame,
  "Men",    "Pre",  "Pre-fatherhood",  "men",
  "Men",    "Post", "Post-fatherhood", "men",
  "Women",  "Pre",  "Pre-motherhood",  "women",
  "Women",  "Post", "Post-motherhood", "women")

.frame_of <- function(f) if (f == "men") men_trans_data else women_trans_data

o1_rows <- list()
for (i in seq_len(nrow(o1_cells))) {
  cc <- o1_cells[i, ]
  dd <- .frame_of(cc$frame) %>% dplyr::filter(period == cc$data_period)
  for (est in c("OLS", "FE")) {
    r <- .b2_aligned(dd, est)
    o1_rows[[length(o1_rows) + 1]] <- if (is.null(r)) {
      tibble::tibble(sex = cc$sex, period = cc$period, estimator = est,
                     beta2_x100 = NA_real_, se_x100 = NA_real_, p = NA_real_,
                     stars = "", N = NA_integer_, N_persons = NA_integer_)
    } else {
      tibble::tibble(sex = cc$sex, period = cc$period, estimator = est,
                     beta2_x100 = r$b2 * 100, se_x100 = r$se * 100, p = r$p,
                     stars = stars_fn(r$p), N = r$n, N_persons = r$n_persons)
    }
  }
}
A3_aligned <- dplyr::bind_rows(o1_rows)

cat("\n  Like-for-like beta2 (x100) by sex x period x estimator:\n\n")
print(as.data.frame(A3_aligned %>%
                      dplyr::mutate(beta2_x100 = round(beta2_x100, 4),
                                    se_x100 = round(se_x100, 4), p = round(p, 4))),
      row.names = FALSE)

# Aligned gender gaps: men-minus-women within the SAME estimator and period
.gap_row <- function(est, per) {
  m <- A3_aligned %>% dplyr::filter(sex == "Men",   estimator == est, period == per)
  w <- A3_aligned %>% dplyr::filter(sex == "Women", estimator == est, period == per)
  if (nrow(m) == 0 || nrow(w) == 0 || is.na(m$beta2_x100) || is.na(w$beta2_x100))
    return(NULL)
  gap <- m$beta2_x100 - w$beta2_x100
  se  <- sqrt(m$se_x100^2 + w$se_x100^2)
  z   <- gap / se
  p   <- 2 * pnorm(-abs(z))
  tibble::tibble(estimator = est, period = per,
                 his_b2 = m$beta2_x100, her_b2 = w$beta2_x100,
                 gap = gap, se = se, z = z, p = p, stars = stars_fn(p))
}
A3_gaps <- dplyr::bind_rows(
  .gap_row("OLS", "Pre"),  .gap_row("OLS", "Post"),
  .gap_row("FE",  "Pre"),  .gap_row("FE",  "Post"))

if (nrow(A3_gaps) > 0) {
  cat("\n  Aligned gender gap (his - hers), same estimator both sides:\n\n")
  print(as.data.frame(A3_gaps %>%
                        dplyr::mutate(dplyr::across(c(his_b2, her_b2, gap, se, z, p), ~ round(.x, 4)))),
        row.names = FALSE)
}

write.csv(A3_aligned, "se_ledger_A3_aligned.csv", row.names = FALSE)
write.csv(A3_gaps,    "se_ledger_A3_aligned_gaps.csv", row.names = FALSE)
cat("\n  Wrote se_ledger_A3_aligned.csv and se_ledger_A3_aligned_gaps.csv\n")

# Ready-to-paste table-note caveat (the non-computational resolution path).
a3_note <- paste0(
  "Note on estimator alignment (Table A3). In the columns labelled \"OLS\", ",
  "men's and women's education--experience cross-partials are both estimated ",
  "by within-period pooled OLS with wave fixed effects; in the columns ",
  "labelled \"FE\", both are estimated by within-person fixed effects with ",
  "person and wave fixed effects. Within each estimator the his-minus-hers ",
  "contrast is therefore like-for-like, so the gender difference cannot be ",
  "attributed to a difference in estimator between the two sides. The two ",
  "estimators identify the cross-partial off different variation -- pooled ",
  "OLS off cross-person differences within the pre (or post) cell, the ",
  "fixed-effects estimator off within-person change -- and the gender ",
  "contrast is small under both, so the conclusion does not depend on the ",
  "estimator choice. Standard errors are clustered at the individual level; ",
  "the gender-gap row reports the his-minus-hers difference with its ",
  "two-sided p-value. Sample: parenthood transitioners in the wage sample, ",
  "HILDA Waves 12--24.")
writeLines(strwrap(a3_note, width = 90), "tableA3_estimator_note.txt")
cat("  Wrote tableA3_estimator_note.txt (table-note caveat, paste-ready)\n")

# =============================================================================
# O-2.  SELECTION-CORRECTED WITHIN-COUPLE betaX GENDER GAP
# =============================================================================
# Issue (revision_plan.md O-2): the within-couple betaX comparison is defined
# only on DUAL-EARNER couples (both partners in the wage sample). Couples where
# the mother exits employment entirely -- arguably where the penalty is largest
# -- drop out. The pre-birth check already in the paper answers the first-order
# version; this implements the deeper correction the plan describes: a selection
# correction applied to the betaX gender-gap comparison ITSELF.
#
# Two complementary treatments are produced:
#   (a) Heckman-type IMR control. A couple-wave dual-earner selection probit
#       gives an inverse Mills ratio that we add to the linear within-couple
#       betaX his/her regressions; the corrected his-minus-hers gap is compared
#       against the uncorrected gap.
#   (b) Lee (2009) bounds. Under the assumption that tightening selection only
#       moves couples OUT of the dual-earner sample, the differential
#       (her-minus-his) exit wedge gives a trimming fraction; trimming the
#       relevant tail of the wives' wage distribution bounds the within-couple
#       betaX gap. This mirrors CS_event_study.R's event-study Lee bounds,
#       transposed onto the betaX comparison.
#
# All betaX regressions use the equation-(4) Mincer form: experience as the
# return of interest, schooling among the controls, NO education x experience
# interaction in the control set (consistent with tab:mincer_gap).
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n")
cat("  O-2. Selection-corrected within-couple betaX gender gap\n")
cat(strrep("=", 78), "\n")

# Couple-wave universe where BOTH partners are at least employment-observed.
# A couple-wave is "dual-earner" (in-sample for betaX) iff both are in the wage
# sample with valid wages; otherwise it is selected out.
cu <- couples %>%
  dplyr::mutate(
    dual_earner = as.integer(his_in_wage == 1 & her_in_wage == 1 &
                               !is.na(his_wage) & !is.na(her_wage)),
    # couple-level selection drivers (exclusion-style: affect dual-earner
    # status through the wife's employment margin, not her wage conditional on
    # working). num_children_under15 carried from the female partner where
    # available; fall back to has_children.
    his_married_num = as.numeric(his_married),
    her_married_num = as.numeric(her_married)) %>%
  dplyr::filter(!is.na(his_educ), !is.na(her_educ),
                !is.na(his_exp),  !is.na(her_exp))

# pull a child-burden driver from panel if present
if ("num_children_under15" %in% names(panel)) {
  kid <- panel %>%
    dplyr::transmute(person_id, wave,
                     kids_u15 = dplyr::coalesce(num_children_under15, 0L)) %>%
    dplyr::distinct(person_id, wave, .keep_all = TRUE)
  cu <- cu %>%
    dplyr::left_join(kid, by = c("person_id" = "person_id", "wave" = "wave")) %>%
    dplyr::mutate(kids_u15 = dplyr::coalesce(kids_u15, 0L))
} else {
  cu <- cu %>% dplyr::mutate(kids_u15 = has_children)
}

cat(sprintf("\n  Couple-waves: %s | dual-earner: %s (%.1f%%)\n",
            format(nrow(cu), big.mark = ","),
            format(sum(cu$dual_earner), big.mark = ","),
            100 * mean(cu$dual_earner)))

# --- (a) Heckman-type IMR from a couple-level dual-earner selection probit ----
cat("\n  (a) Heckman-type IMR correction on the within-couple betaX gap\n")

sel_probit <- tryCatch(
  glm(dual_earner ~ her_educ + her_exp + her_age_sq + her_married_num +
        kids_u15 + his_educ + his_exp + his_age_sq + as.factor(wave),
      data = cu, family = binomial("probit")),
  error = function(e) { cat("    probit FAILED:", e$message, "\n"); NULL })

selcorr_tbl <- NULL
if (!is.null(sel_probit)) {
  cu$xb  <- predict(sel_probit, newdata = cu, type = "link")
  cu$imr <- dnorm(cu$xb) / pmax(pnorm(cu$xb), 1e-6)
  
  # dual-earner stacked panel WITH the couple-level IMR carried to both rows,
  # equation-(4) Mincer form (no educ x exp interaction in controls)
  de <- cu %>% dplyr::filter(dual_earner == 1)
  his_d <- de %>% dplyr::transmute(couple_id, wave, has_children, imr,
                                   ln_wage = his_wage, educ = his_educ, exp = his_exp,
                                   hours = his_hours, age_sq = his_age_sq,
                                   married = his_married_num, fulltime = his_fulltime,
                                   is_female = 0L)
  her_d <- de %>% dplyr::transmute(couple_id, wave, has_children, imr,
                                   ln_wage = her_wage, educ = her_educ, exp = her_exp,
                                   hours = her_hours, age_sq = her_age_sq,
                                   married = her_married_num, fulltime = her_fulltime,
                                   is_female = 1L)
  stk <- dplyr::bind_rows(his_d, her_d) %>% dplyr::filter(complete.cases(.))
  
  # linear within-couple betaX, his and her, with vs without IMR.
  # couple + wave FE absorb couple heterogeneity (the "within-couple" object);
  # treatment = own experience; schooling among controls; NO educ x exp term.
  .betaX_lin <- function(d, with_imr) {
    rhs <- "exp + educ + hours + age_sq + married + fulltime"
    if (with_imr) rhs <- paste(rhs, "+ imr")
    fml <- as.formula(paste0("ln_wage ~ ", rhs, " | couple_id + wave"))
    m <- tryCatch(feols(fml, data = d, cluster = ~couple_id, notes = FALSE),
                  error = function(e) NULL)
    if (is.null(m) || !"exp" %in% names(coef(m))) return(NULL)
    ct <- summary(m)$coeftable["exp", ]
    list(b = unname(ct[1]), se = unname(ct[2]), p = unname(ct[4]))
  }
  .gap <- function(his, her) {
    if (is.null(his) || is.null(her)) return(NULL)
    g  <- his$b - her$b
    se <- sqrt(his$se^2 + her$se^2)
    z  <- g / se; p <- 2 * pnorm(-abs(z))
    list(his = his$b, her = her$b, gap = g, se = se, z = z, p = p)
  }
  
  his_raw  <- .betaX_lin(stk %>% dplyr::filter(is_female == 0), FALSE)
  her_raw  <- .betaX_lin(stk %>% dplyr::filter(is_female == 1), FALSE)
  his_corr <- .betaX_lin(stk %>% dplyr::filter(is_female == 0), TRUE)
  her_corr <- .betaX_lin(stk %>% dplyr::filter(is_female == 1), TRUE)
  
  g_raw  <- .gap(his_raw,  her_raw)
  g_corr <- .gap(his_corr, her_corr)
  
  mkrow <- function(tag, g) if (is.null(g)) NULL else tibble::tibble(
    spec = tag,
    his_betaX_x100 = g$his * 100, her_betaX_x100 = g$her * 100,
    gap_x100 = g$gap * 100, se_x100 = g$se * 100,
    z = g$z, p = g$p, stars = stars_fn(g$p))
  selcorr_tbl <- dplyr::bind_rows(
    mkrow("Within-couple betaX, uncorrected", g_raw),
    mkrow("Within-couple betaX, + IMR (Heckman)", g_corr))
  
  if (!is.null(selcorr_tbl)) {
    cat("\n    Within-couple betaX gender gap, raw vs IMR-corrected:\n\n")
    print(as.data.frame(selcorr_tbl %>%
                          dplyr::mutate(dplyr::across(c(his_betaX_x100, her_betaX_x100, gap_x100,
                                                        se_x100, z, p), ~ round(.x, 4)))),
          row.names = FALSE)
    write.csv(selcorr_tbl, "mincer_gap_betaX_selcorr.csv", row.names = FALSE)
    cat("\n    Wrote mincer_gap_betaX_selcorr.csv\n")
    if (!is.null(g_corr)) { rec("bx_gap_w_imr", g_corr$gap * 100)
      rec("bx_gap_w_imr_p", g_corr$p) }
  }
}

# --- (b) Lee (2009) bounds on the within-couple betaX gap ---------------------
# Selection wedge: among couples observed in a wave, the share whose wives are
# NOT in the dual-earner sample but whose husbands are. That differential exit
# is the trimming proportion p. We trim fraction q = p / (employed share of
# wives) from the relevant tail of wives' wages and re-estimate her betaX; the
# his side is left untrimmed (men's dual-earner exit is the reference, as in
# the event-study Lee bounds). Trimming the TOP of wives' wages assumes the
# marginal exiters were high-wage (worst case for the gap -> upper bound on
# the his-minus-hers gap); trimming the BOTTOM is the lower bound.
#
# Inference: a couple-block bootstrap (resampling whole couples with
# replacement) recomputes the entire bound each replication, including the
# data-estimated trimming fraction q, and an Imbens-Manski (2004) confidence
# interval is formed for the partially identified gap. The IM interval covers
# the PARAMETER (not the identified set), so it is the object to quote when
# stating whether the selection-pessimistic endpoint is distinguishable from
# zero. Set PARTO_BOOT_B to change the replication count (default 500).
cat("\n  (b) Lee (2009) trimming bounds on the within-couple betaX gap\n")

lee_betaX_tbl <- NULL
# differential exit: husband in wage sample, wife not, among couple-waves
# where the husband is in the wage sample (so the wedge is the wife's extra exit)
husb_in <- cu %>% dplyr::filter(his_in_wage == 1, !is.na(his_wage))
if (nrow(husb_in) > 100) {
  
  # ---------------------------------------------------------------------------
  # Bound estimator as a pure function of a couple-wave frame `dat` (a slice of
  # `cu`). Re-estimates EVERYTHING that depends on the data, including the
  # trimming fraction q, so the bootstrap propagates uncertainty in q as well
  # as in the regression coefficients. Returns observed gap and [lo, hi] in
  # x100 units, plus q. NULL on failure.
  # ---------------------------------------------------------------------------
  .betaX_gap_cfe <- function(stk) {
    need <- c("ln_wage", "exp", "educ", "hours", "age_sq", "married",
              "fulltime", "is_female", "couple_id", "wave")
    stk <- stk[stats::complete.cases(stk[, need]), , drop = FALSE]
    m <- tryCatch(feols(
      ln_wage ~ exp * is_female + educ + hours + age_sq + married + fulltime |
        couple_id + wave,
      data = stk, cluster = ~couple_id, notes = FALSE),
      error = function(e) NULL)
    if (is.null(m) || !"exp" %in% names(coef(m))) return(NULL)
    b      <- coef(m)
    his_bx <- unname(b["exp"])
    her_bx <- unname(b["exp"]) +
      (if ("exp:is_female" %in% names(b)) unname(b["exp:is_female"]) else 0)
    list(his = his_bx, her = her_bx, gap = his_bx - her_bx)
  }
  
  .lee_bound <- function(dat) {
    hi_set <- dat %>% dplyr::filter(his_in_wage == 1, !is.na(his_wage))
    if (nrow(hi_set) < 50) return(NULL)
    p_exit <- mean(hi_set$her_in_wage == 0 | is.na(hi_set$her_wage), na.rm = TRUE)
    es_w   <- mean(hi_set$her_in_wage == 1 & !is.na(hi_set$her_wage), na.rm = TRUE)
    qq     <- if (es_w > 0) min(p_exit / es_w, 0.999) else 0
    
    de <- dat %>% dplyr::filter(dual_earner == 1)
    if (nrow(de) < 50) return(NULL)
    his_d <- de %>% dplyr::transmute(couple_id, wave, is_female = 0L,
                                     ln_wage = his_wage, exp = his_exp, educ = his_educ,
                                     hours = his_hours, age_sq = his_age_sq,
                                     married = his_married_num, fulltime = his_fulltime,
                                     trim_wage = NA_real_)
    her_d <- de %>% dplyr::transmute(couple_id, wave, is_female = 1L,
                                     ln_wage = her_wage, exp = her_exp, educ = her_educ,
                                     hours = her_hours, age_sq = her_age_sq,
                                     married = her_married_num, fulltime = her_fulltime,
                                     trim_wage = her_wage)
    
    base_obs <- .betaX_gap_cfe(dplyr::bind_rows(his_d, her_d))
    if (is.null(base_obs)) return(NULL)
    
    her_valid <- her_d %>% dplyr::filter(!is.na(trim_wage))
    if (qq > 0 && nrow(her_valid) > 50) {
      lo_cut <- quantile(her_valid$trim_wage, 1 - qq, names = FALSE, type = 7)
      hi_cut <- quantile(her_valid$trim_wage,      qq, names = FALSE, type = 7)
      g_top <- .betaX_gap_cfe(dplyr::bind_rows(
        his_d, her_d %>% dplyr::filter(trim_wage <= lo_cut)))
      g_bot <- .betaX_gap_cfe(dplyr::bind_rows(
        his_d, her_d %>% dplyr::filter(trim_wage >= hi_cut)))
    } else {
      g_top <- base_obs; g_bot <- base_obs
    }
    go <- base_obs$gap * 100
    ga <- if (!is.null(g_top)) g_top$gap * 100 else go
    gb <- if (!is.null(g_bot)) g_bot$gap * 100 else go
    list(observed = go,
         lo = min(go, ga, gb, na.rm = TRUE),
         hi = max(go, ga, gb, na.rm = TRUE),
         q  = qq,
         his = base_obs$his * 100, her = base_obs$her * 100)
  }
  
  pt <- .lee_bound(cu)
  if (is.null(pt)) {
    cat("    Point bound failed; Lee block skipped.\n")
  } else {
    q <- pt$q
    cat(sprintf("    Differential wife exit -> trim q=%.3f\n", q))
    gap_obs <- pt$observed; lo_gap <- pt$lo; hi_gap <- pt$hi
    his_b <- pt$his / 100; her_b_obs <- pt$her / 100
    
    # -------------------------------------------------------------------------
    # Couple-block bootstrap: resample whole couples with replacement (the
    # cluster of clustering), recompute the bound each replication, and form an
    # Imbens-Manski (2004) confidence interval for the partially identified
    # gap. The IM interval widens each side of the point bound by a one-sided
    # critical value scaled by that side's bootstrap SE, with the CM constant
    # chosen so coverage is for the PARAMETER, not the identified set.
    # -------------------------------------------------------------------------
    B <- as.integer(Sys.getenv("PARTO_BOOT_B", "500"))
    if (is.na(B) || B < 50) B <- 500L
    set.seed(42)
    couple_ids <- unique(cu$couple_id)
    cu_by <- split(seq_len(nrow(cu)), cu$couple_id)   # row indices per couple
    
    boot_lo <- numeric(0); boot_hi <- numeric(0)
    cat(sprintf("    Couple-block bootstrap: B=%d replications\n", B))
    for (b in seq_len(B)) {
      drawn <- sample(couple_ids, length(couple_ids), replace = TRUE)
      idx   <- unlist(cu_by[drawn], use.names = FALSE)
      rb <- tryCatch(.lee_bound(cu[idx, , drop = FALSE]),
                     error = function(e) NULL)
      if (!is.null(rb)) { boot_lo <- c(boot_lo, rb$lo); boot_hi <- c(boot_hi, rb$hi) }
    }
    
    im_ci <- c(NA_real_, NA_real_); se_lo <- NA_real_; se_hi <- NA_real_
    n_ok <- length(boot_lo)
    if (n_ok >= 50) {
      se_lo <- sd(boot_lo); se_hi <- sd(boot_hi)
      delta <- hi_gap - lo_gap
      # Imbens-Manski critical value C solves
      #   Phi( C + delta / max(se_lo, se_hi) ) - Phi( -C ) = 0.95
      target <- 0.95
      cm_root <- function(C)
        pnorm(C + delta / max(se_lo, se_hi, 1e-8)) - pnorm(-C) - target
      C <- tryCatch(uniroot(cm_root, c(0, 10))$root, error = function(e) qnorm(0.975))
      im_ci <- c(lo_gap - C * se_lo, hi_gap + C * se_hi)
      cat(sprintf("    Bootstrap SE(lower)=%.3f  SE(upper)=%.3f  C_IM=%.3f  (%d/%d ok)\n",
                  se_lo, se_hi, C, n_ok, B))
    } else {
      cat(sprintf("    Bootstrap produced too few valid replications (%d); CI omitted.\n", n_ok))
    }
    
    lee_betaX_tbl <- tibble::tibble(
      quantity = "Within-couple betaX gender gap (his - hers), x100",
      observed = round(gap_obs, 4),
      lee_lower = round(lo_gap, 4),
      lee_upper = round(hi_gap, 4),
      im_ci_lower = round(im_ci[1], 4),
      im_ci_upper = round(im_ci[2], 4),
      boot_se_lower = round(se_lo, 4),
      boot_se_upper = round(se_hi, 4),
      trim_q = round(q, 4),
      his_betaX_x100 = round(his_b * 100, 4),
      her_betaX_obs_x100 = round(her_b_obs * 100, 4))
    
    cat("\n    Lee bounds + Imbens-Manski 95% CI on the within-couple betaX gap:\n\n")
    print(as.data.frame(lee_betaX_tbl), row.names = FALSE)
    write.csv(lee_betaX_tbl, "betaX_lee_bounds.csv", row.names = FALSE)
    cat("\n    Wrote betaX_lee_bounds.csv\n")
    rec("bx_gap_w_lee_lo", lo_gap); rec("bx_gap_w_lee_hi", hi_gap)
    rec("bx_gap_w_lee_ci_lo", im_ci[1]); rec("bx_gap_w_lee_ci_hi", im_ci[2])
  }
} else {
  cat("    Too few husband-in-sample couple-waves for Lee bounds; skipped.\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n", strrep("=", 78), "\n")
cat("  PART O COMPLETE\n")
cat(strrep("=", 78), "\n")
cat("  O-1  se_ledger_A3_aligned.csv, se_ledger_A3_aligned_gaps.csv,\n")
cat("       tableA3_estimator_note.txt\n")
cat("  O-2  mincer_gap_betaX_selcorr.csv, betaX_lee_bounds.csv\n")
cat("\n  Both items are appendix / robustness; no headline number changes.\n")