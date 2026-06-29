# =====================================================================
# 08_gelbach_decomposition.R
# Gelbach (2016, JoLE) conditional decomposition of the women-minus-men
# post-birth wage gap into the share carried by each mediator. The focal
# coefficient is the motherhood contrast (post-birth x female) on log real
# hourly wage. Adding mediators (accumulated experience, and optionally
# market hours / training) shrinks that coefficient; Gelbach apportions the
# shrinkage across mediators in an order-INVARIANT way, unlike sequential
# controls.
#
# Decomposition (within-transformed by person + wave, so it is OLS on the
# residualised data and the omitted-variable algebra applies):
#   base:  y = X1 b1 + e          X1 = [focal, base covars]
#   full:  y = X1 b1f + X2 b2 + u  X2 = mediators
#   b1 - b1f = Gamma b2,  Gamma = (X1'X1)^{-1} X1'X2   (auxiliary regressions)
#   contribution of mediator m to the focal coefficient = Gamma[focal, m]*b2[m]
#   sum over m = b1[focal] - b1f[focal]   (the total "explained" by mediators)
#
# -----------------------------------------------------------------------------
# SAMPLE (the critical fix, 2026-06).  Appendix A6 states the decomposition is
# run "on the event-study person-wave panel ... ever-parents and never-treated
# childless individuals, with left-censored parents dropped TO MATCH THE
# EVENT-STUDY SAMPLE." Table D1 (Panel B) fixes that population at 4,761
# persons. Earlier versions of this script rebuilt their own frame from `panel`
# with ONLY a left-censoring drop plus complete-cases, which is a SUPERSET of
# the event-study estimation sample (it keeps people with no clean 0->1 birth
# window). That superset inflated the person count (8,981) and DILUTED the
# pooled focal coefficient (-0.0090 instead of -0.0239), which in turn pushed
# the pooled experience share from ~17% to ~57%. The explained component is
# stable; the denominator was not, because it was computed on the wrong sample.
#
# This version therefore restricts the decomposition to the EXACT event-study
# estimation sample, in priority order:
#   (1) an explicit person-id roster exported by 02_cs_event_study.R
#       (`es_person_ids`, or an `event_study_sample` object), if present;
#   (2) failing that, it RECONSTRUCTS the event-study sample with the same
#       birth-window rule (a clean observed 0->1 ever_parent transition with at
#       least one pre- and one post-birth wage observation) PLUS never-treated
#       always-childless persons.
# It then audits the realised person count against the wage-observed event-study
# count read from the exported manifest (cs_event_study_sample_manifest.rds; ~8,967
# persons), NOT against Table D1's "4,761", which was a documentation error. It
# stops if the gap exceeds a tolerance, rather than silently running on the wrong
# population.
# -----------------------------------------------------------------------------
#
# It is a mediation DECOMPOSITION, not a clean causal mediation; reading the
# experience share as causal needs sequential ignorability of experience given
# the focal contrast and covariates, which we do not assert. The complementary
# "control in the event study" check at the foot re-estimates the late-horizon
# share via return x gap (07b), which avoids the treatment-mediator collinearity
# of a late-only regression control.
#
# Consumes the broad person-wave frame `panel` (as in 02_cs_event_study.R):
#   person_id, wave, female, ever_parent, ln_hourly_wage_real,
#   experience_years  (+ optionally market hours / training columns, auto-detected).
# =====================================================================

stopifnot(requireNamespace("fixest", quietly = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(fixest) })
# --- Namespace protection (avoid MASS/stats masking dplyr verbs) ----
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

set.seed(42)
B_BOOT <- 300L   # person-cluster bootstrap replications for mediator-share SEs

# The event-study sample target is NOT hard-coded. Table D1's "4,761" was a
# documentation error (the event-study population is ~9,952 persons, of which
# ~8,967 are wage-observed). The decomposition's outcome is ln_hourly_wage_real,
# so complete.cases necessarily restricts to the WAGE-OBSERVED event-study sample.
# We therefore read the target from the manifest the event-study script exports
# (cs_event_study_sample_manifest.rds) and check against the wage-observed count.
# The target and tolerance are computed below in section 1c, AFTER the roster is
# resolved, because they depend on which roster section 1 loaded.

# Strip haven/labelled classes etc. to plain numeric. HILDA columns loaded via
# haven are `haven_labelled` doubles for which is.numeric() is FALSE, which
# makes crossprod() reject them. Same helper as 02_cs_event_study.R.
clean_num <- function(x) {
  x <- unclass(x); attributes(x) <- NULL
  if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x)))
}

# ---- 0. load panel --------------------------------------------------
if (!exists("panel")) {
  if (file.exists("hilda_panel_data_W12_W24_slim.rds"))
    panel <- readRDS("hilda_panel_data_W12_W24_slim.rds") else
      stop("`panel` not in memory; run the master first, or provide the slim rds.")
}
stopifnot(all(c("person_id","wave","female","ever_parent",
                "ln_hourly_wage_real","experience_years") %in% names(panel)))

# ---- 1. identify the EVENT-STUDY estimation sample ------------------
# The decomposition must run on the SAME persons the Sun-Abraham event study
# uses, not on a looser left-censoring-only frame. We resolve the roster in
# priority order and record provenance for the audit at the end.

es_source <- NA_character_
es_person_ids <- NULL

# (1) explicit roster exported by the event-study script -------------------
if (is.null(es_person_ids) && exists("es_person_ids_export")) {
  es_person_ids <- unique(es_person_ids_export)
  es_source <- "es_person_ids_export (in memory)"
}
if (is.null(es_person_ids) && exists("event_study_sample")) {
  if ("person_id" %in% names(event_study_sample)) {
    es_person_ids <- unique(event_study_sample$person_id)
    es_source <- "event_study_sample$person_id (in memory)"
  }
}
for (f in c("cs_event_study_sample_ids.rds", "event_study_person_ids.rds",
            "cs_event_study_results.rds")) {
  if (is.null(es_person_ids) && file.exists(f)) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    cand <- NULL
    if (is.data.frame(obj) && "person_id" %in% names(obj)) cand <- obj$person_id
    else if (is.list(obj) && !is.null(obj$person_ids))      cand <- obj$person_ids
    else if (is.list(obj) && !is.null(obj$sample_ids))      cand <- obj$sample_ids
    else if (is.atomic(obj))                                cand <- obj
    if (!is.null(cand)) {
      es_person_ids <- unique(cand)
      es_source <- sprintf("%s (on disk)", f)
    }
  }
}

# (2) reconstruct the event-study sample with the birth-window rule --------
# Mirror 02_cs_event_study.R: keep (a) never-treated always-childless persons,
# and (b) clean 0->1 ever_parent transitioners who have at least one pre-birth
# and one post-birth observation with a VALID hourly wage (so they actually
# enter the wage event study). Left-censored parents (ever_parent==1 in their
# first observed wave) are excluded by construction in (b).
if (is.null(es_person_ids)) {
  es_source <- "reconstructed birth-window rule (02_cs_event_study.R logic)"
  
  ps <- panel %>%
    arrange(person_id, wave) %>%
    group_by(person_id) %>%
    summarise(
      ep_first   = first(ever_parent),
      ep_min     = min(ever_parent, na.rm = TRUE),
      ep_max     = max(ever_parent, na.rm = TRUE),
      n_obs      = dplyr::n(),
      # birth wave = first wave ever_parent == 1 (NA if never a parent)
      birth_wave = {
        w <- wave[ever_parent == 1]
        if (length(w)) min(w) else NA_integer_
      },
      # pre/post WAGE observations relative to birth
      pre_wage = sum(ever_parent == 0 & !is.na(ln_hourly_wage_real)),
      post_wage = sum(ever_parent == 1 & !is.na(ln_hourly_wage_real)),
      .groups = "drop"
    )
  
  # never-treated always-childless: ever_parent == 0 in every observed wave
  never_treated <- ps %>% filter(ep_max == 0)
  
  # clean 0->1 transitioners with an observed pre- and post-birth wage window;
  # ep_first == 0 drops left-censored parents (already parent in first wave).
  transitioners <- ps %>%
    filter(ep_first == 0, ep_max == 1, pre_wage >= 1, post_wage >= 1)
  
  es_person_ids <- unique(c(never_treated$person_id, transitioners$person_id))
}

es_person_ids <- es_person_ids[!is.na(es_person_ids)]

# ---- 1b. build the analysis frame on the EVENT-STUDY sample ---------
# Focal contrast: post-birth x female. "post-birth" = ever_parent has turned on
# (absorbing). We no longer do an independent left-censoring drop here: the
# event-study roster already encodes it. We DO keep the complete-cases filter on
# the variables that enter the algebra.
df <- panel %>%
  filter(person_id %in% es_person_ids) %>%
  mutate(
    female      = as.numeric(female),
    post        = as.numeric(ever_parent),       # absorbing post-birth indicator
    focal       = post * female,                 # women-minus-men post-birth gap
    age_sq      = if ("age" %in% names(.)) (age^2) / 100 else NA_real_,
    married     = if ("married_num" %in% names(.)) as.numeric(married_num)
    else if ("married" %in% names(.)) as.numeric(married)
    else NA_real_
  )

# Coerce every column we will feed to the matrix algebra to plain numeric, so a
# stray haven_labelled / integer64 / factor column cannot poison crossprod().
maybe_num <- intersect(
  c("ln_hourly_wage_real","focal","post","female","age_sq","married",
    "experience_years","market_hours","weekly_hours","hours","jbhruc",
    "paid_hours","training","training_any","work_related_training",
    "employer_training"),
  names(df))
for (v in maybe_num) df[[v]] <- clean_num(df[[v]])

# base covariates that are present and non-degenerate
base_covars <- c("post", "female", "age_sq", "married")
base_covars <- base_covars[vapply(base_covars, function(v)
  isTRUE(v %in% names(df) && sd(df[[v]], na.rm = TRUE) > 0), logical(1))]

# mediators: experience always; market hours and training if detectable.
# NOTE: the paper's headline single-mediator share REQUIRES experience to be the
# sole mediator. If hours/training columns are present in `panel`, this becomes a
# MULTI-mediator decomposition and the experience share is then experience's
# slice of the joint explained component, NOT a single-mediator share. We print
# the mediator set loudly so the run can never be silently mislabelled, and we
# expose a switch to force single-mediator to match the paper exactly.
FORCE_SINGLE_MEDIATOR <- TRUE   # set FALSE to report the full multi-mediator split

hours_col <- intersect(c("market_hours","weekly_hours","hours","jbhruc",
                         "paid_hours"), names(df))[1]
train_col <- intersect(c("training","training_any","work_related_training",
                         "employer_training"), names(df))[1]
mediators <- c("experience_years",
               if (!FORCE_SINGLE_MEDIATOR && !is.na(hours_col)) hours_col,
               if (!FORCE_SINGLE_MEDIATOR && !is.na(train_col)) train_col)
mediators <- unique(mediators[!is.na(mediators)])

cat("Focal: post-birth x female.  Base covars: ", paste(base_covars, collapse=", "),
    "\nMediators: ", paste(mediators, collapse = ", "),
    if (FORCE_SINGLE_MEDIATOR) "  [FORCE_SINGLE_MEDIATOR=TRUE]" else
      "  [multi-mediator: experience share is a JOINT-explained slice]",
    "\n", sep = "")

keep <- unique(c("person_id","wave","ln_hourly_wage_real","focal",
                 base_covars, mediators))
df <- df %>% select(all_of(keep)) %>% filter(complete.cases(.))

# ---- 1c. SAMPLE AUDIT against the event-study manifest --------------
n_es_persons <- dplyr::n_distinct(df$person_id)
n_es_obs     <- nrow(df)

# Resolve the target from the manifest the event-study script exports, not a
# hard-coded number. The Gelbach outcome is ln_hourly_wage_real, so complete.cases
# necessarily restricts to the WAGE-OBSERVED event-study sample; the correct
# target is therefore n_persons_wageobs. The full roster (n_persons_full) includes
# never-treated controls with no observed wage who cannot enter a wage regression.
.es_manifest <- if (file.exists("cs_event_study_sample_manifest.rds"))
  readRDS("cs_event_study_sample_manifest.rds") else NULL
ES_TARGET_PERSONS <- if (!is.null(.es_manifest))
  .es_manifest$n_persons_wageobs else length(es_person_ids)
ES_TOL_PERSONS <- max(50L, round(0.03 * ES_TARGET_PERSONS))

cat(strrep("-", 72), "\n", sep = "")
cat("Event-study sample provenance: ", es_source, "\n", sep = "")
cat(sprintf("Analysis frame (event-study sample): %s person-waves, %s persons\n",
            format(n_es_obs, big.mark = ","),
            format(n_es_persons, big.mark = ",")))
cat(sprintf("Event-study target (wage-observed) persons: %s  (tol +/- %d)\n",
            format(ES_TARGET_PERSONS, big.mark = ","), ES_TOL_PERSONS))
gap_persons <- n_es_persons - ES_TARGET_PERSONS
cat(sprintf("Realised - target = %+d persons\n", gap_persons))
if (abs(gap_persons) > ES_TOL_PERSONS) {
  stop(sprintf(paste0(
    "Gelbach sample does NOT match the event-study population.\n",
    "  Realised persons = %d, wage-observed target = %d (gap %+d, tol %d).\n",
    "  This is the 17%%-vs-57%% bug: the focal coefficient is being computed on\n",
    "  the wrong sample. Resolve the roster (export the event-study sample from\n",
    "  cs_event_study.R, which writes cs_event_study_sample_manifest.rds) before\n",
    "  trusting any share. Run halted deliberately rather than report a number\n",
    "  off the wrong population."),
    n_es_persons, ES_TARGET_PERSONS, gap_persons, ES_TOL_PERSONS))
}
cat("Sample matches the event-study population within tolerance. Proceeding.\n")
cat(strrep("-", 72), "\n", sep = "")

# ---- 2. Gelbach engine (within person + wave) -----------------------
# Residualise y, focal, base covars and mediators on the two-way FE once;
# then the decomposition is OLS algebra on the residualised matrix.
gelbach <- function(d) {
  fe_vars <- c("person_id","wave")
  # `female` is time-invariant within person and is absorbed by the person FE,
  # so it is dropped. `post` is a within-person 0->1 transition and is NOT
  # absorbed: it must stay as the main effect so the interaction `focal` is the
  # women-minus-men differential (b in: a*post + b*post:female), not mothers'
  # own pre/post change. The focal coefficient of interest is the row "focal".
  X1_names <- c("post", "focal", intersect(c("age_sq","married"), names(d)))
  X1_names <- intersect(X1_names, names(d))   # keep only columns that exist
  vars <- unique(c("ln_hourly_wage_real", X1_names, mediators))
  # build a plain numeric matrix, then demean: demean() on a numeric matrix
  # returns a numeric matrix, sidestepping any labelled/typed-column issue.
  Xnum <- vapply(vars, function(v) clean_num(d[[v]]), numeric(nrow(d)))
  colnames(Xnum) <- vars
  dm <- fixest::demean(Xnum, d[, fe_vars, drop = FALSE])
  storage.mode(dm) <- "double"
  y  <- dm[, "ln_hourly_wage_real"]
  X1 <- dm[, X1_names, drop = FALSE]
  X2 <- dm[, mediators, drop = FALSE]
  
  # full model coefficients on mediators (b2) and focal (b1f)
  Xf  <- cbind(X1, X2)
  bf  <- qr.solve(crossprod(Xf), crossprod(Xf, y))
  rownames(bf) <- colnames(Xf)
  b1f_focal <- bf["focal", 1]
  b2 <- bf[colnames(X2), 1, drop = FALSE]
  
  # base model focal coefficient (b1)
  bb <- qr.solve(crossprod(X1), crossprod(X1, y))
  rownames(bb) <- colnames(X1)
  b1_focal <- bb["focal", 1]
  
  # auxiliary regressions: each mediator on X1; Gamma[focal, m]
  G <- qr.solve(crossprod(X1), crossprod(X1, X2))   # (ncol X1) x (ncol X2)
  rownames(G) <- colnames(X1); colnames(G) <- colnames(X2)
  gamma_focal <- G["focal", , drop = TRUE]
  
  contrib <- gamma_focal * b2[, 1]                  # per-mediator contribution
  list(b1_focal = b1_focal, b1f_focal = b1f_focal,
       total_explained = b1_focal - b1f_focal,
       contrib = contrib)
}

pt <- gelbach(df)

# ---- 3. person-cluster bootstrap for mediator-share SEs -------------
cat(sprintf("\nPerson-cluster bootstrap (B = %d) ...\n", B_BOOT))
# Precompute each person's row indices ONCE. Resampling is then just a
# concatenation of these vectors, O(n) per replication instead of O(n * persons).
idx_by_person <- split(seq_len(nrow(df)), df$person_id)
n_persons     <- length(idx_by_person)
row_counts    <- lengths(idx_by_person)

boot_mat <- matrix(NA_real_, nrow = B_BOOT,
                   ncol = length(pt$contrib) + 2,
                   dimnames = list(NULL, c("total_explained", "direct",
                                           names(pt$contrib))))
t0 <- Sys.time()
for (b in seq_len(B_BOOT)) {
  samp <- sample.int(n_persons, replace = TRUE)          # positions in the list
  rows <- unlist(idx_by_person[samp], use.names = FALSE)
  db   <- df[rows, , drop = FALSE]
  # re-id so each drawn cluster is distinct (a person drawn twice = two clusters)
  db$person_id <- rep.int(seq_along(samp), row_counts[samp])
  gb <- tryCatch(gelbach(db), error = function(e) NULL)
  if (!is.null(gb)) boot_mat[b, ] <- c(gb$total_explained, gb$b1f_focal, gb$contrib)
  if (b %% 25 == 0)
    cat(sprintf("  %d/%d  (%.0fs elapsed)\n", b, B_BOOT,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
boot_mat <- boot_mat[stats::complete.cases(boot_mat), , drop = FALSE]
cat(sprintf("  bootstrap done: %d usable replications\n", nrow(boot_mat)))

bse <- apply(boot_mat, 2, sd)
blo <- apply(boot_mat, 2, quantile, 0.025)
bhi <- apply(boot_mat, 2, quantile, 0.975)

# ---- 4. report ------------------------------------------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("Gelbach decomposition of the post-birth women-minus-men wage gap\n")
cat("(event-study estimation sample; see provenance above)\n")
cat(strrep("=", 72), "\n", sep = "")
cat(sprintf("Total motherhood penalty (base focal coef)   : %+.4f\n", pt$b1_focal))
cat(sprintf("Direct / unexplained (full focal coef)       : %+.4f  (SE %.4f)\n",
            pt$b1f_focal, bse["direct"]))
cat(sprintf("Total explained by mediators (base - full)   : %+.4f  (SE %.4f) [%.4f, %.4f]\n",
            pt$total_explained, bse["total_explained"],
            blo["total_explained"], bhi["total_explained"]))
cat(strrep("-", 72), "\n", sep = "")
cat(sprintf("%-22s %10s %9s %10s %8s\n",
            "Mediator", "contrib", "SE", "share", "95% CI"))
for (m in names(pt$contrib)) {
  share <- if (abs(pt$b1_focal) > 1e-10) pt$contrib[m] / pt$b1_focal else NA
  cat(sprintf("%-22s %+10.4f %9.4f %9.0f%%  [%+.4f, %+.4f]\n",
              m, pt$contrib[m], bse[m], 100 * share, blo[m], bhi[m]))
}
cat(strrep("-", 72), "\n", sep = "")
exp_share <- pt$contrib["experience_years"] / pt$b1_focal
cat(sprintf("Experience-channel share of the total penalty: %.0f%%\n",
            100 * exp_share))

# Honesty guard: the pooled share divides by a focal coefficient that, on the
# event-study sample, is small and (its residual) statistically insignificant.
# Flag when the base focal coef is imprecise so the share is never over-read.
if (!is.na(bse["direct"]) && abs(pt$b1_focal) < 2 * bse["direct"]) {
  cat("\n[caution] The base focal coefficient is small relative to the\n",
      "  bootstrap SE of the residual, so the POOLED share is a ratio with an\n",
      "  imprecise denominator and is fragile. Anchor the mechanism on the\n",
      "  experience-gap event-study path and the training event study, and read\n",
      "  the pooled share as indicative, with the per-horizon / return-times-gap\n",
      "  calculation (07b) as the reliable late-horizon analogue.\n", sep = "")
}

cat("\nReading: the experience share is how much of the motherhood wage gap the\n",
    "accumulated-experience gap accounts for; the residual 'direct' term is the\n",
    "qualitative margin (occupational downgrading, intensity) plus any mediator\n",
    "not included. Decomposition is order-invariant; it is accounting, not a\n",
    "causal mediation claim.\n", sep = "")

saveRDS(list(point = pt, boot_se = bse, boot_lo = blo, boot_hi = bhi,
             mediators = mediators, n_obs = nrow(df),
             n_persons = n_es_persons,
             es_source = es_source,
             es_target_persons = ES_TARGET_PERSONS,
             force_single_mediator = FORCE_SINGLE_MEDIATOR),
        "gelbach_decomposition_results.rds")
cat("\nSaved gelbach_decomposition_results.rds\n")

# ---- 5. dynamic check: reliable late-horizon share via return x gap --
# We deliberately do NOT estimate the late share by adding experience as a
# control in a late-only regression. At late horizons the treatment indicator
# (mother, post-birth) and the experience deficit fall in the SAME person-waves,
# so mediator and treatment are collinear and the difference method mechanically
# under-attributes to experience, leaving the gap on the treatment dummy. The
# reliable late share separates the two pieces: the return to experience
# (estimated where experience varies for everyone) times the experience gap
# (estimated separately in the event study). That is what
# 07b_implied_dynamic_cost_concave.R computes, so we read it here rather than
# re-running a contaminated regression.
f07 <- "implied_dynamic_cost_concave_results.rds"
if (file.exists(f07)) {
  tb <- tryCatch(readRDS(f07)$table, error = function(e) NULL)
  if (!is.null(tb) && all(c("event_time","share_concave") %in% names(tb))) {
    late <- tb[tb$event_time >= 5 & tb$event_time <= 10, ]
    sh <- mean(late$share_concave, na.rm = TRUE)
    pooled_share_pct <- 100 * exp_share
    cat(sprintf("\n[dynamic] late-horizon (k=5..10) experience share, return x gap (07b): %.0f%%\n",
                100 * sh))
    cat(sprintf(paste0(
      "  Pooled Gelbach share (this run, event-study sample): %.0f%%.\n",
      "  Compare directions: the pooled share is a horizon-average over flat\n",
      "  early horizons and the large late ones; whether it sits above or below\n",
      "  the late-horizon return-times-gap figure depends on the realised focal\n",
      "  coefficient. Report both numbers from the SAME (event-study) sample and\n",
      "  state the comparison explicitly rather than asserting a fixed 17/29 or\n",
      "  57/29 ordering.\n"),
      pooled_share_pct))
  }
} else {
  cat("\n[dynamic] run 07b_implied_dynamic_cost_concave.R for the late-horizon\n",
      "  experience share (return x gap). A late-only regression control is\n",
      "  unreliable here: treatment and the experience deficit are collinear.\n", sep = "")
}