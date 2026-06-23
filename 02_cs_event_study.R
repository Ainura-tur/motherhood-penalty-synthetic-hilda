# =====================================================================
# cs_event_study.R
# Heterogeneity-robust staggered event study of the child penalty on
# WAGES and EXPERIENCE, via Sun & Abraham (2021), fixest::sunab.
# (Used in place of the Callaway-Sant'Anna `did` engine, which fails in
#  its internal standardization on this extract; same design and
#  identifying assumption, robust to heterogeneous timing.)
#
# Run after within_couples_MASTER_hh.R with the broad person-wave frame
# `panel` in the workspace. Treatment timing = first birth, defined as in
# the master's section 6.1 (the wave the absorbing ever_parent goes 0->1).
# Never-treated = always-childless; left-censored parents are dropped.
# Endpoints are BINNED at |event time| = K (default 8); the binned cells
# pool the thinning long-horizon observations.
# =====================================================================
stopifnot(requireNamespace("fixest", quietly = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(fixest) })
# --- Namespace protection (avoid MASS/stats masking dplyr verbs) ----
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

K <- 8L   # bin |event time| >= K into the endpoints; widen/narrow as needed

if (!exists("panel")) {
  if (file.exists("hilda_panel_data_W12_W24_slim.rds"))
    panel <- readRDS("hilda_panel_data_W12_W24_slim.rds") else
      stop("`panel` not in memory; run the master first, or provide the slim rds.")
}
stopifnot(all(c("person_id","wave","female","ever_parent",
                "ln_hourly_wage_real","experience_years") %in% names(panel)))
has_lw  <- "longitudinal_weight" %in% names(panel)
has_emp <- "employed" %in% names(panel)   # extensive-margin (employment) outcome

# --- first-birth timing: wave ever_parent first turns 0 -> 1 --------
g_tbl <- panel %>% group_by(person_id) %>% arrange(wave) %>%
  summarise(ep_first = first(ever_parent), ep_last = last(ever_parent),
            trans = { d_ep <- diff(ever_parent)            # NA-safe 0->1 jump
            idx  <- which(!is.na(d_ep) & d_ep > 0)
            if (length(idx)) wave[idx[1] + 1] else NA_integer_ },
            .groups = "drop") %>%
  mutate(cohort = case_when(
    !is.na(trans)                ~ as.double(trans),
    ep_first == 0 & ep_last == 0 ~ 10000,
    TRUE                         ~ NA_real_))

cs <- panel %>% inner_join(select(g_tbl, person_id, cohort), by = "person_id") %>%
  filter(!is.na(cohort))

clean_num <- function(x) { x <- unclass(x); attributes(x) <- NULL
if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x))) }
keep_num <- c("wave","female","ln_hourly_wage_real","experience_years",
              if (has_emp) "employed",
              if (has_lw) "longitudinal_weight")
for (v in keep_num) cs[[v]] <- clean_num(cs[[v]])
cs <- cs %>% mutate(id = as.integer(factor(person_id)),
                    wave = as.double(wave), cohort = as.double(cohort))
cs <- as.data.frame(cs)
cat(sprintf("Sample: %d persons (%d never-treated, %d treated); waves %d-%d; K=%d\n",
            n_distinct(cs$id), n_distinct(cs$id[cs$cohort == 10000]),
            n_distinct(cs$id[cs$cohort < 10000]),
            as.integer(min(cs$wave)), as.integer(max(cs$wave)), K))

# bin relative time >= K and <= -K into the endpoints (pools thinning cells)
BIN <- stats::setNames(list((-200):(-K), K:200), c(as.character(-K), as.character(K)))

# --- estimator: Sun-Abraham event study (binned endpoints) ----------
fit <- function(d, y, w = NULL) {
  rhs <- "sunab(cohort, wave, bin.rel = BIN)"
  f <- as.formula(paste0(y, " ~ ", rhs, " | id + wave"))
  m <- tryCatch(feols(f, data = d, weights = w, cluster = ~ id),
                error = function(e) NULL)
  if (is.null(m))                                  # fallback: no binning
    m <- feols(as.formula(paste0(y, " ~ sunab(cohort, wave) | id + wave")),
               data = d, weights = w, cluster = ~ id)
  m
}
es_tab <- function(m) {
  ct <- coef(m); s <- se(m); nm <- names(ct)
  rt <- suppressWarnings(as.numeric(sub(".*::(-?\\d+).*", "\\1", nm)))
  k <- !is.na(rt) & grepl("::", nm)
  data.frame(event_time = rt[k], att = round(unname(ct[k]), 4),
             se = round(unname(s[k]), 4))[order(rt[k]), ]
}
att_overall <- function(m) tryCatch({
  a <- summary(m, agg = "att")$coeftable
  a[grep("ATT|att", rownames(a))[1], 1:2] }, error = function(e) c(NA, NA))

# ---- headline: UNWEIGHTED -----------------------------------------
res <- list()
outcomes <- c("ln_hourly_wage_real", "experience_years",
              if (has_emp) "employed")
for (sx in c("women","men")) {
  d <- filter(cs, female == ifelse(sx == "women", 1, 0))
  for (y in outcomes) {
    m <- fit(d, y); tag <- paste(sx, y); a <- att_overall(m)
    cat("\n== ", tag, " ==\n", sep = "")
    cat(sprintf("overall post-birth ATT = %.4f (SE %.4f)\n", a[1], a[2]))
    et <- es_tab(m); print(et); res[[tag]] <- et
  }
}
gap <- function(a, b) { e <- intersect(a$event_time, b$event_time)
data.frame(event_time = e,
           gap = round(a$att[match(e, a$event_time)] - b$att[match(e, b$event_time)], 4),
           se  = round(sqrt(a$se[match(e, a$event_time)]^2 +
                              b$se[match(e, b$event_time)]^2), 4)) }
cat("\n== child penalty (women - men), ln hourly wage ==\n")
print(gap(res[["women ln_hourly_wage_real"]], res[["men ln_hourly_wage_real"]]))
cat("\n== differential experience (women - men) ==\n")
print(gap(res[["women experience_years"]], res[["men experience_years"]]))
if (has_emp) {
  cat("\n== child penalty (women - men), employment (extensive margin) ==\n")
  print(gap(res[["women employed"]], res[["men employed"]]))
}
saveRDS(res, "cs_event_study_results.rds")

# ---- attrition robustness: LONGITUDINAL-WEIGHTED wage path ---------
if (has_lw) {
  cat("\n#### attrition robustness: longitudinal-weighted (lnwtrp>0) wage ATT ####\n")
  for (sx in c("women","men")) {
    d <- filter(cs, female == ifelse(sx == "women", 1, 0),
                longitudinal_weight > 0, !is.na(longitudinal_weight))
    a <- att_overall(fit(d, "ln_hourly_wage_real", w = ~longitudinal_weight))
    cat(sprintf("  %-6s weighted wage ATT = %.4f (SE %.4f)\n", sx, a[1], a[2]))
  }
  cat("  (compare with the unweighted ATTs above; large divergence flags attrition.)\n")
} else cat("\n(no longitudinal_weight column found; skipping attrition robustness)\n")

cat("\nNotes:\n",
    sprintf(" * endpoints binned at |event time| >= %d; those bins pool the thinning\n", K),
    "   long-horizon cells, whose per-period counts and precision fall off.\n",
    " * event_time < 0 coefficients test parallel pre-trends; read them in the\n",
    "   women-minus-men gap, where life-cycle trends cancel and the near-birth\n",
    "   pre-period is flat.\n",
    " * wages are observed only for the employed, so the post-birth wage path is a\n",
    "   lower bound; the employment (extensive-margin) event study above bounds the\n",
    "   selection that drives this, women's employment falls sharply at first birth.\n",
    " * headline is unweighted (a causal ATT); the longitudinal-weighted rerun is an\n",
    "   attrition check only, not the main estimate.\n")

# =====================================================================
# When sourced from MASTER_hh.R with options(cs_estimation_only = TRUE),
# estimation stops above; the figure below is built only on a standalone
# run (or by RUN_FIGURES_hh.R).
# =====================================================================
if (!isTRUE(getOption("cs_estimation_only"))) {
  
  # =====================================================================
  # fig_event_study.R
  # Builds figures/fig_event_study.pdf, the women-minus-men event-study
  # paths for log wage and experience, from cs_event_study_results.rds
  # (saved by cs_event_study.R). Endpoints binned at |event time| = K to
  # match Table~\ref{tab:event_study}.
  # =====================================================================
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
  # --- Namespace protection (avoid MASS/stats masking dplyr verbs) ----
  select <- dplyr::select;       filter    <- dplyr::filter
  mutate <- dplyr::mutate;       slice     <- dplyr::slice
  recode <- dplyr::recode;       rename    <- dplyr::rename
  summarise <- dplyr::summarise; summarize <- dplyr::summarize
  arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag
  
  K <- 8L
  res <- readRDS("cs_event_study_results.rds")
  
  
  gap <- function(w, m) {                      # women minus men, with binning at +/-K
    e <- intersect(w$event_time, m$event_time)
    g  <- w$att[match(e, w$event_time)] - m$att[match(e, m$event_time)]
    v  <- w$se[match(e, w$event_time)]^2 + m$se[match(e, m$event_time)]^2
    eb <- pmax(pmin(e, K), -K)                 # bin tails into +/-K
    tapply(seq_along(eb), eb, function(ix) {   # inverse-variance pool within bin
      wgt <- 1 / v[ix]
      c(att = sum(g[ix] * wgt) / sum(wgt), se = sqrt(1 / sum(wgt)))
    }) -> agg
    data.frame(event_time = as.numeric(names(agg)),
               att = sapply(agg, `[`, "att"), se = sapply(agg, `[`, "se"))
  }
  
  dat <- bind_rows(
    transform(gap(res[["women ln_hourly_wage_real"]], res[["men ln_hourly_wage_real"]]),
              outcome = "Log hourly wage"),
    transform(gap(res[["women experience_years"]],    res[["men experience_years"]]),
              outcome = "Experience (years)"),
    if (all(c("women employed","men employed") %in% names(res)))
      transform(gap(res[["women employed"]], res[["men employed"]]),
                outcome = "Employment (pp)")) %>%
    mutate(lo = att - 1.96 * se, hi = att + 1.96 * se,
           outcome = factor(outcome, levels = c("Log hourly wage", "Experience (years)",
                                                "Employment (pp)")))
  
  p <- ggplot(dat, aes(event_time, att)) +
    facet_wrap(~outcome, scales = "free_y") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = 2, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line() + geom_point(size = 1.4) +
    labs(x = "Years from first birth (women minus men)", y = "ATT") +
    theme_minimal(base_size = 11)
  
  dir.create("figures", showWarnings = FALSE)
  ggsave("figures/fig_event_study.pdf", p, width = 12, height = 3.6)
  cat("wrote figures/fig_event_study.pdf\n")
  
}  # end if (!isTRUE(getOption("cs_estimation_only")))