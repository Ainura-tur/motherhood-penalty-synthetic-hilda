# =====================================================================
# training_event_study.R
# Heterogeneity-robust staggered event study of TRAINING around the
# first birth, via Sun & Abraham (2021), fixest::sunab.
#
# Design mirrors CS_event_study.R exactly: treatment = wave the absorbing
# ever_parent indicator turns 0 -> 1; never-treated = always-childless;
# left-censored parents dropped; endpoints binned at |event time| = K;
# women and men estimated separately; women-minus-men gap reported.
#
# Outcomes (linear probability / counts):
#   training_paid_time   employer-time job-related training (jttopwt+jttpewt>0)
#   training_own_time    own-time job-related training      (jttopot+jttpeot>0)
#   total_training_days  days of training, last 12 months   (jttdays, NA->0)
#   total_courses        number of courses, last 12 months  (jttnum,  NA->0)
#   jttpewt_any / jttpeot_any   PLACE-OF-EMPLOYMENT-ONLY variants, used
#       automatically when present. The amended loader (HILDA_LOADIND_IV.R)
#       now derives and carries both. HILDA coding confirmed: 1 = Yes,
#       2 = No, available from wave 7; all indicators use == 1 (the former
#       "> 0" rule misclassified 2 = No as training).
#
# Two estimands per outcome, run for both sexes:
#   (a) UNCONDITIONAL: non-employed coded 0 (the loader's coalesce), so the
#       path is the TOTAL effect including employment exit at birth. This is
#       the policy-relevant incidence path.
#   (b) EMPLOYED-ONLY: conditional on employment, subject to the same
#       post-birth selection caveat as wages (positive selection => the
#       conditional training drop is a lower bound). Comparison object for
#       Blundell, Costa-Dias, Goll & Meghir (2021).
#
# Run after MASTER_hh.R with the broad person-wave frame `panel` in the
# workspace, or standalone with the slim rds. Training variables exist
# from the W15-W24 loader; W12-W14 rows (if present) are dropped from the
# training sample with a message.
# Writes training_event_study_results.rds; figure only on standalone run
# (or set options(training_estimation_only = TRUE) when sourcing).
# =====================================================================
stopifnot(requireNamespace("fixest", quietly = TRUE))
suppressPackageStartupMessages({ library(dplyr); library(fixest) })
# --- Namespace protection (avoid MASS/stats masking dplyr verbs) ----
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

K <- 8L   # bin |event time| >= K into the endpoints

# Work on a LOCAL frame `tr_data`; never overwrite the master's `panel`.
# The master's slim W12-W24 panel has no training columns, so when sourced
# after the master we fall through to the training rds automatically.
.tr_names <- c("training_paid_time","training_own_time","total_training_days",
               "total_courses","jttpewt_any","jttpeot_any")
tr_data <- NULL
if (exists("panel") && any(.tr_names %in% names(panel))) {
  tr_data <- panel
} else if (file.exists("hilda_panel_data_training.rds")) {
  if (exists("panel"))
    cat("In-memory `panel` carries no training variables;",
        "using hilda_panel_data_training.rds instead.\n")
  tr_data <- readRDS("hilda_panel_data_training.rds")
} else if (file.exists("hilda_panel_data_W12_W24_slim.rds")) {
  tr_data <- readRDS("hilda_panel_data_W12_W24_slim.rds")
} else
  stop("No training data: run HILDA_LOADIND_MASTER.R first to build ",
       "hilda_panel_data_training.rds.")
stopifnot(all(c("person_id","wave","female","ever_parent","employed") %in% names(tr_data)))

# --- outcome set: use what the panel actually carries ----------------
tr_core <- c("training_paid_time","training_own_time",
             "total_training_days","total_courses")
tr_pe   <- c("jttpewt_any","jttpeot_any")           # pe-only, if loader amended
tr_outcomes <- c(intersect(tr_core, names(tr_data)),
                 intersect(tr_pe,   names(tr_data)))
if (!length(tr_outcomes))
  stop("No training variables found in `panel`. Check the loader's final select.")
cat("Training outcomes found:", paste(tr_outcomes, collapse = ", "), "\n")
if (!all(tr_pe %in% names(tr_data)))
  cat("(place-of-employment-only items not in panel; using pooled indicators.\n",
      " To add them, see the loader amendment in this script's header.)\n")

# --- training sample: waves where the training module is observed ----
# The W15-W24 loader supplies the training vars; if the in-memory panel is
# W12-W24, the early waves carry no training info. Keep waves where at least
# one training outcome is non-degenerate.
tr_wave_ok <- tr_data %>% group_by(wave) %>%
  summarise(ok = any(sapply(across(all_of(tr_outcomes)), function(x) any(!is.na(x) & x != 0))),
            .groups = "drop")
tr_waves <- tr_wave_ok$wave[tr_wave_ok$ok]
cat(sprintf("Training module observed in waves %d-%d; dropping other waves.\n",
            min(tr_waves), max(tr_waves)))
tpanel <- tr_data %>% filter(wave %in% tr_waves)

# --- first-birth timing: identical to CS_event_study.R ---------------
# (cohort defined on the FULL panel so timing is not censored by the
#  training-wave restriction)
g_tbl <- tr_data %>% group_by(person_id) %>% arrange(wave) %>%
  summarise(ep_first = first(ever_parent), ep_last = last(ever_parent),
            trans = { d_ep <- diff(ever_parent)
            idx  <- which(!is.na(d_ep) & d_ep > 0)
            if (length(idx)) wave[idx[1] + 1] else NA_integer_ },
            .groups = "drop") %>%
  mutate(cohort = case_when(
    !is.na(trans)                ~ as.double(trans),
    ep_first == 0 & ep_last == 0 ~ 10000,
    TRUE                         ~ NA_real_))

cs <- tpanel %>% inner_join(select(g_tbl, person_id, cohort), by = "person_id") %>%
  filter(!is.na(cohort))

clean_num <- function(x) { x <- unclass(x); attributes(x) <- NULL
if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x))) }
for (v in c("wave","female","employed", tr_outcomes)) cs[[v]] <- clean_num(cs[[v]])
cs <- cs %>% mutate(id = as.integer(factor(person_id)),
                    wave = as.double(wave), cohort = as.double(cohort))
cs <- as.data.frame(cs)
cat(sprintf("Sample: %d persons (%d never-treated, %d treated); waves %d-%d; K=%d\n",
            n_distinct(cs$id), n_distinct(cs$id[cs$cohort == 10000]),
            n_distinct(cs$id[cs$cohort < 10000]),
            as.integer(min(cs$wave)), as.integer(max(cs$wave)), K))

BIN <- stats::setNames(list((-200):(-K), K:200), c(as.character(-K), as.character(K)))

fit <- function(d, y) {
  f <- as.formula(paste0(y, " ~ sunab(cohort, wave, bin.rel = BIN) | id + wave"))
  m <- tryCatch(feols(f, data = d, cluster = ~ id), error = function(e) NULL)
  if (is.null(m))
    m <- feols(as.formula(paste0(y, " ~ sunab(cohort, wave) | id + wave")),
               data = d, cluster = ~ id)
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
gap <- function(a, b) { e <- intersect(a$event_time, b$event_time)
data.frame(event_time = e,
           gap = round(a$att[match(e, a$event_time)] - b$att[match(e, b$event_time)], 4),
           se  = round(sqrt(a$se[match(e, a$event_time)]^2 +
                              b$se[match(e, b$event_time)]^2), 4)) }

# --- estimation: unconditional and employed-only ---------------------
res <- list()
for (cond in c("unconditional","employed_only")) {
  cat("\n", strrep("#", 70), "\n#  ", toupper(cond), "\n", strrep("#", 70), "\n", sep = "")
  base <- if (cond == "employed_only") filter(cs, employed == 1) else cs
  for (sx in c("women","men")) {
    d <- filter(base, female == ifelse(sx == "women", 1, 0))
    for (y in tr_outcomes) {
      m <- tryCatch(fit(d, y), error = function(e) NULL)
      if (is.null(m)) { cat(sprintf("  %s %s %s: not estimable\n", cond, sx, y)); next }
      tag <- paste(cond, sx, y); a <- att_overall(m)
      cat(sprintf("\n== %s ==\noverall post-birth ATT = %.4f (SE %.4f)\n", tag, a[1], a[2]))
      et <- es_tab(m); print(et); res[[tag]] <- et
    }
  }
  for (y in tr_outcomes) {
    w <- res[[paste(cond, "women", y)]]; m <- res[[paste(cond, "men", y)]]
    if (is.null(w) || is.null(m)) next
    cat(sprintf("\n== child penalty (women - men), %s [%s] ==\n", y, cond))
    print(gap(w, m))
  }
}
saveRDS(res, "training_event_study_results.rds")
cat("\nSaved training_event_study_results.rds\n")

cat("\nNotes:\n",
    " * UNCONDITIONAL paths code the non-employed as 0 (the loader's coalesce),\n",
    "   so they bundle employment exit with the training margin: the total\n",
    "   incidence effect of the first birth on training receipt.\n",
    " * EMPLOYED-ONLY paths isolate the training margin among the employed but\n",
    "   inherit post-birth positive selection (cf. the wage lower-bound logic);\n",
    "   a conditional drop is therefore a lower bound on the true drop.\n",
    " * Read pre-birth coefficients as the parallel-trends test, preferably in\n",
    "   the women-minus-men gap where life-cycle training trends cancel.\n",
    " * Interpretation hook: Blundell et al. (2021) argue training partially\n",
    "   rebuilds post-birth earnings capacity; a post-birth training DEFICIT for\n",
    "   mothers means the rebuild channel is under-used precisely where the\n",
    "   experience channel bites, reinforcing the returnship policy lever.\n")

# =====================================================================
# Figure: women-minus-men training incidence paths (standalone run only)
# =====================================================================
if (!isTRUE(getOption("training_estimation_only"))) {
  suppressPackageStartupMessages({ library(ggplot2) })
  pick <- intersect(c("training_paid_time","training_own_time"), tr_outcomes)
  dat <- dplyr::bind_rows(lapply(pick, function(y) {
    g <- gap(res[[paste("unconditional","women",y)]],
             res[[paste("unconditional","men",y)]])
    transform(g, outcome = y)
  })) %>%
    mutate(lo = gap - 1.96 * se, hi = gap + 1.96 * se,
           event_time = pmax(pmin(event_time, K), -K))
  p <- ggplot(dat, aes(event_time, gap)) +
    facet_wrap(~outcome) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = 2, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line() + geom_point(size = 1.4) +
    labs(x = "Years from first birth (women minus men)",
         y = "ATT, training incidence") +
    theme_minimal(base_size = 11)
  dir.create("figures", showWarnings = FALSE)
  ggsave("figures/fig_training_event_study.pdf", p, width = 9, height = 3.6)
  cat("wrote figures/fig_training_event_study.pdf\n")
}