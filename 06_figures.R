# ============================================================================
# RUN_FIGURES.R  --  Generates and saves EVERY paper figure (except CDAG)
#                    into figures/, under the filenames the .tex expects.
#
# RUN AFTER within_couples_MASTER_ANALYSISprof.R, IN THE SAME R SESSION.
# The figure code below was lifted verbatim from the master and depends on the
# master's in-memory data objects (wage, panel, stacked, couples, dml_cells_rf,
# event_hw, ce_hw, beta2/transition/robustness/afb frames, ...). Those must be
# present in the workspace (they are, once the master has run / the A-D state
# files are loaded).
#
# Paper figures written (matches \\includegraphics in cost_hh_spec.tex):
#   fig2_event_study_beta2, fig3_employment_birth,
#   fig6_robustness_forest  -> fig7_robustness_forest,
#   figC1_gender_gap_vanishes -> fig5C1_gender_gap_vanishes   (via .fig_map)
#   fig_iv_coefplot, fig_hw_within_fe, fig_event_study        (built directly)
# All other legacy figure blocks (fig1/4/5/7/8, C2-C6, HW1-HW5) still run for
# their printed verification tables, but their saves are skipped by the
# masked ggsave()/save_fig() because they are not in .fig_map.
#   CDAG.pdf -> generated separately (not here)
# ============================================================================

## Do NOT clear the workspace: the figure code reuses the master's objects.
suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(patchwork); library(scales); library(grid)
  library(lmtest); library(sandwich); library(fixest)
})

# Namespace protection (avoid car/MASS/plyr masking dplyr verbs)
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

# Resolve the project folder the same way as MASTER_hh.R: script location
# when source()d or run via Rscript, then the known Windows folder.
.this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile),
                       error = function(e) NA_character_)
if (is.na(.this_file)) {
  .arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.arg) > 0)
    .this_file <- tryCatch(normalizePath(sub("^--file=", "", .arg[1])),
                           error = function(e) NA_character_)
}
if (!is.na(.this_file)) {
  setwd(dirname(.this_file))
} else if (dir.exists("C:/Users/tur277/OneDrive - CSIRO/Desktop/WG")) {
  setwd("C:/Users/tur277/OneDrive - CSIRO/Desktop/WG")
}
cat("Working directory:", getwd(), "\n")

# ── Console log: mirror all output to a timestamped file ──────────────────────
# The verification tables printed below (Fig 3 employment table, three-way
# beta2 table, etc.) are needed to populate/check the paper, so they are
# written to run_logs/06_figures_<timestamp>.log as well as the console.
# Opened AFTER the setwd() above so the log lands in the project folder
# regardless of the directory the script was launched from.
logdir <- file.path(getwd(), "run_logs")
if (!dir.exists(logdir)) dir.create(logdir, recursive = TRUE)
.figlog_path <- file.path(logdir, format(Sys.time(), "06_figures_%Y%m%d_%H%M%S.log"))
.figlog_con  <- file(.figlog_path, open = "wt")
sink(.figlog_con, split = TRUE)
sink(.figlog_con, type = "message")
close_fig_log <- function() {
  # message sinks do not stack: a single sink(type = "message") restores stderr.
  # (Do NOT loop on sink.number(type="message"): it returns 2 (= stderr) when
  # no diversion is active, so a "> 0" loop never terminates.)
  try(sink(type = "message"), silent = TRUE)
  while (sink.number() > 0) sink()
  try(close(.figlog_con), silent = TRUE)
  message(sprintf("Console log written to: %s", .figlog_path))
}
options(error = function() {
  traceback(2); close_fig_log()
  if (!interactive()) quit(status = 1, save = "no")
})
cat(sprintf("Logging console output to: %s\n", .figlog_path))


# ----------------------------------------------------------------------------
# Restore the master's PART A-D workspace if the data objects are not already
# live in this session. The master writes full save.image() checkpoints; the
# PART D checkpoint holds wage, panel, stacked, couples, dml_cells_rf, event_hw,
# ce_hw, etc. The figure code below recomputes its own frames (wbw, ev_b2, rob,
# afb_data) from these base objects, so the A-D checkpoint is all that is needed.
# ----------------------------------------------------------------------------
if (!exists("wage") || !exists("panel") || !exists("stacked") || !exists("couples") || !exists("dml_cells_rf")) {
  .ckpts <- c("master_state_after_partD.RData", "master_state_after_partC.RData",
              "master_state_after_partB.RData", "master_state_after_partA.RData")
  .hit <- .ckpts[file.exists(.ckpts)]
  if (length(.hit)) {
    load(.hit[1], envir = .GlobalEnv)
    cat(sprintf("  restored workspace from %s\n", .hit[1]))
  } else {
    cat("[WARN] no master_state_after_part*.RData found in wd; figures will fail without workspace data.\n")
  }
}

# Namespace protection -- RE-APPLIED here, after the checkpoint load, because the
# restored workspace reinstates the master's (incomplete) masking and can leave
# car::recode shadowing dplyr::recode.
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

# Output directory: prefer config.R's fig_dir (output/figures). When run from
# run_all.R this is already set; the guard keeps a standalone run working and
# falls back to a local output/figures/ if config.R is unavailable.
if (!exists("fig_dir")) {
  if (file.exists("config.R")) source("config.R") else
    fig_dir <- file.path(getwd(), "output", "figures")
}
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
outdir <- fig_dir

# ----------------------------------------------------------------------------
# House-style helpers (verbatim from the master)
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# Save redirection: every figure goes to figures/<paper-name>.pdf.
# These wrappers mask ggsave()/save_fig() so the lifted blocks run unchanged.
# PNG outputs and non-paper figures are skipped automatically.
# ----------------------------------------------------------------------------
.fig_map <- c(
  # Only the figures used in the current paper are mapped. Unmapped figure
  # blocks still run but their save is skipped (save_one returns NULL).
  # fig_iv_coefplot, fig_hw_within_fe and fig_event_study save directly to
  # figures/ at the end of this script (they bypass this map).
  "fig2_event_study_beta2"    = "fig2_event_study_beta2",
  "fig3_employment_birth"     = "fig3_employment_birth",
  "fig6_robustness_forest"    = "fig7_robustness_forest",
  "figC1_gender_gap_vanishes" = "fig5C1_gender_gap_vanishes"
)

# Render with cairo_pdf when available so unicode glyphs (Delta, minus sign,
# beta, en-dash) embed correctly; fall back to pdf() if cairo is unavailable.
.pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

save_one <- function(plot, base, w, h) {
  paper <- unname(.fig_map[base])
  if (is.na(paper)) return(invisible(NULL))          # not a paper figure -> skip
  ggplot2::ggsave(file.path(outdir, paste0(paper, ".pdf")),
                  plot, width = w, height = h, units = "in", device = .pdf_device)
  cat(sprintf("  figures/%s.pdf\n", paper))
}

# Mask ggsave() used by the lifted blocks (raw calls).
ggsave <- function(filename, plot = ggplot2::last_plot(),
                   width = 13/2.54, height = 8/2.54, ...) {
  if (grepl("\\.png$", filename)) return(invisible(NULL))   # skip PNG
  save_one(plot, sub("\\.pdf$", "", basename(filename)), width, height)
}

# Mask the master's save_fig(p, name, w, h)
save_fig <- function(p, name, w = 13/2.54, h = 8/2.54) save_one(p, name, w, h)

# Neutralise the master's "copy to home dir" loop (figures already in figures/)
file.copy <- function(...) invisible(TRUE)

# ----------------------------------------------------------------------------
# Workspace guard: warn (do not stop) if a key data object is missing
# ----------------------------------------------------------------------------
.need <- c("wage","panel","stacked","couples","dml_cells_rf","event_hw","ce_hw")
.missing <- .need[!vapply(.need, exists, logical(1))]
if (length(.missing))
  cat(sprintf("\n[WARN] missing workspace objects: %s\n  -> run the master in this session first; some figures may error.\n",
              paste(.missing, collapse = ", ")))

cat("\n=== Building paper figures into output/figures/ ===\n")

# ----------------------------------------------------------------------------
# Hoist PART E housework intermediates (event_hw, ce_hw, ...) from the saved
# bundle if they are not already in the session. The master writes them to
# housework_analysis.rds; the figure code below needs event_hw and ce_hw.
# ----------------------------------------------------------------------------
if (!exists("event_hw") || !exists("ce_hw")) {
  if (file.exists("housework_analysis.rds")) {
    .hw_state <- readRDS("housework_analysis.rds")
    for (.nm in names(.hw_state)) if (!exists(.nm)) assign(.nm, .hw_state[[.nm]])
    cat("  loaded event_hw / ce_hw from housework_analysis.rds\n")
  } else {
    cat("[WARN] housework_analysis.rds not found; HW figures (HW2/HW3/HW4) will be skipped.\n")
  }
}

# ----------------------------------------------------------------------------
# Training event study (paper Figure 5). Mirror the CS handling: source
# 03_training_event_study.R estimation-only to (re)build its results rds if it is
# missing, then build the paper-style two-panel figure HERE from the rds, so it
# matches Figure 5 (theme_paper, "A. Employer-funded" / "B. Own-time" panels)
# rather than the script's rough standalone diagnostic. Figure rebuilt each run.
# ----------------------------------------------------------------------------
if (!file.exists("training_event_study_results.rds") && file.exists("03_training_event_study.R")) {
  cat("  training_event_study_results.rds not found; sourcing 03_training_event_study.R (estimation only)...\n")
  options(training_estimation_only = TRUE)   # estimation only; paper-style figure built below
  source("03_training_event_study.R")
  options(training_estimation_only = NULL)
}
if (file.exists("training_event_study_results.rds")) {
  .tr_res <- readRDS("training_event_study_results.rds")
  .tr_K   <- 8L
  .tr_gap <- function(w, m) {                  # women minus men, inverse-variance pooled at +/-K
    e  <- intersect(w$event_time, m$event_time)
    g  <- w$att[match(e, w$event_time)] - m$att[match(e, m$event_time)]
    v  <- w$se[match(e, w$event_time)]^2 + m$se[match(e, m$event_time)]^2
    eb <- pmax(pmin(e, .tr_K), -.tr_K)
    agg <- tapply(seq_along(eb), eb, function(ix) {
      wgt <- 1 / v[ix]
      c(att = sum(g[ix] * wgt) / sum(wgt), se = sqrt(1 / sum(wgt)))
    })
    data.frame(event_time = as.numeric(names(agg)),
               att = sapply(agg, `[`, "att"), se = sapply(agg, `[`, "se"))
  }
  .tr_panels <- c(training_paid_time = "A. Employer-funded training",
                  training_own_time  = "B. Own-time training")
  .tr_have <- names(.tr_panels)[vapply(names(.tr_panels), function(y)
    all(c(paste("unconditional women", y), paste("unconditional men", y)) %in% names(.tr_res)),
    logical(1))]
  if (length(.tr_have)) {
    tr_dat <- dplyr::bind_rows(lapply(.tr_have, function(y)
      transform(.tr_gap(.tr_res[[paste("unconditional women", y)]],
                        .tr_res[[paste("unconditional men",   y)]]),
                outcome = .tr_panels[[y]]))) %>%
      mutate(lo = att - 1.96 * se, hi = att + 1.96 * se,
             outcome = factor(outcome, levels = unname(.tr_panels)))
    p_tr <- ggplot(tr_dat, aes(event_time, att)) +
      facet_wrap(~outcome) +
      geom_hline(yintercept = 0,    linetype = "dashed", colour = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#7570B3", alpha = 0.15) +
      geom_line(colour = "#7570B3", linewidth = 0.9) +
      geom_point(colour = "#7570B3", size = 1.8) +
      labs(x = "Years from first birth", y = "ATT, training incidence (women \u2212 men)") +
      theme_paper
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(file.path(outdir, "fig_training_event_study.pdf"), p_tr,
                    width = 13/2.54, height = 6/2.54, units = "in")
    cat("  figures/fig_training_event_study.pdf (paper style)\n")
  } else {
    cat("[WARN] training outcomes not in results rds; fig_training_event_study.pdf not rebuilt.\n")
  }
} else {
  cat("[WARN] training_event_study_results.rds unavailable; fig_training_event_study.pdf not rebuilt.\n")
}

cat("\n[A] Housework event-study figures (HW2/HW3/HW4) ...\n")
if (exists("event_hw") && exists("ce_hw")) {
  
  # --- Fig HW2: Housework event study (individual) ---
  cat("--- Fig HW2: Housework around birth ---\n")
  
  event_hw_long <- event_hw %>%
    select(t, gender, Housework = mean_hw, Market = mean_market) %>%
    pivot_longer(cols = c(Housework, Market), names_to = "type", values_to = "hours")
  
  fig_hw2 <- ggplot(event_hw_long,
                    aes(x = t, y = hours, colour = type, linetype = type)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_line(linewidth = 1) +
    facet_wrap(~gender, nrow = 1) +
    scale_colour_manual(
      values = c("Market" = "#2166AC", "Housework" = "#B2182B"),
      labels = c("Market" = "Market hours", "Housework" = "Housework hours")
    ) +
    scale_linetype_manual(
      values = c("Market" = "solid", "Housework" = "dashed"),
      labels = c("Market" = "Market hours", "Housework" = "Housework hours")
    ) +
    scale_x_continuous(breaks = -3:4) +
    labs(
      title    = "Market Hours and Housework Hours Around First Birth",
      subtitle = NULL,
      x = "Years relative to first birth", y = "Hours per week",
      colour = NULL, linetype = NULL,
      caption  = "Source: HILDA Waves 12\u201324. Vertical dotted line = birth event."
    ) +
    theme_paper +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.4, "cm"),
          legend.key.size  = unit(0.3, "cm"),
          legend.text      = element_text(size = 5.5),
          strip.text       = element_text(size = 6.5, face = "bold"))
  
  save_fig(fig_hw2, "figHW2_hw_event_study", w = 14/2.54, h = 9/2.54)
  
  
  # --- Fig HW3: Within-couple housework reallocation around birth ---
  cat("--- Fig HW3: Within-couple housework around birth ---\n")
  
  ce_hw_long <- ce_hw %>%
    select(t, his_hw, her_hw) %>%
    pivot_longer(-t, names_to = "who", values_to = "hw") %>%
    mutate(gender = ifelse(who == "his_hw", "Men", "Women"))
  
  p_hw3a <- ggplot(ce_hw_long, aes(x = t, y = hw, colour = gender)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    scale_colour_manual(values = pal_gender) +
    scale_x_continuous(breaks = -3:4) +
    labs(title = "A. Housework hours", y = "Hours/week", x = NULL) +
    theme_paper + theme(legend.position = c(0.22, 0.88),
                        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
                        legend.key.size = unit(0.3, "cm"),
                        legend.text = element_text(size = 5.5))
  
  p_hw3b <- ggplot(ce_hw, aes(x = t, y = hw_gap)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_line(colour = "#7570B3", linewidth = 1) +
    geom_point(colour = "#7570B3", size = 2.5) +
    scale_x_continuous(breaks = -3:4) +
    labs(title = "B. Housework gap (hers \u2212 his)", y = "Hours/week", x = NULL) +
    theme_paper
  
  p_hw3c <- ggplot(ce_hw, aes(x = t, y = her_hw_share)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey60") +
    geom_line(colour = "#B2182B", linewidth = 1) +
    geom_point(colour = "#B2182B", size = 2.5) +
    scale_x_continuous(breaks = -3:4) +
    labs(title = "C. HW share (%)",
         y = "%", x = "Years from birth") +
    theme_paper
  
  fig_hw3 <- (p_hw3a | p_hw3b | p_hw3c) +
    plot_annotation(
      title = "Within-Couple Housework Reallocation Around First Birth",
      subtitle = NULL,
      caption = "Source: HILDA Waves 12\u201324. Matched heterosexual couples with observed birth. The housework gap widens at birth and does not reverse.",
      theme = theme(plot.title = element_text(size = 7.5, face = "bold"),
                    plot.subtitle = element_text(size = 6, colour = "grey40"),
                    plot.caption = element_text(size = 6, colour = "grey50")))
  
  save_fig(fig_hw3, "figHW3_couple_hw_birth", w = 18/2.54, h = 10/2.54)
  
  
  # --- Fig HW4: Complete time budget around birth (stacked area) ---
  cat("--- Fig HW4: Complete time budget around birth ---\n")
  
  budget_event <- event_hw %>%
    select(t, gender, market = mean_market, housework = mean_hw) %>%
    pivot_longer(cols = c(market, housework), names_to = "type", values_to = "hours") %>%
    mutate(type = factor(ifelse(type == "market", "Market", "Housework"),
                         levels = c("Market", "Housework")))
  
  fig_hw4 <- ggplot(budget_event, aes(x = t, y = hours, fill = type)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_area(alpha = 0.7, position = "stack") +
    facet_wrap(~gender) +
    scale_fill_manual(values = c("Market" = "#4393C3", "Housework" = "#D6604D")) +
    scale_x_continuous(breaks = -3:4) +
    labs(
      title = "Total Productive Hours Around First Birth: Market + Housework",
      subtitle = NULL,
      x = "Years relative to first birth", y = "Hours per week", fill = NULL,
      caption = "Source: HILDA Waves 12\u201324. Women's market hours drop but housework rises only partially; men show minimal change."
    ) +
    theme_paper +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.key.size  = unit(0.3, "cm"),
      legend.text      = element_text(size = 5.5),
      legend.margin    = margin(2, 4, 2, 4)
    )
  
  save_fig(fig_hw4, "figHW4_total_budget_birth", w = 14/2.54, h = 9/2.54)
  
  
  
} else cat("  SKIPPED HW2/HW3/HW4 (event_hw/ce_hw unavailable)\n")

cat("\n[B] Main + couple figures (fig1-8, C1/C4/C5/C6) ...\n")
# #############################################################################
# FIGURE 1: WAVE-BY-WAVE β₂ BY GROUP
# #############################################################################

# Paper Appendix Figure A1 (fig8_wave_by_wave_beta2.pdf)
# Sample: Wage sample, rolling OLS by wave and group.
cat("--- Figure 1: Wave-by-wave \u03b2\u2082 ---\n")

wbw <- tibble()
for (w in 12:24) {
  for (g in c("childless_men_ever", "fathers_ever",
              "never_mothers", "mothers_ever")) {
    d <- wage %>%
      filter(wave == w, group_ever == g) %>%
      mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
             exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
             educ_exp_c = educ_c * exp_c)
    
    if (nrow(d) < 30) next
    
    m <- tryCatch(
      lm(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
           age_sq + married_num, data = d),
      error = function(e) NULL)
    
    if (!is.null(m) && "educ_exp_c" %in% names(coef(m))) {
      ct <- coeftest(m, vcov = vcovHC(m, "HC1"))
      wbw <- bind_rows(wbw,
                       tibble(wave = w, year = 2000 + w, group = g,
                              beta2 = ct["educ_exp_c", 1], se = ct["educ_exp_c", 2]))
    }
  }
}
wbw$group_label <- factor(label_group(wbw$group),
                          levels = c("Childless men", "Fathers", "Never-mothers", "Mothers (ever)"))

fig1 <- ggplot(wbw, aes(x = year, y = beta2 * 100,
                        colour = group_label, fill = group_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = (beta2 - 1.96*se)*100,
                  ymax = (beta2 + 1.96*se)*100), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = pal4) +
  scale_fill_manual(values = pal4) +
  scale_x_continuous(breaks = 2012:2024) +
  labs(title = expression(paste("Education\u2013Experience Complementarity (",
                                beta[2], ") by Year and Group")),
       subtitle = NULL,
       x = NULL, y = expression(hat(beta)[2] %*% 100),
       caption = "Source: HILDA Waves 12\u201324, ages 30\u201349.") +
  theme_paper +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        legend.position = c(0.5, 0.92),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5.5))

ggsave( file.path(outdir, "fig1_wave_by_wave_beta2.pdf"), fig1, width = 13/2.54, height = 10/2.54, units = "in")
ggsave( file.path(outdir, "fig1_wave_by_wave_beta2.png"), fig1, width = 13/2.54, height = 10/2.54, dpi = 500, units = "in")

# #############################################################################
# FIGURE 2: EVENT STUDY — β₂ AROUND FIRST BIRTH
# #############################################################################

# Paper Figure 1 (fig1_event_study_beta2.pdf): "Men's Complementarity
#   Collapses at First Birth; Women's Reflects Selection"
# Sample: Wage-sample persons with observed first birth during panel.
#   3-wave rolling OLS. (Sample size: re-verify against current extract.)
cat("--- Figure 2: Event study \u03b2\u2082 around birth ---\n")

# First birth event per person
birth_events <- panel %>%
  filter(had_birth_this_wave == 1) %>%
  group_by(person_id) %>%
  summarise(birth_wave = min(wave), .groups = "drop")

event_panel <- wage %>%
  inner_join(birth_events, by = "person_id") %>%
  mutate(t = wave - birth_wave) %>%
  filter(t >= -4, t <= 5)

num_rows    <- format(nrow(event_panel), big.mark = ",")
num_persons <- format(n_distinct(event_panel$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

# Rolling 3-wave window β₂ at each event-time
ev_b2 <- tibble()
for (sex in 0:1) {
  sex_label <- ifelse(sex == 0, "Men", "Women")
  for (tt in -4:5) {
    win <- c(max(-4, tt - 1), tt, min(5, tt + 1))
    d <- event_panel %>%
      filter(female == sex, t %in% win) %>%
      mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
             exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
             educ_exp_c = educ_c * exp_c)
    
    if (nrow(d) < 40) next
    
    m <- tryCatch(
      feols(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
              age_sq + married_num | wave,
            data = d, cluster = ~person_id, notes = FALSE),
      error = function(e) NULL)
    
    if (!is.null(m) && "educ_exp_c" %in% names(coef(m))) {
      s <- summary(m)$coeftable
      ev_b2 <- bind_rows(ev_b2,
                         tibble(t = tt, gender = sex_label, N = nrow(d),
                                beta2 = s["educ_exp_c", 1], se = s["educ_exp_c", 2]))
    }
  }
}

fig2 <- ggplot(ev_b2, aes(x = t, y = beta2 * 100,
                          colour = gender, fill = gender)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
  geom_ribbon(aes(ymin = (beta2 - 1.96*se)*100,
                  ymax = (beta2 + 1.96*se)*100), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = pal_gender) +
  scale_fill_manual(values = pal_gender) +
  scale_x_continuous(breaks = -4:5,
                     labels = c(paste0("t\u2212", 4:1), "Birth", paste0("t+", 1:5))) +
  annotate("text", x = -3.2,
           y = max(ev_b2$beta2 * 100, na.rm = TRUE) + 0.03,
           label = "Pre-birth", size = 2.2, colour = "grey40", vjust = 0) +
  annotate("text", x = 2.5,
           y = max(ev_b2$beta2 * 100, na.rm = TRUE) + 0.03,
           label = "Post-birth", size = 2.2, colour = "grey40", vjust = 0) +
  labs(title = expression(paste("Education\u2013Experience Complementarity (",
                                beta[2], ") Around First Birth")),
       subtitle = NULL,
       x = "Years relative to first birth",
       y = expression(hat(beta)[2] %*% 100),
       caption = "Source: HILDA Waves 12\u201324. Sample restricted to individuals with observed first birth during panel.") +
  theme_paper +
  theme(legend.position = c(0.92, 0.92),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5.5))

ggsave( file.path(outdir, "fig2_event_study_beta2.pdf"), fig2, width = 13/2.54, height = 10/2.54, units = "in")
ggsave( file.path(outdir, "fig2_event_study_beta2.png"), fig2, width = 13/2.54, height = 10/2.54, dpi = 500, units = "in")


# #############################################################################
# FIGURE 3: EMPLOYMENT / FULL-TIME / HOURS AROUND BIRTH
# #############################################################################

# Paper Figure 2 (fig2_employment_birth.pdf): "Women's Employment,
#   Full-Time Rates, and Hours Collapse at Birth; Men's Are Flat"
# Sample: Broad panel (incl. non-employed), genuine transitions only.
#   Expected: ~2,641 genuine transitions (~1,246 men + ~1,395 women).
cat("--- Figure 3: Employment around birth ---\n")

# Derive birth event wave from ever_parent transitions.
# This is more robust than had_birth_this_wave, which may be NA/absent
# depending on how the panel RDS was constructed.
# Identify genuine transitions: person must have been OBSERVED as childless
# (ever_parent == 0) in at least one wave that precedes their first parent wave.
# This avoids the lag(default=0) artefact that flags people who enter the panel
# already as parents (in any wave 15–24) as false transitions.
# This exactly mirrors the transition logic in Part A and should yield ~2,641
# genuine transitions (~1,246 men + ~1,395 women).
birth_events_fig3 <- panel %>%
  filter(!is.na(ever_parent)) %>%
  group_by(person_id) %>%
  filter(
    any(ever_parent == 0L) &
      any(ever_parent == 1L) &
      suppressWarnings(min(wave[ever_parent == 0L])) <
      suppressWarnings(min(wave[ever_parent == 1L]))
  ) %>%
  summarise(birth_wave = min(wave[ever_parent == 1L]), .groups = "drop")

cat(sprintf("  birth_events_fig3: %d genuine within-panel transitions\n",
            nrow(birth_events_fig3)))
cat(sprintf("  (expected ~2,641: ~1,246 men + ~1,395 women from Part A)\n"))

birth_emp <- panel %>%
  inner_join(birth_events_fig3, by = "person_id") %>%
  mutate(t = wave - birth_wave) %>%
  filter(t >= -3, t <= 4) %>%
  group_by(t, female) %>%
  summarise(
    emp   = mean(employed == 1, na.rm = TRUE) * 100,
    ft    = mean(fulltime == 1, na.rm = TRUE) * 100,
    hours = mean(hours_worked_clean, na.rm = TRUE),
    N     = n(),
    .groups = "drop"
  ) %>%
  mutate(gender = ifelse(female == 0, "Men", "Women"))

cat(sprintf("  birth_emp rows: %d\n", nrow(birth_emp)))

# Print table for console verification
if (nrow(birth_emp) > 0) {
  cat(sprintf("\n  %3s  %-6s  %5s  %8s  %8s  %8s\n",
              "t", "Gender", "N", "Emp%", "FT%", "Hours"))
  cat("  ", strrep("-", 47), "\n")
  for (i in 1:nrow(birth_emp)) {
    r <- birth_emp[i, ]
    cat(sprintf("  %+2d   %-6s  %5d  %7.1f%%  %7.1f%%  %7.1f\n",
                r$t, r$gender, r$N, r$emp, r$ft, r$hours))
  }
  cat("\n")
} else {
  cat("  WARNING: birth_emp is empty — Figure 3 will not be generated.\n")
  cat("  Check that panel contains ever_parent transitions.\n")
}

if (nrow(birth_emp) > 0) {
  p3a <- ggplot(birth_emp, aes(x = t, colour = gender)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_line(aes(y = emp), linewidth = 1) +
    geom_point(aes(y = emp), size = 2) +
    scale_colour_manual(values = pal_gender) +
    labs(title = "A. Employment rate", y = "%", x = NULL) +
    theme_paper + theme(legend.position = "none")
  
  p3b <- ggplot(birth_emp, aes(x = t, colour = gender)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_line(aes(y = ft), linewidth = 1) +
    geom_point(aes(y = ft), size = 2) +
    scale_colour_manual(values = pal_gender) +
    labs(title = "B. Full-time rate", y = "%", x = NULL) +
    theme_paper + theme(legend.position = "none")
  
  p3c <- ggplot(birth_emp, aes(x = t, colour = gender)) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_line(aes(y = hours), linewidth = 1) +
    geom_point(aes(y = hours), size = 2) +
    scale_colour_manual(values = pal_gender) +
    labs(title = "C. Weekly hours", y = "Hours", x = "Years from birth") +
    theme_paper
  
  fig3 <- (p3a | p3b | p3c) +
    plot_annotation(
      title = "Labour Market Outcomes Around First Birth",
      subtitle = "Broad sample (employed + non-employed). Vertical dotted line = birth event.",
      caption = "Source: HILDA Waves 12\u201324.",
      theme = theme(plot.title = element_text(size = 7.5, face = "bold"),
                    plot.subtitle = element_text(size = 6, colour = "grey40"),
                    plot.caption = element_text(size = 6, colour = "grey50")))
  
  tryCatch({
    ggsave(file.path(outdir, "fig3_employment_birth.pdf"), fig3, width = 18/2.54, height = 10/2.54, units = "in")
    ggsave(file.path(outdir, "fig3_employment_birth.png"), fig3, width = 18/2.54, height = 10/2.54, dpi = 500, units = "in")
    for (ext in c(".pdf", ".png")) {
      ff <- paste0("fig3_employment_birth", ext)
      if (file.exists(ff)) file.copy(ff, file.path(outdir, ff), overwrite = TRUE)
    }
    cat("  ✓ fig3_employment_birth saved\n")
  }, error = function(e) cat("  ERROR saving fig3:", e$message, "\n"))
} # end if nrow > 0



# #############################################################################
# FIGURE 4: THREE-WAY WOMEN COEFFICIENT PLOT
# #############################################################################

cat("--- Figure 4: Three-way women coefplot ---\n")

tw_data <- tibble()
for (g in c("never_mothers", "active_mothers", "post_mothers")) {
  d <- wage %>%
    filter(group_3way_women == g) %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c)
  
  for (est_info in list(
    list(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
         age_sq + married_num | wave, "OLS"),
    list(ln_hourly_wage_real ~ educ_exp_c + experience_years +
         age_sq + married_num | person_id + wave, "FE"))) {
    
    m <- tryCatch(feols(est_info[[1]], data = d, cluster = ~person_id,
                        notes = FALSE), error = function(e) NULL)
    if (!is.null(m) && "educ_exp_c" %in% names(coef(m))) {
      s <- summary(m)$coeftable
      tw_data <- bind_rows(tw_data,
                           tibble(group = g, estimator = est_info[[2]], N = nrow(d),
                                  beta2 = s["educ_exp_c", 1], se = s["educ_exp_c", 2]))
    }
  }
}

tw_data$group_label <- factor(
  case_when(tw_data$group == "never_mothers"  ~ "Never-mothers",
            tw_data$group == "active_mothers" ~ "Active mothers",
            tw_data$group == "post_mothers"   ~ "Post-mothers"),
  levels = c("Never-mothers", "Active mothers", "Post-mothers"))

cat("\n  Three-way women \u03b2\u2082 estimates:\n")
cat(sprintf("  %-18s  %-5s  %8s  %10s\n", "Group", "Est.", "N", "b2 x100"))
cat("  ", strrep("-", 50), "\n")
for (i in 1:nrow(tw_data)) {
  r <- tw_data[i, ]
  p_val <- 2 * pnorm(-abs(r$beta2 / r$se))
  cat(sprintf("  %-18s  %-5s  %8s  %+.4f%%%s\n",
              as.character(r$group_label), r$estimator,
              format(r$N, big.mark = ","),
              r$beta2 * 100, stars_fn(p_val)))
}

fig4 <- ggplot(tw_data,
               aes(x = group_label, y = beta2 * 100, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = (beta2 - 1.96*se)*100,
                      ymax = (beta2 + 1.96*se)*100),
                  position = position_dodge(width = 0.5),
                  size = 0.7, linewidth = 0.9) +
  scale_colour_manual(values = c("OLS" = "#2166AC", "FE" = "#B2182B")) +
  labs(
    title = expression(paste("Three-Way Women Split: ", beta[2],
                             " (OLS vs Individual FE)")),
    subtitle = "Post-mothers = ever-parent with no current dependent children. 95% CI.",
    x = NULL, y = expression(hat(beta)[2] %*% 100),
    colour = NULL,
    caption = "Source: HILDA Waves 12\u201324, ages 30\u201349. Absorbing definition."
  ) +
  theme_paper +
  theme(legend.position = c(0.5, 0.92),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5.5))

ggsave( file.path(outdir, "fig4_threeway_coefplot.pdf"), fig4, width = 13/2.54, height = 8/2.54, units = "in")
ggsave( file.path(outdir, "fig4_threeway_coefplot.png"), fig4, width = 13/2.54, height = 8/2.54, dpi = 500, units = "in")


# #############################################################################
# FIGURE 5: EMPLOYMENT RATES BY GROUP ACROSS WAVES
# #############################################################################

cat("--- Figure 5: Employment across waves ---\n")

emp_wave <- panel %>%
  filter(!is.na(group_ever)) %>%
  group_by(year, group_ever) %>%
  summarise(emp = mean(employed == 1, na.rm = TRUE) * 100,
            N = n(), .groups = "drop") %>%
  mutate(group_label = factor(label_group(group_ever),
                              levels = c("Childless men", "Fathers", "Never-mothers", "Mothers (ever)")))


cat("\n  Employment rate by group and year:\n")
cat(sprintf("  %4s  %-20s  %7s  %6s\n", "Year", "Group", "N", "Emp%"))
cat("  ", strrep("-", 45), "\n")
for (i in 1:nrow(emp_wave)) {
  r <- emp_wave[i, ]
  cat(sprintf("  %4d  %-20s  %7d  %5.1f%%\n",
              r$year, as.character(r$group_label), r$N, r$emp))
}
fig5 <- ggplot(emp_wave, aes(x = year, y = emp, colour = group_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = pal4) +
  scale_x_continuous(breaks = 2012:2024) +
  scale_y_continuous(limits = c(60, 100)) +
  labs(title = "Employment Rate by Parenthood Status Over Time",
       subtitle = "Broad sample (ages 30\u201349). Absorbing (ever-parent) definition.",
       x = NULL, y = "Employment rate (%)",
       caption = "Source: HILDA Waves 12\u201324.") +
  theme_paper +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

ggsave( file.path(outdir, "fig5_employment_waves.pdf"), fig5, width = 13/2.54, height = 8/2.54, units = "in")
ggsave( file.path(outdir, "fig5_employment_waves.png"), fig5, width = 13/2.54, height = 8/2.54, dpi = 500, units = "in")


# #############################################################################
# FIGURE 6: ROBUSTNESS FOREST PLOT
# #############################################################################

cat("--- Figure 6: Robustness forest plot ---\n")

# Estimate FE β₂ across specifications for childless men + never-mothers
est_fe <- function(data, gvar, gval, spec_label) {
  d <- data %>%
    filter(.data[[gvar]] == gval, in_wage_sample == 1) %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c) %>%
    filter(complete.cases(ln_hourly_wage_real, educ_exp_c, experience_years,
                          age_sq, married_num))
  if (nrow(d) < 80) return(NULL)
  m <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL)
  if (is.null(m) || !("educ_exp_c" %in% names(coef(m)))) return(NULL)
  s <- summary(m)$coeftable
  tibble(spec = spec_label, group = gval, N = nrow(d),
         beta2 = s["educ_exp_c", 1], se = s["educ_exp_c", 2])
}

specs <- list(
  list(panel, "Baseline"),
  list(panel %>% filter(age >= 25 & age <= 54), "Age 25\u201354"),
  list(panel %>% filter(age >= 35 & age <= 49), "Age 35\u201349"),
  list(panel %>% filter(fulltime == 1), "Full-time only"),
  list(panel %>% filter(fulltime == 0, employed == 1), "Part-time only"),
  list(panel %>% filter(!(wave %in% c(20, 21))), "Excl. COVID")
)

# Tighter trim
panel_trim <- panel %>%
  filter(in_wage_sample == 1) %>%
  group_by(female) %>%
  mutate(p01 = quantile(hourly_wage_clean, 0.01, na.rm = TRUE),
         p99 = quantile(hourly_wage_clean, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(hourly_wage_clean >= p01, hourly_wage_clean <= p99)
specs <- c(specs, list(list(panel_trim, "Trim 1/99")))

# No trim
panel_notrim <- panel %>%
  filter(employed == 1, !is.na(hourly_wage_clean), hourly_wage_clean > 0) %>%
  dplyr::select(-dplyr::any_of("deflator")) %>%   # drop carried deflator; avoids deflator.x/.y on join
  left_join(wpi_data, by = "wave") %>%
  mutate(ln_hourly_wage_real = log(hourly_wage_clean * deflator), in_wage_sample = 1L)
specs <- c(specs, list(list(panel_notrim, "No trim")))

rob <- tibble()
for (sp in specs) {
  for (g in c("childless_men_ever", "never_mothers")) {
    r <- est_fe(sp[[1]], "group_ever", g, sp[[2]])
    if (!is.null(r)) rob <- bind_rows(rob, r)
  }
}

rob$group_label <- factor(
  ifelse(rob$group == "childless_men_ever", "Childless men", "Never-mothers"),
  levels = c("Childless men", "Never-mothers"))
rob$spec <- factor(rob$spec, levels = rev(unique(rob$spec)))

cat("\n  Robustness FE b2 across specifications:\n")
cat(sprintf("  %-22s  %-20s  %7s  %10s\n",
            "Spec", "Group", "N", "FE b2 x100"))
cat("  ", strrep("-", 65), "\n")
for (i in 1:nrow(rob)) {
  r <- rob[i, ]
  p_val <- 2 * pnorm(-abs(r$beta2 / r$se))
  cat(sprintf("  %-22s  %-20s  %7d  %+.4f%%%s\n",
              as.character(r$spec), as.character(r$group_label),
              r$N, r$beta2 * 100, stars_fn(p_val)))
}

fig6 <- ggplot(rob, aes(x = beta2 * 100, y = spec, colour = group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = (beta2 - 1.96*se)*100,
                      xmax = (beta2 + 1.96*se)*100),
                  position = position_dodge(width = 0.6),
                  size = 0.5, linewidth = 0.7) +
  scale_colour_manual(values = c("Childless men" = "#2166AC",
                                 "Never-mothers" = "#B2182B")) +
  labs(
    title = expression(paste("Robustness: FE ", beta[2],
                             " Across Specifications")),
    subtitle = NULL,
    x = expression(hat(beta)[2] %*% 100), y = NULL,
    caption = "Source: HILDA Waves 12\u201324."
  ) +
  theme_paper +
  theme(legend.position = "right",
        panel.grid.major.y = element_line(colour = "grey92"))

ggsave( file.path(outdir, "fig6_robustness_forest.pdf"), fig6, width = 13/2.54, height = 9/2.54, units = "in")
ggsave( file.path(outdir, "fig6_robustness_forest.png"), fig6, width = 13/2.54, height = 9/2.54, dpi = 500, units = "in")


# #############################################################################
# FIGURE 7: AGE AT FIRST BIRTH GRADIENT
# #############################################################################

cat("--- Figure 7: Age-at-first-birth gradient ---\n")

# BUG FIX: compute centering reference means from the FULL mothers sample
# before the sub-group loop, so all AFB groups share a common reference.
mothers <- wage %>%
  filter(group_ever == "mothers_ever", !is.na(approx_age_at_first_birth)) %>%
  mutate(afb_group = case_when(
    approx_age_at_first_birth < 25 ~ "Early (<25)",
    approx_age_at_first_birth < 32 ~ "Mid (25\u201331)",
    TRUE ~ "Late (32+)"
  ))
educ_mean_ref_fig <- mean(mothers$educ_years,      na.rm = TRUE)
exp_mean_ref_fig  <- mean(mothers$experience_years, na.rm = TRUE)

afb_data <- tibble()
for (afb in c("Early (<25)", "Mid (25\u201331)", "Late (32+)")) {
  d <- mothers %>%
    filter(afb_group == afb) %>%
    mutate(educ_c     = educ_years      - educ_mean_ref_fig,
           exp_c      = experience_years - exp_mean_ref_fig,
           educ_exp_c = educ_c * exp_c)
  
  m <- tryCatch(
    feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
            age_sq + married_num | person_id + wave,
          data = d, cluster = ~person_id, notes = FALSE),
    error = function(e) NULL)
  
  if (!is.null(m) && "educ_exp_c" %in% names(coef(m))) {
    s <- summary(m)$coeftable
    afb_data <- bind_rows(afb_data,
                          tibble(afb_group = afb, N = nrow(d), N_id = n_distinct(d$person_id),
                                 mean_afb = mean(d$approx_age_at_first_birth),
                                 mean_educ = mean(d$educ_years, na.rm = TRUE),
                                 beta2 = s["educ_exp_c", 1], se = s["educ_exp_c", 2]))
  }
}

# Add never-mothers benchmark
nm_d <- wage %>%
  filter(group_ever == "never_mothers") %>%
  mutate(educ_c     = educ_years      - mean(educ_years,      na.rm = TRUE),
         exp_c      = experience_years - mean(experience_years, na.rm = TRUE),
         educ_exp_c = educ_c * exp_c)
nm_m <- feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
                age_sq + married_num | person_id + wave,
              data = nm_d, cluster = ~person_id, notes = FALSE)
nm_b2 <- summary(nm_m)$coeftable["educ_exp_c", 1]
nm_se <- summary(nm_m)$coeftable["educ_exp_c", 2]

afb_data$afb_group <- factor(afb_data$afb_group,
                             levels = c("Early (<25)", "Mid (25\u201331)", "Late (32+)"))

cat("\n  AFB gradient FE b2 by age at first birth:\n")
cat(sprintf("  %-14s  %7s %7s  %8s  %7s  %10s\n",
            "AFB group", "N", "N_id", "Mean AFB", "Mean Ed", "FE b2 x100"))
cat("  ", strrep("-", 65), "\n")
for (i in 1:nrow(afb_data)) {
  r <- afb_data[i, ]
  p_val <- 2 * pnorm(-abs(r$beta2 / r$se))
  cat(sprintf("  %-14s  %7d %7d  %8.1f  %7.1f  %+.4f%%%s\n",
              as.character(r$afb_group), r$N, r$N_id,
              r$mean_afb, r$mean_educ,
              r$beta2 * 100, stars_fn(p_val)))
}
cat(sprintf("  %-14s  (reference)  FE b2 = %+.4f%%%s\n",
            "Never-mothers", nm_b2 * 100,
            stars_fn(2 * pnorm(-abs(nm_b2 / nm_se)))))

fig7 <- ggplot(afb_data, aes(x = afb_group, y = beta2 * 100)) +
  # Never-mothers benchmark band
  annotate("rect", xmin = 0.4, xmax = 3.6,
           ymin = (nm_b2 - 1.96*nm_se)*100, ymax = (nm_b2 + 1.96*nm_se)*100,
           fill = "#B2182B", alpha = 0.08) +
  geom_hline(yintercept = nm_b2 * 100, linetype = "dashed", colour = "#B2182B",
             linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(ymin = (beta2 - 1.96*se)*100,
                      ymax = (beta2 + 1.96*se)*100),
                  colour = "#D6604D", size = 0.8, linewidth = 1) +
  # N and Educ labels: to the LEFT of each point, at the point y-value
  geom_text(aes(y = beta2 * 100,
                label = sprintf("N=%s\nEduc=%.1f yr",
                                format(N_id, big.mark = ","), mean_educ)),
            hjust = 1.25, size = 1.9, colour = "grey40", lineheight = 0.85) +
  # Never-mothers label: inside the pink CI band at left, below dashed line
  annotate("text", x = 0.5, y = (nm_b2 - 0.4*(nm_b2 - (nm_b2 - 1.96*nm_se))) * 100,
           label = "Never-\nmothers", colour = "#B2182B",
           size = 1.9, hjust = 0, vjust = 0.5, fontface = "italic",
           lineheight = 0.85) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.22))) +
  labs(
    title = expression(paste("Motherhood Timing and Complementarity: FE ",
                             beta[2], " by Age at First Birth")),
    subtitle = NULL,
    x = "Age at first birth", y = expression(hat(beta)[2] %*% 100),
    caption = "Source: HILDA Waves 12\u201324, ages 30\u201349. Red dashed line = never-mothers benchmark (95% CI band)."
  ) +
  theme_paper

ggsave( file.path(outdir, "fig7_age_first_birth.pdf"), fig7, width = 13/2.54, height = 10/2.54, units = "in")
ggsave( file.path(outdir, "fig7_age_first_birth.png"), fig7, width = 13/2.54, height = 10/2.54, dpi = 500, units = "in")


# #############################################################################
# FIGURE 8: COMPOSITION DRIFT
# #############################################################################

cat("--- Figure 8: Composition drift ---\n")

comp <- wage %>%
  filter(!is.na(group_ever)) %>%
  group_by(year, group_ever) %>%
  summarise(
    N = n(),
    educ = mean(educ_years, na.rm = TRUE),
    exp  = mean(experience_years, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group_label = factor(label_group(group_ever),
                              levels = c("Childless men", "Fathers", "Never-mothers", "Mothers (ever)")))

comp_long <- comp %>%
  pivot_longer(cols = c(educ, exp), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(
    ifelse(variable == "educ", "Mean education (years)", "Mean experience (years)"),
    levels = c("Mean education (years)", "Mean experience (years)")))

fig8 <- ggplot(comp_long, aes(x = year, y = value, colour = group_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~variable, scales = "free_y") +
  scale_colour_manual(values = pal4) +
  scale_x_continuous(breaks = seq(2012, 2024, 2)) +
  labs(title = "Sample Composition Over Time",
       subtitle = "Wage sample, absorbing (ever-parent) definition.",
       x = NULL, y = NULL,
       caption = "Source: HILDA Waves 12\u201324, ages 30\u201349.") +
  theme_paper +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

ggsave( file.path(outdir, "fig8_composition_drift.pdf"), fig8, width = 14/2.54, height = 8/2.54, units = "in")
ggsave( file.path(outdir, "fig8_composition_drift.png"), fig8, width = 14/2.54, height = 8/2.54, dpi = 500, units = "in")


# COPY ALL TO OUTPUTS

cat("\n--- Figures written to output/figures/ ---\n")
figs <- list.files(pattern = "^fig[0-9].*\\.(pdf|png)$")
for (f in figs) {
  ok <- file.copy(f, file.path(outdir, f), overwrite = TRUE)
  if (ok) cat(sprintf("  \u2713 %s\n", f))
}

# Save intermediate data for reproducibility
saveRDS(wbw, "wave_by_wave_beta2.rds")
saveRDS(ev_b2, "event_study_beta2.rds")
saveRDS(rob, "robustness_forest_data.rds")
saveRDS(afb_data, "age_first_birth_beta2.rds")
write_csv(ev_b2, "event_study_beta2.csv")
write_csv(rob, "robustness_forest_data.csv")

cat("\n")
cat(strrep("=", 70), "\n")
cat("  ALL 8 FIGURES COMPLETE\n")
cat(strrep("=", 70), "\n")


# COUPLE_FIGURES.R
# Couple-level figures for journal submission
#
# OUTPUT (6 figures, PDF + PNG):
#   Fig C1: Pooled vs within-couple β₂ (gender gap vanishes)
#   Fig C2: Four-cell DML β₂ with/without partner controls
#   Fig C3: Hours reallocation by parenthood (couple-level)
#   Fig C4: Within-couple hours and wage dynamics around birth
#   Fig C5: Causal forest variable importance (own vs partner)
#   Fig C6: CATE-experience gradient by group
#
# Requires: hilda_panel_data_extended.rds, tidyverse, fixest, grf (for CF)


pal_4cell <- c(
  "Childless men"   = "#2166AC",
  "Fathers"         = "#92C5DE",
  "Childless women" = "#B2182B",
  "Mothers"         = "#F4A582"
)




# 0. LOAD + BUILD COUPLE PANEL

cat(strrep("=", 70), "\n")
cat("  COUPLE FIGURES\n")
cat(strrep("=", 70), "\n\n")


# WPI deflators available for no-trim panel rebuild

num_rows    <- format(nrow(panel), big.mark = ",")
num_persons <- format(n_distinct(panel$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))

# --- Build couple panel ---
cat("\n--- Building couple panel ---\n")

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
  select(person_id, partner_id = partner_id_clean, wave, year, female, ever_parent,
         group_ever, ln_hourly_wage_real, educ_years, experience_years,
         hours_worked_clean, employed, fulltime, in_wage_sample, married, age_sq,
         had_birth_this_wave)

p2 <- panel %>%
  select(person_id, wave, female, ever_parent, ln_hourly_wage_real,
         educ_years, experience_years, hours_worked_clean, employed,
         fulltime, in_wage_sample, married, age_sq)

couples_raw <- p1 %>%
  inner_join(p2, by = c("partner_id" = "person_id", "wave" = "wave"),
             suffix = c("_f", "_m")) %>%
  filter(female_f != female_m)  # heterosexual only

# Standardise: her/his
couples <- couples_raw %>%
  mutate(
    her_wage  = ifelse(female_f == 1, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
    his_wage  = ifelse(female_f == 0, ln_hourly_wage_real_f, ln_hourly_wage_real_m),
    her_hours = ifelse(female_f == 1, hours_worked_clean_f, hours_worked_clean_m),
    his_hours = ifelse(female_f == 0, hours_worked_clean_f, hours_worked_clean_m),
    her_educ  = ifelse(female_f == 1, educ_years_f, educ_years_m),
    his_educ  = ifelse(female_f == 0, educ_years_f, educ_years_m),
    her_exp   = ifelse(female_f == 1, experience_years_f, experience_years_m),
    his_exp   = ifelse(female_f == 0, experience_years_f, experience_years_m),
    her_employed = ifelse(female_f == 1, employed_f, employed_m),
    his_employed = ifelse(female_f == 0, employed_f, employed_m),
    her_fulltime = ifelse(female_f == 1, fulltime_f, fulltime_m),
    his_fulltime = ifelse(female_f == 0, fulltime_f, fulltime_m),
    her_in_wage = ifelse(female_f == 1, in_wage_sample_f, in_wage_sample_m),
    his_in_wage = ifelse(female_f == 0, in_wage_sample_f, in_wage_sample_m),
    has_children = as.integer(ever_parent_f == 1 | ever_parent_m == 1),
    couple_id = ifelse(person_id < partner_id,
                       paste0(person_id, "_", partner_id),
                       paste0(partner_id, "_", person_id))
  )

n_couples <- n_distinct(couples$couple_id)
cat(sprintf("  Couples: %s, couple-waves: %s\n",
            format(n_couples, big.mark = ","),
            format(nrow(couples), big.mark = ",")))

# Dual-earner stacked for β₂
dual <- couples %>%
  filter(his_in_wage == 1, her_in_wage == 1,
         !is.na(his_wage), !is.na(her_wage))

his_stack <- dual %>%
  transmute(couple_id, wave, year, has_children,
            ln_wage = his_wage, educ = his_educ, exp = his_exp,
            hours = his_hours, is_female = 0L,
            partner_educ = her_educ, partner_exp = her_exp,
            partner_hours = her_hours)

her_stack <- dual %>%
  transmute(couple_id, wave, year, has_children,
            ln_wage = her_wage, educ = her_educ, exp = her_exp,
            hours = her_hours, is_female = 1L,
            partner_educ = his_educ, partner_exp = his_exp,
            partner_hours = his_hours)

stacked <- bind_rows(his_stack, her_stack) %>%
  mutate(educ_c = educ - mean(educ, na.rm = TRUE),
         exp_c  = exp - mean(exp, na.rm = TRUE),
         educ_exp_c = educ_c * exp_c,
         person_couple = paste0(couple_id, "_", is_female))

num_rows    <- format(nrow(stacked), big.mark = ",")
num_persons <- format(n_distinct(stacked$person_id), big.mark = ",")
cat(paste(num_rows, "obs |", num_persons, "individuals\n"))


# #############################################################################
# FIG C1: POOLED vs WITHIN-COUPLE β₂ (GENDER GAP VANISHES)
# #############################################################################

cat("\n--- Fig C1: Pooled vs within-couple β₂ ---\n")

# Estimate β₂ four ways: pooled OLS, pooled FE, couple OLS, couple FE
c1_data <- tibble()

for (sex in 0:1) {
  sex_lab <- ifelse(sex == 0, "Men", "Women")
  
  # A) Individual-level (from main wage sample, NOT couple-restricted)
  d_ind <- wage %>%
    filter(female == sex) %>%
    mutate(educ_c = educ_years - mean(educ_years, na.rm = TRUE),
           exp_c  = experience_years - mean(experience_years, na.rm = TRUE),
           educ_exp_c = educ_c * exp_c)
  
  m_ols <- feols(ln_hourly_wage_real ~ educ_years + educ_exp_c + experience_years +
                   age_sq + married | wave,
                 data = d_ind, cluster = ~person_id, notes = FALSE)
  m_fe  <- feols(ln_hourly_wage_real ~ educ_exp_c + experience_years +
                   age_sq + married | person_id + wave,
                 data = d_ind, cluster = ~person_id, notes = FALSE)
  
  s_ols <- summary(m_ols)$coeftable
  s_fe  <- summary(m_fe)$coeftable
  
  c1_data <- bind_rows(c1_data,
                       tibble(gender = sex_lab, spec = "Individual OLS",
                              beta2 = s_ols["educ_exp_c", 1], se = s_ols["educ_exp_c", 2]),
                       tibble(gender = sex_lab, spec = "Individual FE",
                              beta2 = s_fe["educ_exp_c", 1], se = s_fe["educ_exp_c", 2])
  )
  
  # B) Couple-stacked (with partner controls)
  d_couple <- stacked %>% filter(is_female == sex)
  
  m_c_ols <- feols(ln_wage ~ educ + educ_exp_c + exp + hours +
                     partner_educ + partner_exp + partner_hours | wave,
                   data = d_couple, cluster = ~couple_id, notes = FALSE)
  m_c_fe  <- feols(ln_wage ~ educ + educ_exp_c + exp + hours +
                     partner_educ + partner_exp + partner_hours | couple_id + wave,
                   data = d_couple, cluster = ~couple_id, notes = FALSE)
  
  s_c_ols <- summary(m_c_ols)$coeftable
  s_c_fe  <- summary(m_c_fe)$coeftable
  
  c1_data <- bind_rows(c1_data,
                       tibble(gender = sex_lab, spec = "Couple OLS\n(+ partner controls)",
                              beta2 = s_c_ols["educ_exp_c", 1], se = s_c_ols["educ_exp_c", 2]),
                       tibble(gender = sex_lab, spec = "Couple FE\n(+ partner controls)",
                              beta2 = s_c_fe["educ_exp_c", 1], se = s_c_fe["educ_exp_c", 2])
  )
}

c1_data$spec <- factor(c1_data$spec,
                       levels = c("Individual OLS", "Individual FE",
                                  "Couple OLS\n(+ partner controls)", "Couple FE\n(+ partner controls)"))

# Print for verification
cat("\n")
for (i in 1:nrow(c1_data)) {
  r <- c1_data[i,]
  cat(sprintf("  %-6s  %-30s  \u03b2\u2082=%+.4f%%  SE=%.4f\n",
              r$gender, gsub("\n", " ", r$spec), r$beta2*100, r$se*100))
}

# Compute gender gap for annotation
gap_ind_ols <- diff(c1_data$beta2[c1_data$spec == "Individual OLS"]) * 100
gap_ind_fe  <- diff(c1_data$beta2[c1_data$spec == "Individual FE"]) * 100
gap_cpl_fe  <- diff(c1_data$beta2[c1_data$spec == "Couple FE\n(+ partner controls)"]) * 100

fig_c1 <- ggplot(c1_data,
                 aes(x = spec, y = beta2 * 100, colour = gender)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = (beta2 - 1.96*se)*100,
                      ymax = (beta2 + 1.96*se)*100),
                  position = position_dodge(width = 0.5),
                  size = 0.7, linewidth = 0.9) +
  # Gap annotations
  annotate("segment", x = 0.75, xend = 1.25,
           y = max(c1_data$beta2[c1_data$spec == "Individual OLS"])*100,
           yend = max(c1_data$beta2[c1_data$spec == "Individual OLS"])*100,
           colour = "grey40", linewidth = 0.3) +
  annotate("text", x = 1, vjust = 0,
           y = (max(c1_data$beta2[c1_data$spec == "Individual OLS"]) +
                  1.96 * max(c1_data$se[c1_data$spec == "Individual OLS"])) * 100 + 0.008,
           label = sprintf("Gap: %.3f%%", abs(gap_ind_ols)),
           size = 2.2, colour = "grey30", hjust = 0.5) +
  annotate("text", x = 4, y = max(c1_data$beta2*100) + 0.02,
           label = sprintf("Gap: %.3f%%", abs(gap_cpl_fe)),
           size = 2.2, colour = "grey30") +
  scale_colour_manual(values = pal_gender) +
  labs(
    title = expression(paste("The Gender Gap in ", beta[2],
                             " Narrows Within Couples")),
    subtitle = "Individual-level estimates show a gender gap; couple FE with partner controls absorbs it.",
    x = NULL, y = expression(hat(beta)[2] %*% 100),
    caption = "Source: HILDA Waves 12\u201324. Couple FE absorbs between-couple heterogeneity."
  ) +
  theme_paper +
  theme(axis.text.x = element_text(size = 6))

ggsave( file.path(outdir, "figC1_gender_gap_vanishes.pdf"), fig_c1, width = 13/2.54, height = 9/2.54, units = "in")
ggsave( file.path(outdir, "figC1_gender_gap_vanishes.png"), fig_c1, width = 13/2.54, height = 9/2.54, dpi = 500, units = "in")


# #############################################################################
# FIG C2: FOUR-CELL β₂ WITH/WITHOUT PARTNER CONTROLS
# #############################################################################

cat("\n--- Fig C2: Four-cell \u03b2\u2082 with/without partner controls ---\n")

c2_data <- tibble()
for (sex in 0:1) {
  for (par in 0:1) {
    lab <- case_when(
      sex == 0 & par == 0 ~ "Childless men",
      sex == 0 & par == 1 ~ "Fathers",
      sex == 1 & par == 0 ~ "Childless women",
      sex == 1 & par == 1 ~ "Mothers"
    )
    d <- stacked %>% filter(is_female == sex, has_children == par)
    if (nrow(d) < 100) next
    
    # Without partner controls
    m_own <- tryCatch(
      feols(ln_wage ~ educ + educ_exp_c + exp + hours | wave,
            data = d, cluster = ~couple_id, notes = FALSE),
      error = function(e) NULL)
    
    # With partner controls
    m_part <- tryCatch(
      feols(ln_wage ~ educ + educ_exp_c + exp + hours +
              partner_educ + partner_exp + partner_hours | wave,
            data = d, cluster = ~couple_id, notes = FALSE),
      error = function(e) NULL)
    
    for (m_info in list(list(m_own, "Own controls"), list(m_part, "+ Partner controls"))) {
      m <- m_info[[1]]; est_lab <- m_info[[2]]
      if (!is.null(m) && "educ_exp_c" %in% names(coef(m))) {
        s <- summary(m)$coeftable
        c2_data <- bind_rows(c2_data,
                             tibble(group = lab, controls = est_lab, N = nrow(d),
                                    beta2 = s["educ_exp_c", 1], se = s["educ_exp_c", 2]))
      }
    }
  }
}

c2_data$group <- factor(c2_data$group,
                        levels = c("Childless men", "Fathers", "Childless women", "Mothers"))

fig_c2 <- ggplot(c2_data,
                 aes(x = group, y = beta2 * 100, colour = controls)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = (beta2 - 1.96*se)*100,
                      ymax = (beta2 + 1.96*se)*100),
                  position = position_dodge(width = 0.5),
                  size = 0.6, linewidth = 0.8) +
  scale_colour_manual(values = c("Own controls" = "#4393C3",
                                 "+ Partner controls" = "#D6604D")) +
  labs(
    title = bquote("Partner Controls and " * beta[2] * " by Gender \u00d7 Parenthood"),
    subtitle = "Adding partner education, experience, and hours. Couple-stacked OLS with wave FE.",
    x = NULL, y = expression(hat(beta)[2] %*% 100),
    colour = NULL,
    caption = "Source: HILDA Waves 12\u201324. Dual-earner couples only."
  ) +
  theme_paper +
  theme(legend.position = "right")

ggsave( file.path(outdir, "figC2_four_cell_partner.pdf"), fig_c2, width = 13/2.54, height = 8/2.54, units = "in")
ggsave( file.path(outdir, "figC2_four_cell_partner.png"), fig_c2, width = 13/2.54, height = 8/2.54, dpi = 500, units = "in")


# #############################################################################
# FIG C3: HOURS REALLOCATION BY PARENTHOOD
# #############################################################################

cat("\n--- Fig C3: Hours reallocation ---\n")

hours_data <- couples %>%
  filter(!is.na(his_hours), !is.na(her_hours)) %>%
  group_by(has_children) %>%
  summarise(
    his_hours = mean(his_hours, na.rm = TRUE),
    her_hours = mean(her_hours, na.rm = TRUE),
    his_ft_pct = mean(his_fulltime == 1, na.rm = TRUE) * 100,
    her_ft_pct = mean(her_fulltime == 1, na.rm = TRUE) * 100,
    his_emp = mean(his_employed == 1, na.rm = TRUE) * 100,
    her_emp = mean(her_employed == 1, na.rm = TRUE) * 100,
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(status = ifelse(has_children == 0, "Childless couples", "Parent couples"))

cat("\n  Within-couple hours by parenthood status:\n")
cat(sprintf("  %-18s  %7s  %7s  %7s  %7s  %7s\n",
            "Status", "N", "HisHrs", "HerHrs", "HisFT%", "HerFT%"))
cat("  ", strrep("-", 62), "\n")
for (i in 1:nrow(hours_data)) {
  r <- hours_data[i, ]
  cat(sprintf("  %-18s  %7d  %7.1f  %7.1f  %6.1f%%  %6.1f%%\n",
              r$status, r$N, r$his_hours, r$her_hours,
              r$his_ft_pct, r$her_ft_pct))
}

# Reshape for plotting
h_long <- hours_data %>%
  select(status, his_hours, her_hours, his_ft_pct, her_ft_pct) %>%
  pivot_longer(-status, names_to = "var", values_to = "val") %>%
  mutate(
    gender = ifelse(grepl("^his", var), "Men", "Women"),
    measure = ifelse(grepl("hours", var), "Weekly hours", "Full-time rate (%)")
  )

p3a <- h_long %>% filter(measure == "Weekly hours") %>%
  ggplot(aes(x = status, y = val, fill = gender)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", val)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 2.2) +
  scale_fill_manual(values = pal_gender) +
  labs(title = "A. Mean weekly hours", y = "Hours", x = NULL) +
  coord_cartesian(ylim = c(0, 50)) +
  theme_paper + theme(legend.position = "none")

p3b <- h_long %>% filter(measure == "Full-time rate (%)") %>%
  ggplot(aes(x = status, y = val, fill = gender)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", val)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 2.2) +
  scale_fill_manual(values = pal_gender) +
  labs(title = "B. Full-time rate", y = "%", x = NULL) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_paper

fig_c3 <- (p3a | p3b) +
  plot_annotation(
    title = "Within-Couple Hours Reallocation by Parenthood Status",
    subtitle = "Men's hours barely change; women's hours and FT rates drop sharply in parent couples.",
    caption = "Source: HILDA Waves 12\u201324. Matched heterosexual couples.",
    theme = theme(plot.title = element_text(size = 7.5, face = "bold"),
                  plot.subtitle = element_text(size = 6, colour = "grey40"),
                  plot.caption = element_text(size = 6, colour = "grey50")))

ggsave( file.path(outdir, "figC3_hours_reallocation.pdf"), fig_c3, width = 13/2.54, height = 8/2.54, units = "in")
ggsave( file.path(outdir, "figC3_hours_reallocation.png"), fig_c3, width = 13/2.54, height = 8/2.54, dpi = 500, units = "in")


# #############################################################################
# FIG C4: WITHIN-COUPLE DYNAMICS AROUND BIRTH
# #############################################################################

cat("\n--- Fig C4: Within-couple dynamics around birth ---\n")

# Find couples with birth events
couple_births <- couples %>%
  filter(had_birth_this_wave == 1) %>%
  group_by(couple_id) %>%
  summarise(birth_wave = min(wave), .groups = "drop")

cat(sprintf("  Couples with birth during panel: %s\n",
            format(nrow(couple_births), big.mark = ",")))

couple_event <- couples %>%
  inner_join(couple_births, by = "couple_id") %>%
  mutate(t = wave - birth_wave) %>%
  filter(t >= -3, t <= 4)

ce_summary <- couple_event %>%
  group_by(t) %>%
  summarise(
    N = n(),
    his_hours = mean(his_hours, na.rm = TRUE),
    her_hours = mean(her_hours, na.rm = TRUE),
    hours_gap = mean(his_hours - her_hours, na.rm = TRUE),
    his_ft = mean(his_fulltime == 1, na.rm = TRUE) * 100,
    her_ft = mean(her_fulltime == 1, na.rm = TRUE) * 100,
    her_emp = mean(her_employed == 1, na.rm = TRUE) * 100,
    his_emp = mean(his_employed == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n  Within-couple dynamics around birth:\n")
cat(sprintf("  %3s  %5s  %7s  %7s  %8s  %7s  %7s  %7s  %7s\n",
            "t", "N", "HisHrs", "HerHrs", "HrsGap",
            "HisFT%", "HerFT%", "HisEmp", "HerEmp"))
cat("  ", strrep("-", 78), "\n")

for (i in 1:nrow(ce_summary)) {
  r <- ce_summary[i, ]
  cat(sprintf(
    "  %+2d   %5d  %7.1f  %7.1f  %+8.1f  %6.1f%%  %6.1f%%  %6.1f%%  %6.1f%%\n",
    r$t, r$N, r$his_hours, r$her_hours, r$hours_gap,
    r$his_ft, r$her_ft, r$his_emp, r$her_emp
  ))
}
# Panel A: hours
ce_hours <- ce_summary %>%
  select(t, his_hours, her_hours) %>%
  pivot_longer(-t, names_to = "who", values_to = "hours") %>%
  mutate(gender = ifelse(who == "his_hours", "Men", "Women"))

p4a <- ggplot(ce_hours, aes(x = t, y = hours, colour = gender)) +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = pal_gender) +
  scale_x_continuous(breaks = -3:4) +
  labs(title = "A. Weekly hours", y = "Hours", x = NULL) +
  theme_paper + theme(legend.position = "none")

# Panel B: hours gap
p4b <- ggplot(ce_summary, aes(x = t, y = hours_gap)) +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
  geom_line(colour = "#7570B3", linewidth = 1) +
  geom_point(colour = "#7570B3", size = 2.5) +
  scale_x_continuous(breaks = -3:4) +
  labs(title = "B. Hours gap (his \u2212 hers)", y = "Hours", x = NULL) +
  theme_paper

# Panel C: her employment + FT
ce_emp <- ce_summary %>%
  select(t, her_emp, her_ft) %>%
  pivot_longer(-t, names_to = "measure", values_to = "pct") %>%
  mutate(measure = ifelse(measure == "her_emp", "Employment", "Full-time"))

p4c <- ggplot(ce_emp, aes(x = t, y = pct, linetype = measure)) +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
  geom_line(colour = "#B2182B", linewidth = 1) +
  geom_point(colour = "#B2182B", size = 2) +
  scale_x_continuous(breaks = -3:4) +
  labs(title = "C. Emp. & FT rate",
       y = "%", x = "Years from birth") +
  theme_paper + theme(legend.position = "bottom")

fig_c4 <- (p4a | p4b | p4c) +
  plot_annotation(
    title = "Within-Couple Labour Market Dynamics Around First Birth",
    subtitle = "Matched couples observed before and after birth. Vertical line = birth event.",
    caption = "Source: HILDA Waves 12\u201324. Matched heterosexual couples with observed birth during panel.",
    theme = theme(plot.title = element_text(size = 7.5, face = "bold"),
                  plot.subtitle = element_text(size = 6, colour = "grey40"),
                  plot.caption = element_text(size = 6, colour = "grey50")))

ggsave( file.path(outdir, "figC4_couple_birth_dynamics.pdf"), fig_c4, width = 18/2.54, height = 10/2.54, units = "in")
ggsave( file.path(outdir, "figC4_couple_birth_dynamics.png"), fig_c4, width = 18/2.54, height = 10/2.54, dpi = 500, units = "in")


###########################################################################
# FIG C5: TRIPLE INTERACTION — IMPLIED β₂ BY CELL
###########################################################################

cat("\n--- Fig C5: Triple interaction implied β₂ ---\n")

###########################################################################
# 1. OLS (STACKED)
###########################################################################

triple_ols <- feols(
  ln_wage ~ educ_exp_c * is_female * has_children + educ + exp + hours | wave,
  data = stacked,
  cluster = ~couple_id,
  notes = FALSE
)

###########################################################################
# 2. DELTA-METHOD IMPLIED β₂
###########################################################################

get_implied_delta <- function(mod, label) {
  
  b  <- coef(mod)
  V  <- vcov(mod)
  
  # helper to safely extract coefficient index
  idx <- function(nm) match(nm, names(b), nomatch = 0)
  
  coef_vec <- function(weights) {
    w <- rep(0, length(b))
    names(w) <- names(b)
    w[names(weights)] <- weights
    w
  }
  
  cells <- list(
    
    "Childless men" =
      coef_vec(c("educ_exp_c" = 1)),
    
    "Childless women" =
      coef_vec(c("educ_exp_c" = 1,
                 "educ_exp_c:is_female" = 1)),
    
    "Fathers" =
      coef_vec(c("educ_exp_c" = 1,
                 "educ_exp_c:has_children" = 1)),
    
    "Mothers" =
      coef_vec(c("educ_exp_c" = 1,
                 "educ_exp_c:is_female" = 1,
                 "educ_exp_c:has_children" = 1,
                 "educ_exp_c:is_female:has_children" = 1))
  )
  
  map_dfr(names(cells), function(g) {
    
    a <- cells[[g]]
    
    beta_hat <- sum(a * b)
    
    # delta-method SE
    se_hat <- sqrt( t(a) %*% V %*% a )
    
    tibble(
      group = g,
      beta2 = beta_hat,
      se = as.numeric(se_hat),
      estimator = label
    )
  })
}

ols_cells <- get_implied_delta(triple_ols, "OLS (stacked)")

cat("\n  Implied b2 by cell (OLS stacked, delta-method SEs):\n")
cat(sprintf("  %-18s  %10s  %8s\n", "Group", "b2 x100", "SE x100"))
cat("  ", strrep("-", 42), "\n")
for (i in 1:nrow(ols_cells)) {
  r <- ols_cells[i, ]
  p_val <- 2 * pnorm(-abs(r$beta2 / r$se))
  cat(sprintf("  %-18s  %+.4f%%%s  %.4f%%\n",
              r$group, r$beta2 * 100, stars_fn(p_val), r$se * 100))
}

###########################################################################
# 3. DML (RF) — replace Couple FE
###########################################################################

library(dplyr)
library(purrr)
library(tibble)

# Supplementary figC5 (OLS vs DML-RF implied beta2 by cell). dml_cells_rf is
# built in MASTER_hh.R Part D and saved only in master_state_after_partD.RData,
# so guard the whole block and skip cleanly when it is not in the workspace.
if (exists("dml_cells_rf") && length(dml_cells_rf) > 0) {
  
  # Debug: check what dml_cells_rf contains
  cat(sprintf("\n  dml_cells_rf has %d entries: %s\n",
              length(dml_cells_rf), paste(names(dml_cells_rf), collapse = ", ")))
  if (length(dml_cells_rf) > 0) {
    cat(sprintf("  First entry structure: theta=%s, se=%s\n",
                paste(class(dml_cells_rf[[1]]$theta), collapse="/"),
                paste(class(dml_cells_rf[[1]]$se), collapse="/")))
  }
  
  # Build DML tibble — force scalars with as.numeric()[1]
  dml_cells <- imap_dfr(dml_cells_rf, function(x, group_name) {
    if (is.null(x)) return(tibble())
    tibble(
      group     = group_name,
      beta2     = as.numeric(x$theta)[1],
      se        = as.numeric(x$se)[1],
      estimator = "DML (RF)"
    )
  })
  
  cat(sprintf("\n  DML cells: %d rows, columns: %s\n",
              nrow(dml_cells), paste(names(dml_cells), collapse = ", ")))
  
  cat("\n  DML (RF) implied β₂ by cell:\n")
  cat(sprintf("  %-18s  %10s  %8s\n", "Group", "β₂ x100", "SE x100"))
  cat("  ", strrep("-", 42), "\n")
  for (i in seq_len(nrow(dml_cells))) {
    r <- dml_cells[i, ]
    p_val <- 2 * pnorm(-abs(r$beta2 / r$se))
    cat(sprintf("  %-18s  %+.4f%%%s  %.4f%%\n",
                r$group, r$beta2 * 100, stars_fn(p_val), r$se * 100))
  }
  
  ###########################################################################
  # 4. COMBINE DATA
  ###########################################################################
  
  # Ensure both have exactly the same columns for bind_rows
  ols_for_plot <- ols_cells %>% select(group, beta2, se, estimator)
  dml_for_plot <- dml_cells %>% select(group, beta2, se, estimator)
  
  cat(sprintf("\n  OLS rows: %d, DML rows: %d\n", nrow(ols_for_plot), nrow(dml_for_plot)))
  cat(sprintf("  OLS beta2 range: [%.4f, %.4f]\n", min(ols_for_plot$beta2)*100, max(ols_for_plot$beta2)*100))
  if (nrow(dml_for_plot) > 0) {
    cat(sprintf("  DML beta2 range: [%.4f, %.4f]\n", min(dml_for_plot$beta2)*100, max(dml_for_plot$beta2)*100))
  }
  
  c5_data <- bind_rows(ols_for_plot, dml_for_plot)
  
  cat(sprintf("  Combined: %d rows, any NA in beta2: %s\n",
              nrow(c5_data), any(is.na(c5_data$beta2))))
  
  c5_data$group <- factor(
    c5_data$group,
    levels = c("Childless men","Childless women","Fathers","Mothers")
  )
  
  ###########################################################################
  # 5. DYNAMIC ANNOTATION VALUES
  ###########################################################################
  
  # OLS: gender gap p-value from triple interaction model
  extract_p <- function(model, coef_name) {
    # Get the summary object
    s <- summary(model, infer = TRUE)
    
    # Determine where the coefficient table is
    tab <- if (!is.null(s$coeftable)) {
      s$coeftable
    } else if (!is.null(s$coefficients) && is.matrix(s$coefficients)) {
      s$coefficients
    } else if (!is.null(s$coefficients) && is.numeric(s$coefficients)) {
      # If it's a named numeric vector (like your felm summary sometimes)
      stop("Cannot extract p-values: coefficient table not found; only point estimates exist.")
    } else {
      stop("No coefficient table found.")
    }
    
    # If column names exist, find p-value column
    if (!is.null(colnames(tab))) {
      p_col <- grep("^Pr|P\\|>t\\|", colnames(tab))[1]
    } else {
      # fallback: p-value is usually the last column
      p_col <- ncol(tab)
    }
    
    # Return the p-value for the requested coefficient
    if (!(coef_name %in% rownames(tab))) {
      stop(paste0("Coefficient '", coef_name, "' not found in model."))
    }
    
    tab[coef_name, p_col]
  }
  
  p_gender_ols <- extract_p(triple_ols, "educ_exp_c:is_female")
  
  # DML: triple difference p-value (fatherhood penalty vs motherhood penalty)
  # Fatherhood penalty = Fathers θ − Childless men θ
  # Motherhood penalty = Mothers θ − Childless women θ
  # Triple = fatherhood penalty − motherhood penalty
  if (all(c("Childless men", "Fathers", "Childless women", "Mothers") %in%
          names(dml_cells_rf))) {
    fp_theta <- dml_cells_rf[["Fathers"]]$theta - dml_cells_rf[["Childless men"]]$theta
    fp_se    <- sqrt(dml_cells_rf[["Fathers"]]$se^2 + dml_cells_rf[["Childless men"]]$se^2)
    mp_theta <- dml_cells_rf[["Mothers"]]$theta - dml_cells_rf[["Childless women"]]$theta
    mp_se    <- sqrt(dml_cells_rf[["Mothers"]]$se^2 + dml_cells_rf[["Childless women"]]$se^2)
    triple_diff <- fp_theta - mp_theta
    triple_se   <- sqrt(fp_se^2 + mp_se^2)
    z_triple    <- triple_diff / triple_se
    p_triple_dml <- 2 * pnorm(-abs(z_triple))
  } else {
    p_triple_dml <- NA_real_
  }
  
  cat(sprintf("\n  Annotation p-values: OLS gender gap = %.4f, DML triple = %.4f\n",
              p_gender_ols, p_triple_dml))
  
  # Print full combined data for verification
  cat("\n  Combined c5_data for plot:\n")
  print(c5_data)
  
  
  ###########################################################################
  # 6. PLOT
  ###########################################################################
  
  # Compute annotation y-position safely
  annot_y <- max(c5_data$beta2 * 100, na.rm = TRUE) + 0.05
  
  fig_c5 <- ggplot(
    c5_data,
    aes(x = group, y = beta2 * 100, colour = estimator)
  ) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               colour = "grey60") +
    
    geom_pointrange(
      aes(ymin = (beta2 - 1.96*se)*100,
          ymax = (beta2 + 1.96*se)*100),
      position = position_dodge(width = 0.5),
      size = 0.6,
      linewidth = 0.8
    ) +
    
    scale_colour_manual(
      values = c(
        "OLS (stacked)" = "#2166AC",
        "DML (RF)"      = "#B2182B"
      )
    ) +
    
    # Gender gap annotation (OLS)
    annotate(
      "text",
      x = 1.5,
      y = annot_y,
      label = sprintf("Gender gap within childless couples\n(OLS): p=%.3f", p_gender_ols),
      size = 2.0,
      colour = "#2166AC",
      hjust = 0.5,
      lineheight = 0.9
    ) +
    
    # Triple penalty annotation (DML)
    annotate(
      "text",
      x = 3.5,
      y = annot_y,
      label = sprintf("Symmetric parenthood penalty across genders\n(DML): p=%.3f", p_triple_dml),
      size = 2.0,
      colour = "#B2182B",
      hjust = 0.5,
      lineheight = 0.9
    ) +
    
    labs(
      title = expression(
        paste("Triple Interaction: Implied ", beta[2],
              " by Gender × Parenthood (Within Couples)")
      ),
      subtitle = expression(
        paste("From: ", ln(w) == beta[2] %.%
                "(educ×exp) × female × parent + controls")
      ),
      x = NULL,
      y = expression("Implied " * hat(beta)[2] %*% 100),
      colour = NULL,
      caption = "Source: HILDA Waves 12–24. Couple-stacked dual-earner panel."
    ) +
    
    theme_paper +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 5.5)
    )
  
  ###########################################################################
  # 7. SAVE
  ###########################################################################
  
  ggsave(file.path(outdir, "figC5_triple_interaction.pdf"),
         fig_c5,
         width = 13/2.54,
         height = 9/2.54,
         units = "in")
  
  ggsave(file.path(outdir, "figC5_triple_interaction.png"),
         fig_c5,
         width = 13/2.54,
         height = 9/2.54,
         dpi = 500,
         units = "in")
  
} else {
  cat("\n[WARN] dml_cells_rf not in workspace; skipping supplementary figC5_triple_interaction.\n",
      "       It is built in MASTER_hh.R Part D and saved in master_state_after_partD.RData;\n",
      "       run the master through Part D (and re-save that checkpoint) to include it.\n", sep = "")
}



# #############################################################################
# FIG C6: COUPLE-LEVEL WAGE TRAJECTORIES (educ × exp profile)
# #############################################################################

cat("\n--- Fig C6: Predicted wage profiles ---\n")

# =============================================================================
# FIG 15 (C6): SCATTERPLOT — Education × Experience Wage Profiles
#
# Replaces the predicted-line version with a scatterplot showing raw data
# + linear fits for high vs low education, on a uniform y-axis.
#
# PREREQUISITES: Run after the couple-stacked panel ('stacked') is built
# in HILDA_within_couples_MASTER.R (around line 2400). This code block
# can replace the existing FIG C6 block (lines 4991–5090).
# =============================================================================

cat("\n--- Fig C6 (scatterplot): Education × experience wage profiles ---\n")

# ── 1. Classify workers as high/low education (median split within group) ────

scatter_data <- stacked %>%
  filter(!is.na(ln_wage), !is.na(exp), !is.na(educ)) %>%
  mutate(
    group = case_when(
      is_female == 0 & has_children == 0 ~ "Childless men",
      is_female == 0 & has_children == 1 ~ "Fathers",
      is_female == 1 & has_children == 0 ~ "Childless women",
      is_female == 1 & has_children == 1 ~ "Mothers"
    )
  ) %>%
  filter(!is.na(group))

# Within each group, define high/low education relative to group mean
scatter_data <- scatter_data %>%
  group_by(group) %>%
  mutate(
    educ_group_mean = mean(educ, na.rm = TRUE),
    educ_label = ifelse(
      educ >= educ_group_mean + 1.5,
      "High educ (mean + 3 yr)",
      ifelse(educ <= educ_group_mean - 1.5,
             "Low educ (mean \u2212 3 yr)",
             NA_character_)
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(educ_label))  # Keep only the tails, not the middle

cat(sprintf("  Scatterplot observations: %s (high + low education tails)\n",
            format(nrow(scatter_data), big.mark = ",")))
cat(sprintf("  By group:\n"))
scatter_data %>%
  count(group, educ_label) %>%
  mutate(label = sprintf("    %-20s %-25s N=%s", group, educ_label,
                         format(n, big.mark = ","))) %>%
  pull(label) %>%
  cat(sep = "\n")
cat("\n")

# ── 2. Residualise wages (remove wave FE) for cleaner visual ─────────────────

# Quick within-wave demeaning so profiles are comparable across waves
scatter_data <- scatter_data %>%
  group_by(group, wave) %>%
  mutate(ln_wage_demean = ln_wage - mean(ln_wage, na.rm = TRUE)) %>%
  ungroup()

# ── 3. Set factor levels ─────────────────────────────────────────────────────

scatter_data$group <- factor(scatter_data$group,
                             levels = c("Childless men", "Childless women",
                                        "Fathers", "Mothers"))

scatter_data$educ_label <- factor(scatter_data$educ_label,
                                  levels = c("Low educ (mean \u2212 3 yr)",
                                             "High educ (mean + 3 yr)"))

# ── 4. Plot ──────────────────────────────────────────────────────────────────

fig_c6 <- ggplot(scatter_data,
                 aes(x = exp, y = ln_wage_demean, colour = educ_label)) +
  # Raw data as transparent points
  geom_point(alpha = 0.04, size = 0.3, shape = 16) +
  # Linear fits with SE bands
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.8, alpha = 0.2) +
  # Uniform y-axis across all panels
  facet_wrap(~group, nrow = 1) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  # Colours
  scale_colour_manual(
    values = c("Low educ (mean \u2212 3 yr)"  = "#2166AC",
               "High educ (mean + 3 yr)" = "#B2182B"),
    labels = c("Low educ (mean \u2212 3 yr)"  = "Low educ (mean \u2212 3 yr)",
               "High educ (mean + 3 yr)" = "High educ (mean + 3 yr)")
  ) +
  labs(
    title    = NULL,
    x        = "Experience (years)",
    y        = "Demeaned ln(hourly wage)",
    colour   = NULL,
    caption  = "Source: HILDA Waves 12\u201324. Couple-stacked panel. Wave-demeaned within group."
  ) +
  theme_paper +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(0.7, "cm"),
    legend.key.size  = unit(0.3, "cm"),
    legend.text      = element_text(size = 6),
    strip.text       = element_text(size = 7, face = "bold"),
    plot.caption     = element_text(size = 5, hjust = 0)
  )

# ── 5. Save ──────────────────────────────────────────────────────────────────

ggsave(file.path(outdir, "fig15C6_wage_profiles.pdf"), fig_c6,
       width = 16/2.54, height = 8/2.54, units = "in")
ggsave(file.path(outdir, "fig15C6_wage_profiles.png"), fig_c6,
       width = 16/2.54, height = 8/2.54, dpi = 500, units = "in")

# ============================================================================
# Locally-built figures (not in the master): IV forest + within-person hours
# ============================================================================
cat("\n[C] IV coefficient plot + within-person hours adjustment ...\n")

## ---- fig_iv_coefplot ----
.iv_specs <- c("OLS", "IV pooled", "IV Reform 1 (preferred)", "IV falsification (<\\$80k)")
if (file.exists("ftbb_reform1_iv_results.rds")) {
  .g_iv <- readRDS("ftbb_reform1_iv_results.rds")$summary
  iv_df <- tibble(
    spec      = .iv_specs,
    theta     = .g_iv$theta_hat,
    se        = .g_iv$se,
    F_first   = .g_iv$fs_F_Sct,
    N         = .g_iv$n_obs,
    preferred = c(FALSE, FALSE, TRUE, FALSE)
  )
  cat("  fig_iv_coefplot: numbers read from ftbb_reform1_iv_results.rds\n")
} else {
  iv_df <- tibble(
    spec      = .iv_specs,
    theta     = c(-0.203,+0.914,-1.030,-0.087),
    se        = c( 0.048, 2.600, 1.010, 5.300),
    F_first   = c(   NA,  12.5,  20.1,  2.14),
    N         = c(18990,18990,  4650,   756),
    preferred = c(FALSE,FALSE, TRUE,  FALSE)
  )
  cat("  fig_iv_coefplot: ftbb_reform1_iv_results.rds not found; using fallback literals\n")
}
iv_df <- iv_df %>% mutate(
  ci_lo = theta - 1.96*se, ci_hi = theta + 1.96*se,
  label = ifelse(is.na(F_first), sprintf("N=%s", format(N, big.mark=",")),
                 sprintf("F=%.1f, N=%s", F_first, format(N, big.mark=","))),
  spec  = factor(spec, levels = rev(spec)))

p_iv <- ggplot(iv_df, aes(x = theta, y = spec, colour = preferred)) +
  geom_vline(xintercept = -0.2, linetype = "dashed", colour = "grey50") +
  annotate("rect", xmin = -0.4, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.3) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.20, linewidth = 0.6) +
  geom_point(size = 2.6) +
  geom_text(aes(x = ci_hi + 0.4, label = label), colour = "grey25", hjust = 0, size = 2.1) +
  scale_colour_manual(values = c(`TRUE`="#B2182B", `FALSE`="grey35"), guide = "none") +
  scale_x_continuous(limits = c(-13,16), breaks = seq(-12,16,4), expand = expansion(mult = c(0.02,0.02))) +
  labs(title = "FTB-B simulated instrument: first stage and imprecise contemporaneous wage response",
       subtitle = NULL,
       x = expression(hat(theta)^{IV}~~"(effect of "*S[ct]*" on log hourly wage)"), y = NULL,
       caption = "HILDA W15–W24. Reform 1: $100k–$150k baseline. 95% CI, cluster: couple. Grey band = OLS scale reference.") +
  theme_paper
if (FALSE) {  # fig_iv_coefplot removed from the paper; block disabled
ggplot2::ggsave(file.path(outdir, "fig_iv_coefplot.pdf"), p_iv, width = 16/2.54, height = 10/2.54, units = "in", device = .pdf_device)
cat("  figures/fig_iv_coefplot.pdf\n")
}

## ---- fig_hw_within_fe ----
## Within-person FE hours adjustment upon parenthood (master post-fix / paper Table 8).
## Hardcoded so the figure is self-contained and matches the table exactly.
hw_df <- tibble(
  sex     = c("Men","Men","Women","Women"),
  outcome = c("Housework hours","Market hours","Housework hours","Market hours"),
  delta   = c( 1.41,  0.04,  7.58, -11.21),
  se      = c( 0.13,  0.27,  0.21,   0.36)
) %>% mutate(ci_lo = delta - 1.96*se, ci_hi = delta + 1.96*se,
             sex = factor(sex, levels = c("Men","Women")),
             outcome = factor(outcome, levels = c("Housework hours","Market hours")))
p_hw <- ggplot(hw_df, aes(x = sex, y = delta, fill = outcome)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.55,
           colour = "grey20", linewidth = 0.35) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.5,
                position = position_dodge(width = 0.7)) +
  geom_text(aes(label = sprintf("%+.2f", delta),
                y = ifelse(delta >= 0, ci_hi + 0.6, ci_lo - 0.8)),
            position = position_dodge(width = 0.7), size = 2.4, colour = "grey15", fontface = "bold") +
  scale_fill_manual(values = c("Housework hours"="#B2182B","Market hours"="#2166AC")) +
  scale_y_continuous(limits = c(-15.5,10), breaks = seq(-15,10,5), expand = expansion(mult = c(0.03,0.05))) +
  labs(title = "Within-person hours adjustment upon parenthood",
       subtitle = "Becoming a parent: women add 7.6 hrs housework, lose 11.2 hrs market; men barely change",
       y = "\u0394 weekly hours (parent vs pre-parent)", x = NULL,
       caption = "HILDA Waves 12-24. Individual + wave FE. Cluster: individual. 95% CI.") +
  theme_paper
if (FALSE) {  # fig_hw_within_fe removed from the paper; block disabled
ggplot2::ggsave(file.path(outdir, "fig_hw_within_fe.pdf"), p_hw, width = 13/2.54, height = 9/2.54, units = "in", device = .pdf_device)
cat("  figures/fig_hw_within_fe.pdf\n")
}

# ============================================================================
# PART D -- Appendix: actual (ehtjb) vs Mincer-potential experience
# ============================================================================
cat("\n[D] Experience-measure comparison (actual vs Mincer potential) ...\n")
if (exists("wage")) {
  exp_all <- wage %>%
    dplyr::filter(!is.na(actual_experience), !is.na(educ_years), !is.na(age)) %>%
    dplyr::mutate(
      # recompute the Mincer proxy inline: the stored potential_experience column
      # is NA for W15-24 (the rows were stacked without it), so use age/educ_years,
      # which are present for all waves. Matches the loading definition.
      pot_exp  = pmax(0, age - educ_years - 5),
      educ_grp = factor(ifelse(educ_years > 12, "Post-secondary",
                               "Secondary or less"),
                        levels = c("Secondary or less", "Post-secondary")),
      gap = actual_experience - pot_exp)
  
  # share with actual > potential (a Mincer ceiling would truncate it), per group
  cap_lab <- exp_all %>%
    dplyr::group_by(educ_grp) %>%
    dplyr::summarise(pct = mean(gap > 0) * 100, .groups = "drop") %>%
    dplyr::mutate(lab = sprintf("%.0f%% above 0\n(Mincer would cap)", pct))
  
  exp_disp <- dplyr::filter(exp_all, gap >= -30, gap <= 20)   # trim tails for display only
  
  p_exp <- ggplot(exp_disp, aes(x = gap)) +
    geom_density(fill = "#B2182B", colour = "#B2182B", alpha = 0.30, linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    facet_wrap(~ educ_grp) +
    geom_text(data = cap_lab, aes(x = 11, y = Inf, label = lab),
              vjust = 1.4, hjust = 0.5, size = 2.3, colour = "grey20",
              lineheight = 0.9, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(-30, 20, 10)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    labs(title = "Actual minus Mincer-potential experience, by education",
         subtitle = "Positive values: actual work history exceeds the proxy, so a Mincer ceiling would truncate it",
         x = "Actual (ehtjb) - Mincer potential (years)", y = "Density") +
    theme_paper
  if (FALSE) ggplot2::ggsave("figures/figA1_experience_measures.pdf", p_exp,  # not a paper figure
                             width = 16/2.54, height = 8/2.54, units = "in")
  cat("  SKIPPED figA1_experience_measures (not a paper figure)\n")
} else cat("  SKIPPED figA1_experience_measures (wage not in workspace)\n")

cat("\nDone. All paper figures (except CDAG) written to output/figures/.\n")

# =====================================================================
# fig_event_study.R
# Figure 3: women-minus-men event-study paths for log wage and
# experience around first birth. Styled to match RUN_FIGURES.R
# (theme_paper, gap-series colour, save size) so it matches the rest of
# the figures in cost_hh_spec.pdf. Endpoints binned at |event time| = K.
# =====================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
# --- Namespace protection (avoid MASS/stats masking dplyr verbs) ----
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

# ---- shared paper style (lifted from RUN_FIGURES.R) ----
theme_paper <- theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    plot.title         = element_text(size = 7.5, face = "bold", hjust = 0),
    plot.subtitle      = element_text(size = 6,   colour = "grey40", hjust = 0),
    plot.caption       = element_text(size = 6,   colour = "grey50", hjust = 0),
    strip.text         = element_text(size = 6.5, face = "bold"),
    axis.title         = element_text(size = 6.5))
gap_col <- "#7570B3"   # gap-series colour, as in the housework-gap panel

K <- 8L

# cs_event_study_results.rds is produced by 02_cs_event_study.R (Sun-Abraham
# event study). If it is missing, build it now rather than aborting here,
# which previously killed the script before fig_event_study.pdf was saved.
if (!file.exists("cs_event_study_results.rds") && file.exists("02_cs_event_study.R")) {
  cat("  cs_event_study_results.rds not found; sourcing 02_cs_event_study.R to build it...\n")
  options(cs_estimation_only = TRUE)   # estimation only; paper-style figure built below
  source("02_cs_event_study.R")
  options(cs_estimation_only = NULL)
}
if (file.exists("cs_event_study_results.rds")) {
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
              outcome = "A. Log hourly wage"),
    transform(gap(res[["women experience_years"]],    res[["men experience_years"]]),
              outcome = "B. Experience (years)")) %>%
    mutate(lo = att - 1.96 * se, hi = att + 1.96 * se,
           outcome = factor(outcome, levels = c("A. Log hourly wage", "B. Experience (years)")))
  
  p <- ggplot(dat, aes(event_time, att)) +
    facet_wrap(~outcome, scales = "free_y") +
    geom_hline(yintercept = 0,    linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey40") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = gap_col, alpha = 0.15) +
    geom_line(colour = gap_col, linewidth = 0.9) +
    geom_point(colour = gap_col, size = 1.8) +
    labs(x = "Years from first birth", y = "ATT (women \u2212 men)") +
    theme_paper
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(outdir, "fig_event_study.pdf"), p, width = 13/2.54, height = 7/2.54, units = "in", device = .pdf_device)
  cat("  figures/fig_event_study.pdf\n")
} else {
  cat("  SKIPPED fig_event_study (cs_event_study_results.rds not found and",
      "02_cs_event_study.R not in working directory)\n")
}

# Remove the local masks installed at the top of this script so the session
# behaves normally afterwards (in particular the file.copy() no-op).
rm(file.copy, ggsave, save_fig)

options(error = NULL)   # remove the abort handler installed at the top
close_fig_log()