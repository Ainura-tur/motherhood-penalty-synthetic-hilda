# =============================================================================
# 00b_make_synthetic_hilda.R
#
# Emits SYNTHETIC versions of the three analysis panels that 00_load_hilda.R
# would build from the restricted HILDA release:
#
#   hilda_panel_data_W12_W24_slim.rds   (Parts A-F)
#   hilda_panel_data_extended.rds       (Parts G/H; waves 15-24)
#   hilda_panel_data_training.rds       (training event study)
#
# Purpose: let the Data Editor (or anyone without HILDA access) run run_all.R
# end to end and confirm the code EXECUTES. The data are random. They reproduce
# the column names, types and structural relationships the pipeline depends on
# (panel shape, absorbing parenthood, couple linkage, an education-by-experience
# wage gradient, FTB-B income/child-age variation), but they do NOT reproduce
# any published number. Every reported result requires the real HILDA release.
#
# Usage (from the repository root, the same folder run_all.R runs from):
#   Rscript 00b_make_synthetic_hilda.R              # default 4000 persons
#   N_PERSONS=8000 Rscript 00b_make_synthetic_hilda.R
#
# This writes the three .rds files into the working directory, exactly where
# 00_load_hilda.R would write them and where the analysis scripts read them.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(tidyr) })

# ── Namespace protection (avoid MASS/stats masking dplyr verbs) ───────────────
select <- dplyr::select;       filter    <- dplyr::filter
mutate <- dplyr::mutate;       slice     <- dplyr::slice
recode <- dplyr::recode;       rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange;     count     <- dplyr::count;  lag <- dplyr::lag

set.seed(20260618)

n_persons <- suppressWarnings(as.integer(Sys.getenv("N_PERSONS", "4000")))
if (is.na(n_persons) || n_persons < 500L) n_persons <- 4000L
WAVES <- 12:24

cat(strrep("=", 70), "\n", sep = "")
cat("  00b_make_synthetic_hilda.R  --  SYNTHETIC DATA (results are not real)\n")
cat(sprintf("  %d persons x %d waves\n", n_persons, length(WAVES)))
cat(strrep("=", 70), "\n", sep = "")

# WPI deflators, base wave 24 = 2024 (identical to 00_load_hilda.R / the master)
wpi <- 149.6 / c(110.9, 114.6, 117.6, 120.4, 123.0, 125.4, 127.9,
                 130.9, 133.7, 135.7, 138.9, 143.7, 149.6)
names(wpi) <- as.character(WAVES)

# ── 1. PERSON-LEVEL ATTRIBUTES ────────────────────────────────────────────────
# person_id must be >= 5 characters once stringified, because the master cleans
# partner_id with `nchar(partner_id) < 5 ~ NA` and links couples by matching
# partner_id against person_id. In 00_load_hilda.R person_id is
# as.character(xwaveid), so person_id, xwaveid and partner_id are all CHARACTER
# and the B1 couple join matches character to character. We mirror that here:
# integer pid only indexes the couple loop; the saved ids are 7-char strings.
pid <- 100000L + seq_len(n_persons)            # integer, internal use only

# Couple structure: pair the first 2*n_couples persons into opposite-sex
# couples; the rest are single. Couple members share a household_id and each
# carries the other's person_id (as a character string) as partner_id.
n_couples <- floor(n_persons * 0.55 / 2)
partner_pid_int <- rep(NA_integer_, n_persons)
household_int <- 500000L + seq_len(n_persons)
female <- rbinom(n_persons, 1, 0.5)
for (k in seq_len(n_couples)) {
  a <- 2L * k - 1L; b <- 2L * k
  partner_pid_int[a] <- pid[b]; partner_pid_int[b] <- pid[a]
  household_int[a] <- household_int[b] <- 600000L + k
  female[a] <- 1L; female[b] <- 0L          # opposite-sex couple
}

# Character ids (person_id == xwaveid, exactly as the loader builds them).
person_chr  <- sprintf("%07d", pid)
partner_chr <- ifelse(is.na(partner_pid_int), NA_character_, sprintf("%07d", partner_pid_int))
household_chr <- sprintf("HH%06d", household_int)

educ_years_p <- pmin(20, pmax(9, round(rnorm(n_persons, 14, 2.4))))
age12_p      <- sample(22:56, n_persons, replace = TRUE)     # age at wave 12
u_person     <- rnorm(n_persons, 0, 0.25)                    # wage heterogeneity
state_p      <- sample(1:8, n_persons, replace = TRUE,
                       prob = c(.32,.26,.20,.07,.10,.02,.01,.02))

# First-birth wave. Mix of: never-parent, left-censored (birth before wave 12),
# and birth within the window (needed for the event studies to have pre/post).
birth_type <- sample(c("never", "censored", "in_window"), n_persons,
                     replace = TRUE, prob = c(0.34, 0.30, 0.36))
first_birth_wave <- rep(NA_integer_, n_persons)
first_birth_wave[birth_type == "censored"]  <-
  sample(2:11, sum(birth_type == "censored"), replace = TRUE)
first_birth_wave[birth_type == "in_window"] <-
  sample(13:22, sum(birth_type == "in_window"), replace = TRUE)   # leaves >=2 post waves
# A possible second child 2-5 waves after the first.
second_birth_wave <- ifelse(!is.na(first_birth_wave) & runif(n_persons) < 0.5,
                            first_birth_wave + sample(2:5, n_persons, replace = TRUE),
                            NA_integer_)

person <- tibble(
  person_id = person_chr, xwaveid = person_chr, female = as.integer(female),
  educ_years = as.numeric(educ_years_p), age12 = age12_p,
  u = u_person, state_code = state_p,
  partner_pid = partner_chr, household_id = household_chr,
  fb = first_birth_wave, sb = second_birth_wave,
  mother_educ_years = pmin(20, pmax(7, round(rnorm(n_persons, 12, 2.5)))),
  father_educ_years = pmin(20, pmax(7, round(rnorm(n_persons, 12, 2.6))))
)

# ── 2. EXPAND TO PERSON x WAVE AND DERIVE ─────────────────────────────────────
panel <- tidyr::crossing(person, wave = WAVES) %>%
  arrange(person_id, wave) %>%
  mutate(
    year = wave + 2000L,
    age  = age12 + (wave - 12L),
    age_sq = age^2,
    deflator = as.numeric(wpi[as.character(wave)]),
    
    # Children: child ages from the (synthetic) birth waves.
    child1_age = ifelse(!is.na(fb) & wave >= fb, wave - fb, NA_integer_),
    child2_age = ifelse(!is.na(sb) & wave >= sb, wave - sb, NA_integer_),
    n_born = (!is.na(child1_age)) + (!is.na(child2_age)),
    num_children_under15 = ((!is.na(child1_age) & child1_age <= 14) +
                              (!is.na(child2_age) & child2_age <= 14)),
    num_children_0to4  = ((!is.na(child1_age) & child1_age <= 4) +
                            (!is.na(child2_age) & child2_age <= 4)),
    num_children_5to9  = ((!is.na(child1_age) & child1_age %in% 5:9) +
                            (!is.na(child2_age) & child2_age %in% 5:9)),
    num_children_10to14 = ((!is.na(child1_age) & child1_age %in% 10:14) +
                             (!is.na(child2_age) & child2_age %in% 10:14)),
    num_children_15to24 = ((!is.na(child1_age) & child1_age %in% 15:24) +
                             (!is.na(child2_age) & child2_age %in% 15:24)),
    num_dependent_children = num_children_under15,
    num_own_resident_children = n_born,
    num_own_nonresident_children = 0L,
    total_own_children = n_born,
    total_children_ever_had = n_born,
    total_children_died = 0L,
    num_step_foster_grand = 0L,
    age_resident_child_1 = child1_age,
    age_resident_child_2 = child2_age,
    age_resident_child_3 = NA_integer_,
    age_nonresident_child_1 = NA_integer_,
    age_nonresident_child_2 = NA_integer_,
    age_youngest_own_child = pmin(child1_age, child2_age, na.rm = TRUE),
    own_res_children_0to4  = num_children_0to4,
    own_res_children_5to14 = num_children_5to9 + num_children_10to14,
    own_res_children_15to24 = num_children_15to24,
    own_res_children_25plus = 0L,
    own_nonres_children_0to4 = 0L, own_nonres_children_5to14 = 0L,
    own_nonres_children_15to24 = 0L, own_nonres_children_25plus = 0L,
    has_nonresident_children_under18 = FALSE,
    has_any_nonresident_children = FALSE,
    
    # Absorbing parenthood (mirrors 00_load_hilda.R logic).
    ever_parent = as.integer(!is.na(fb) & wave >= fb),
    ever_mother = ifelse(female == 1, ever_parent, NA_integer_),
    ever_father = ifelse(female == 0, ever_parent, NA_integer_),
    motherhood_status = case_when(
      female == 0 ~ NA_character_,
      ever_parent == 1 & num_children_under15 > 0  ~ "active_mother",
      ever_parent == 1 & num_children_under15 == 0 ~ "post_mother",
      ever_parent == 0 ~ "never_mother", TRUE ~ NA_character_),
    fatherhood_status = case_when(
      female == 1 ~ NA_character_,
      ever_parent == 1 & num_children_under15 > 0  ~ "active_father",
      ever_parent == 1 & num_children_under15 == 0 ~ "post_father",
      ever_parent == 0 ~ "never_father", TRUE ~ NA_character_),
    age_oldest_child = pmax(age_resident_child_1, age_resident_child_2,
                            age_resident_child_3, na.rm = TRUE),
    approx_age_at_first_birth = ifelse(!is.na(age_oldest_child) & ever_parent == 1,
                                       age - age_oldest_child, NA_real_),
    has_young_children = as.integer(num_children_0to4 > 0),
    active_caregiver   = as.integer(num_dependent_children > 0),
    is_mother = as.integer(female == 1 & num_dependent_children > 0),
    is_father = as.integer(female == 0 & num_dependent_children > 0),
    
    group_ever = case_when(
      female == 0 & ever_father == 1 ~ "fathers_ever",
      female == 0 & ever_father == 0 ~ "childless_men_ever",
      female == 1 & ever_mother == 1 ~ "mothers_ever",
      female == 1 & ever_mother == 0 ~ "never_mothers", TRUE ~ NA_character_),
    group_current = case_when(
      female == 0 & num_children_under15 > 0  ~ "fathers",
      female == 0 & num_children_under15 == 0 ~ "childless_men",
      female == 1 & num_children_under15 > 0  ~ "mothers",
      female == 1 & num_children_under15 == 0 ~ "childless_women", TRUE ~ NA_character_),
    group_3way_women = case_when(
      female == 0 ~ NA_character_,
      motherhood_status == "never_mother"  ~ "never_mothers",
      motherhood_status == "active_mother" ~ "active_mothers",
      motherhood_status == "post_mother"   ~ "post_mothers", TRUE ~ NA_character_),
    
    # Birth events (the wave the absorbing indicator turns on).
    had_birth_this_wave = as.integer(!is.na(fb) & wave == fb),
    event_birth = had_birth_this_wave,
    
    # Experience, employment, hours.
    experience_years = pmax(0, age - educ_years - 6 -
                              ifelse(female == 1 & !is.na(fb) & wave >= fb, 1.2 * (wave - fb), 0)),
    experience_sq = experience_years^2,
    actual_experience = experience_years,
    tenure_employer = pmin(experience_years, pmax(0, round(rnorm(n(), 6, 4)))),
    tenure_occupation = pmin(experience_years, pmax(0, round(rnorm(n(), 8, 5)))),
    
    employed = rbinom(n(), 1,
                      plogis(1.6 - 0.9 * (female == 1 & num_children_0to4 > 0))),
    fulltime = as.integer(employed == 1 &
                            rbinom(n(), 1, plogis(0.8 - 1.4 * (female == 1 & num_children_under15 > 0))) == 1),
    parttime = as.integer(employed == 1 & fulltime == 0),
    married  = as.integer(!is.na(partner_pid)),
    married_num = as.numeric(married),
    fulltime_num = as.numeric(fulltime),
    esdtl_raw = ifelse(employed == 1, ifelse(fulltime == 1, 1L, 2L), 3L),
    
    hours_worked_clean = ifelse(employed == 0, NA_real_,
                                pmax(1, pmin(60, round(ifelse(fulltime == 1, 40, 22) +
                                                         rnorm(n(), 0, 5) - 8 * (female == 1 & num_children_under15 > 0))))),
    
    # Wage with an education-by-experience gradient and a motherhood penalty.
    ln_real_w = 2.55 + 0.060 * educ_years + 0.020 * experience_years +
      0.0015 * educ_years * experience_years -
      0.10 * female - 0.06 * female * ever_parent + u + rnorm(n(), 0, 0.30),
    real_hourly_wage = ifelse(employed == 0, NA_real_,
                              pmin(300, pmax(3, exp(ln_real_w)))),
    hourly_wage_clean = real_hourly_wage / deflator,            # nominal
    ln_hourly_wage = ifelse(!is.na(hourly_wage_clean) & hourly_wage_clean > 0,
                            log(hourly_wage_clean), NA_real_),
    ln_hourly_wage_real = ifelse(!is.na(real_hourly_wage) & real_hourly_wage > 0,
                                 log(real_hourly_wage), NA_real_),
    weekly_earnings = ifelse(employed == 0, 0,
                             hourly_wage_clean * hours_worked_clean),
    annual_income = pmax(0, weekly_earnings * 52 + rnorm(n(), 0, 4000)),
    childcare_cost = ifelse(num_children_0to4 > 0,
                            pmax(0, rnorm(n(), 8000, 4000)), 0)
  )

# Education factors (reuse 00_load_hilda.R's exact case_when and level order).
edu_level_order <- c("Other","Primary_Incomplete","Primary_Complete",
                     "Year_7","Year_8","Year_9","Year_10","Year_11","Year_12",
                     "Certificate_I","Certificate_II","Certificate_III","Certificate_IV",
                     "Diploma_Other","Diploma","Advanced_Diploma","Bachelor","Bachelor_Honours",
                     "Graduate_Diploma_Certificate","Masters","Doctorate")
edu_pathway_order <- c("Primary_Incomplete","Primary_Complete","Year_7","Year_8",
                       "Year_9","Below_Year_10","Year_10","Year_11","Year_12","Certificate_II",
                       "Cert_III_Trade","Cert_III_Technical","Cert_III_Nursing","Certificate_III",
                       "Cert_IV_Trade","Cert_IV_Technical","Certificate_IV","Diploma_Trade",
                       "Diploma_Teaching","Diploma_Nursing","Diploma","Advanced_Diploma",
                       "AdvDiploma_Trade","Bachelor_General","Bachelor_Honours_General",
                       "Bachelor_Nursing","Bachelor_Teaching","GradDipCert_General",
                       "GradDipCert_Nursing","GradDipCert_Teaching","Masters_General",
                       "Masters_Teaching","Doctorate_General","Doctorate_Teaching","Registered_Nurse")

panel <- panel %>%
  mutate(
    educ_sq = educ_years^2,
    education_level = case_when(
      educ_years >= 19.5 ~ "Doctorate",
      educ_years >= 18   ~ "Masters",
      educ_years >= 17   ~ "Graduate_Diploma_Certificate",
      educ_years >= 16.5 ~ "Bachelor_Honours",
      educ_years >= 16   ~ "Bachelor",
      educ_years >= 15   ~ "Advanced_Diploma",
      educ_years >= 14.5 ~ "Diploma",
      educ_years >= 14   ~ "Diploma_Other",
      educ_years >= 13.5 ~ "Certificate_IV",
      educ_years >= 13   ~ "Certificate_III",
      educ_years >= 12.5 ~ "Certificate_II",
      educ_years == 12   ~ "Year_12",
      educ_years > 12 & educ_years < 12.5 ~ "Certificate_I",
      educ_years == 11 ~ "Year_11", educ_years == 10 ~ "Year_10",
      educ_years == 9 ~ "Year_9", TRUE ~ "Other"),
    education_pathway = dplyr::case_when(
      educ_years >= 16 ~ "Bachelor_General",
      educ_years >= 14 ~ "Diploma",
      educ_years >= 13 ~ "Certificate_III",
      educ_years == 12 ~ "Year_12", TRUE ~ "Year_10"),
    education_data_source = "synthetic",
    has_degree        = as.integer(educ_years >= 16),
    has_postsecondary = as.integer(educ_years > 12),
    has_diploma_plus  = as.integer(educ_years >= 14),
    has_postgrad      = as.integer(educ_years >= 17),
    has_nursing_any = 0L, has_teaching_any = 0L, has_business_any = 0L,
    has_trade_any = as.integer(educ_years %in% 13:14), has_technical_any = 0L,
    
    # Job characteristics / quality.
    casual = as.integer(employed == 1 & parttime == 1 & runif(n()) < 0.4),
    permanent = as.integer(employed == 1 & casual == 0),
    public_sector = rbinom(n(), 1, 0.25),
    union_member  = rbinom(n(), 1, 0.20),
    firm_size = sample(c(5, 25, 75, 250, 1000), n(), replace = TRUE),
    has_holiday_leave = as.integer(employed == 1 & permanent == 1),
    has_sick_leave    = has_holiday_leave,
    leave_none = as.integer(employed == 1 & permanent == 0),
    has_standard_hours = rbinom(n(), 1, 0.7),
    has_shift_work     = rbinom(n(), 1, 0.2),
    has_irregular_hours = rbinom(n(), 1, 0.2),
    underemployed = as.integer(parttime == 1 & runif(n()) < 0.2),
    overemployed  = as.integer(fulltime == 1 & runif(n()) < 0.2),
    precarious_score = round(runif(n(), 0, 5), 1),
    precarious_category = cut(precarious_score, c(-Inf, 1.5, 3, Inf),
                              labels = c("Low", "Medium", "High")),
    flexibility_any = rbinom(n(), 1, 0.5),
    treatment = flexibility_any,
    
    job_satisfaction = round(runif(n(), 0, 10)),
    pay_fairness     = round(runif(n(), 0, 10)),
    job_initiative   = round(runif(n(), 0, 10)),
    job_security     = round(runif(n(), 0, 10)),
    job_worry_future = round(runif(n(), 0, 10)),
    job_autonomy     = round(runif(n(), 0, 10)),
    looking_for_new_job = rbinom(n(), 1, 0.2),
    
    occupation_major_code = sample(1:8, n(), replace = TRUE),
    occupation_major = recode(as.character(occupation_major_code),
                              "1"="Managers","2"="Professionals","3"="Technicians_Trades",
                              "4"="Community_Service","5"="Clerical_Admin","6"="Sales",
                              "7"="Machinery_Operators","8"="Labourers"),
    industry_major_code = sample(1:10, n(), replace = TRUE),
    industry_major = paste0("Industry_", industry_major_code),
    
    # Training (only meaningful when employed).
    any_training = as.integer(employed == 1 & runif(n()) < 0.4),
    total_training_hours = ifelse(any_training == 1, round(runif(n(), 1, 60)), 0),
    promotion_training = as.integer(any_training == 1 & runif(n()) < 0.5),
    maintenance_training = as.integer(any_training == 1 & runif(n()) < 0.5),
    total_training_days = ifelse(any_training == 1, round(runif(n(), 1, 15)), 0),
    avg_training_hours  = ifelse(any_training == 1, round(runif(n(), 1, 8), 1), 0),
    total_courses = ifelse(any_training == 1, sample(1:4, n(), replace = TRUE), 0L),
    training_own_time  = as.integer(any_training == 1 & runif(n()) < 0.4),
    training_paid_time = as.integer(any_training == 1 & runif(n()) < 0.6),
    jttpewt_any = training_paid_time,
    jttpeot_any = training_own_time,
    
    # Social capital (created downstream in the real pipeline; included as a
    # harmless superset so direct references never fail).
    club_member = rbinom(n(), 1, 0.3),
    social_freq = sample(1:6, n(), replace = TRUE),
    social_freq_01 = round(social_freq / 6, 2),
    high_social = as.integer(social_freq >= 4),
    networking_index = round(runif(n(), 0, 1), 2),
    
    housework_hours = pmax(0, round(rnorm(n(),
                                          8 + 6 * (female == 1) + 4 * (num_children_under15 > 0), 4))),
    care_intensity = num_children_under15 * 10,
    
    # Parenting attitudes (fertility-module style; 1-7 Likert, NA when not parent).
    parenting_resp = ifelse(ever_parent == 1, sample(1:7, n(), replace = TRUE), NA_integer_),
    parenting_hard = ifelse(ever_parent == 1, sample(1:7, n(), replace = TRUE), NA_integer_),
    parenting_tired = ifelse(ever_parent == 1, sample(1:7, n(), replace = TRUE), NA_integer_),
    parenting_trapped = ifelse(ever_parent == 1, sample(1:7, n(), replace = TRUE), NA_integer_),
    parenting_fairshare = ifelse(ever_parent == 1, sample(1:7, n(), replace = TRUE), NA_integer_),
    
    currently_pregnant = as.integer(female == 1 & !is.na(fb) & wave == fb - 1L),
    partner_pregnant   = as.integer(female == 0 & !is.na(fb) & wave == fb - 1L),
    self_or_partner_pregnant = as.integer(currently_pregnant == 1 | partner_pregnant == 1),
    expected_birth_month = ifelse(self_or_partner_pregnant == 1,
                                  sample(1:12, n(), replace = TRUE), NA_integer_),
    expected_birth_year  = ifelse(self_or_partner_pregnant == 1, year + 1L, NA_integer_),
    
    event_married   = as.integer(wave == 13 & married == 1 & runif(n()) < 0.1),
    event_separated = rbinom(n(), 1, 0.03),
    event_job_change = rbinom(n(), 1, 0.12),
    event_fired = rbinom(n(), 1, 0.03),
    event_moved = rbinom(n(), 1, 0.15),
    got_married_this_wave = event_married,
    separated_this_wave   = event_separated,
    changed_job_this_wave = event_job_change,
    
    other_parent_employment = ifelse(!is.na(partner_pid), rbinom(n(), 1, 0.85), NA_integer_),
    living_with_partner_last = as.integer(!is.na(partner_pid)),
    
    state = recode(as.character(state_code),
                   "1"="NSW","2"="VIC","3"="QLD","4"="SA","5"="WA","6"="TAS","7"="NT","8"="ACT"),
    period = case_when(year %in% 2015:2019 ~ "Pre_COVID",
                       year %in% 2020:2021 ~ "COVID",
                       year %in% 2022:2024 ~ "Post_COVID", TRUE ~ NA_character_),
    
    analysis_weight     = round(rlnorm(n(), log(1), 0.4), 3),
    population_weight    = analysis_weight * 1000,
    longitudinal_weight = round(rlnorm(n(), log(1), 0.4), 3),
    
    in_wage_sample = as.integer(employed == 1 & !is.na(ln_hourly_wage_real) &
                                  !is.na(hours_worked_clean) & hours_worked_clean > 0 & hours_worked_clean <= 100 &
                                  experience_years >= 0 & analysis_weight > 0)
  )

# Interactions (after the building blocks exist).
panel <- panel %>%
  mutate(
    educ_centered = educ_years - mean(educ_years, na.rm = TRUE),
    educ_centered_sq = educ_centered^2,
    educ_female = educ_years * female,
    educ_mother = educ_years * is_mother,
    educ_treatment = educ_years * treatment,
    educ_precarious = educ_years * precarious_score,
    educ_ever_mother = educ_years * ifelse(is.na(ever_mother), 0, ever_mother),
    educ_ever_father = educ_years * ifelse(is.na(ever_father), 0, ever_father)
  )

# Sparse third (older) resident child for ~8% of persons who have a second
# child, assigned at the person level so it is consistent across that person's
# waves. Keeps age_resident_child_3 and the 15-24 resident-child counts from
# being uniformly empty. Non-resident children are not modelled (their flags
# stay FALSE and counts 0), so their age columns remain NA by construction.
panel <- panel %>%
  mutate(
    .c3 = (as.integer(person_id) %% 12L == 0L) & !is.na(age_resident_child_2),
    age_resident_child_3 = ifelse(.c3, pmin(24L, age_resident_child_2 + 2L),
                                  age_resident_child_3),
    own_res_children_15to24   = own_res_children_15to24 + as.integer(.c3),
    num_own_resident_children = num_own_resident_children + as.integer(.c3),
    total_own_children        = total_own_children + as.integer(.c3),
    total_children_ever_had   = total_children_ever_had + as.integer(.c3)
  ) %>% select(-.c3)

# Identifiers and couple keys. person_id and xwaveid are already the same 7-char
# string (as in the loader). partner_id carries the partner's person_id string,
# which survives the master's nchar(.) >= 5 / no-leading-"-" cleaning.
panel <- panel %>%
  mutate(
    xwaveid = person_id,
    partner_id = partner_pid,
    household_id = as.character(household_id)
  )

# Apply the education factor levels exactly as the loader does.
panel$education_level <- factor(panel$education_level, levels = edu_level_order, ordered = TRUE)
panel$education_pathway <- factor(panel$education_pathway, levels = edu_pathway_order, ordered = TRUE)

# ── 3. ORDER COLUMNS AS 00_load_hilda.R WRITES THEM, THEN SAVE ────────────────
col_order <- c(
  "person_id","xwaveid","wave","year","age","age_sq","female","married",
  "ln_hourly_wage","ln_hourly_wage_real","hourly_wage_clean","real_hourly_wage",
  "deflator","weekly_earnings","hours_worked_clean","annual_income","childcare_cost",
  "educ_years","educ_sq","education_level","education_pathway","has_degree",
  "has_postsecondary","has_diploma_plus","has_postgrad","has_nursing_any",
  "has_teaching_any","has_business_any","has_trade_any","has_technical_any",
  "education_data_source","experience_years","experience_sq","actual_experience",
  "tenure_employer","tenure_occupation","employed","fulltime","parttime",
  "married_num","fulltime_num","esdtl_raw","in_wage_sample","casual","permanent",
  "public_sector","union_member","firm_size","has_holiday_leave","has_sick_leave",
  "leave_none","has_standard_hours","has_shift_work","has_irregular_hours",
  "underemployed","overemployed","precarious_score","precarious_category",
  "flexibility_any","treatment","job_satisfaction","pay_fairness","job_initiative",
  "job_security","job_worry_future","job_autonomy","looking_for_new_job",
  "occupation_major","occupation_major_code","industry_major","industry_major_code",
  "any_training","total_training_hours","promotion_training","maintenance_training",
  "total_training_days","avg_training_hours","total_courses","training_own_time",
  "training_paid_time","jttpewt_any","jttpeot_any","club_member","social_freq",
  "social_freq_01","high_social","networking_index","num_dependent_children",
  "num_children_under15","num_children_0to4","num_children_5to9","num_children_10to14",
  "num_children_15to24","has_young_children","active_caregiver","is_mother","is_father",
  "housework_hours","care_intensity","has_nonresident_children_under18",
  "has_any_nonresident_children","total_children_ever_had","num_own_resident_children",
  "num_own_nonresident_children","total_own_children","age_youngest_own_child",
  "total_children_died","num_step_foster_grand","own_res_children_0to4",
  "own_res_children_5to14","own_res_children_15to24","own_res_children_25plus",
  "own_nonres_children_0to4","own_nonres_children_5to14","own_nonres_children_15to24",
  "own_nonres_children_25plus","age_resident_child_1","age_resident_child_2",
  "age_resident_child_3","age_nonresident_child_1","age_nonresident_child_2",
  "ever_parent","ever_mother","ever_father","motherhood_status","fatherhood_status",
  "approx_age_at_first_birth","age_oldest_child","group_ever","group_current",
  "group_3way_women","parenting_resp","parenting_hard","parenting_tired",
  "parenting_trapped","parenting_fairshare","currently_pregnant","partner_pregnant",
  "self_or_partner_pregnant","expected_birth_month","expected_birth_year","event_birth",
  "event_married","event_separated","event_job_change","event_fired","event_moved",
  "had_birth_this_wave","got_married_this_wave","separated_this_wave",
  "changed_job_this_wave","other_parent_employment","living_with_partner_last",
  "household_id","partner_id","mother_educ_years","father_educ_years","state",
  "analysis_weight","population_weight","longitudinal_weight","educ_centered",
  "educ_centered_sq","educ_female","educ_mother","educ_treatment","educ_precarious",
  "educ_ever_mother","educ_ever_father","period"
)

analysis_sample <- panel %>% select(any_of(col_order))
missing_cols <- setdiff(col_order, names(analysis_sample))
if (length(missing_cols))
  stop("Generator is missing columns: ", paste(missing_cols, collapse = ", "))

# Clean the +/-Inf that pmin/pmax(na.rm=TRUE) can produce when both args are NA.
analysis_sample <- analysis_sample %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))

stopifnot(min(analysis_sample$wave) <= 12, max(analysis_sample$wave) >= 24)

saveRDS(analysis_sample, "hilda_panel_data_W12_W24_slim.rds")
cat(sprintf("Saved hilda_panel_data_W12_W24_slim.rds: %s rows x %d cols\n",
            format(nrow(analysis_sample), big.mark = ","), ncol(analysis_sample)))

ext_sample <- filter(analysis_sample, wave >= 15)
saveRDS(ext_sample, "hilda_panel_data_extended.rds")
cat(sprintf("Saved hilda_panel_data_extended.rds: %s rows (waves 15-24)\n",
            format(nrow(ext_sample), big.mark = ",")))

# Training panel: the exact training_vars subset 00_load_hilda.R writes.
training_vars <- c("person_id","xwaveid","wave","year","age","female","married",
                   "ever_parent","employed","hours_worked_clean","ln_hourly_wage_real",
                   "experience_years","analysis_weight","longitudinal_weight","total_training_days",
                   "avg_training_hours","total_courses","training_own_time","training_paid_time",
                   "jttpewt_any","jttpeot_any")
int_vars <- c("wave","year","female","married","ever_parent","employed",
              "total_courses","training_own_time","training_paid_time","jttpewt_any","jttpeot_any")
training_panel <- analysis_sample %>%
  select(any_of(training_vars)) %>%
  mutate(across(any_of(int_vars), ~ as.integer(round(as.numeric(.x)))))
saveRDS(training_panel, "hilda_panel_data_training.rds", compress = "xz")
cat(sprintf("Saved hilda_panel_data_training.rds: %d rows x %d cols\n",
            nrow(training_panel), ncol(training_panel)))

cat("\nSYNTHETIC build complete. These files let run_all.R execute; the numbers\n")
cat("they produce are random and must not be reported.\n")