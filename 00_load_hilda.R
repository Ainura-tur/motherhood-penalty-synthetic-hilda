# Clear workspace
rm(list = ls())
gc()

# =====================================================================
# HILDA_LOADER_UNIFIED.R -- SINGLE LOADER FOR ALL ANALYSIS DATASETS
# Waves 12-24 (2012-2024), one pass over the Combined files, writes:
#   hilda_panel_data_W12_W24_slim.rds  full select, waves 12-24
#                                      (input to MASTER_hh.R Parts A-F)
#   hilda_panel_data_extended.rds      same select, waves 15-24
#                                      (input to MASTER_PART_G / Part H)
#   hilda_panel_data_training.rds      slim training extract, waves 12-24
#                                      (input to training_event_study.R)
# Part G's fallback file hilda_panel_data_W12_W24.rds is no longer
# needed (extended.rds always exists after this run).
#
# Lineage: supersedes HILDA_LOADIND_MASTER.R (W15-24, v5) and
# LOADIND_HILDA_MAIN_wo_partner.R (W15-24, pre-fix). Carries the v5
# corrections: NA-propagating potential_experience and employed;
# biological experience ceiling (age-13); -5 schooling offset;
# ==1 training indicators; column-restricted reads; versioned wave
# checkpoints. The W12-W24 slim panel it writes is the CANONICAL
# rebuild: it is NOT guaranteed byte-identical to the historical slim
# panel (whose loader is not in the repo), so after the first run,
# re-verify the MASTER FINAL VERIFICATION counts against the tex.
# =====================================================================
WAVES <- 12:24   # full analysis window; wave 12 = 2012

# --- Raw-data path comes from config.R (sourced after the rm() above, so the
#     assignment survives). config.R defines data_dir; run from the repo root. ---
if (file.exists("config.R")) source("config.R") else stop(
  "config.R not found in the working directory. Run 00_load_hilda.R from the ",
  "repository root (the folder that holds config.R and the scripts).")



library(tidyverse)
library(haven)
# truncnorm: legacy import, no longer used (no random imputation
# anywhere in this loader; output is fully deterministic)


# =====================================================================
# SETUP: WPI DATA FOR REAL WAGES  ABS A2705194A in 2023 dollars
# =====================================================================

# wpi_data <- data.frame(
#   wave = 15:23,
#   deflator = 143.0 / c(122.7, 125.1, 127.7, 130.6, 133.5, 135.4, 138.5, 143.2, 149.3)
# )


# WPI index (ABS A2705194A), base wave 24 = 149.6. Waves 15-24 are the
# verified values used throughout the project. Waves 12-14 are ESTIMATES
# back-cast at published WPI growth (2012: 3.6%, 2013: 2.9%, 2014: 2.5%)
# -- VERIFY against the ABS series and replace before a production run.
WPI_12_14_VERIFIED <- TRUE   # ABS A2705194A, FY index (year ended June), confirmed
wpi_data <- data.frame(
  wave     = 12:24,
  deflator = 149.6 / c(110.9, 114.6, 117.6,                       # W12-14: ABS-confirmed
                       120.4, 123.0, 125.4, 127.9, 130.9,
                       133.7, 135.7, 138.9, 143.7, 149.6)
)
if (!WPI_12_14_VERIFIED)
  warning("WPI for waves 12-14 are back-cast ESTIMATES; verify against ",
          "ABS A2705194A and set WPI_12_14_VERIFIED <- TRUE.", call. = FALSE)
# =====================================================================
# TWO HELPER FUNCTIONS FOR DIFFERENT VARIABLE TYPES
# =====================================================================

# Function 1: For variables where zero is INVALID (wages, hours, codes)
safe_numeric <- function(data, varname) {
  if(varname %in% names(data)) {
    val <- as.numeric(data[[varname]])
    val[val >= -10 & val <= -1] <- NA  # HILDA missing codes
    return(val)
  } else {
    return(rep(NA_real_, nrow(data)))
  }
}

# Function 2: For variables where zero is VALID (counts, tenure, experience)
safe_numeric_zero <- function(data, varname) {
  if(varname %in% names(data)) {
    val <- as.numeric(data[[varname]])
    val[val >= -10 & val <= -1] <- NA  # HILDA missing codes only
    # Keep zeros as valid values
    return(val)
  } else {
    return(rep(NA_real_, nrow(data)))
  }
}

# =====================================================================
# MAIN PROCESSING FUNCTION WITH ALL VARIABLES
# =====================================================================

process_wave_complete <- function(wave_num, data_dir, wpi_data) {
  
  wave_letter <- letters[wave_num]
  year_val <- 2000 + wave_num
  cat("\n===========================================\n")
  cat("Processing Wave", wave_num, "(", wave_letter, ") - Year", year_val, "\n")
  cat("===========================================\n")
  
  # Find file. HILDA names embed the RELEASE, not the wave:
  # Combined_<waveletter><release>0c.dta (e.g. wave 12 in Release 24 is
  # Combined_l240c.dta). Match any release for this wave letter and take
  # the most recent if several releases coexist in the folder.
  cand <- list.files(data_dir,
                     pattern = paste0("^Combined_", wave_letter, "[0-9]{2,3}c\\.dta$"),
                     full.names = TRUE)
  file_path <- if (length(cand)) sort(cand, decreasing = TRUE)[1] else NULL
  if (!is.null(file_path)) cat("  Found file:", basename(file_path), "\n")
  
  if(is.null(file_path)) {
    cat("  ERROR: No file found for wave", wave_num, "\n")
    return(NULL)
  }
  
  # (Data read moved to AFTER the variable map `v` below, so only the
  #  columns this script actually uses are loaded from the ~6,000-column
  #  Combined file. See "Load data (column-restricted)".)
  
  # Get WPI deflator for this wave
  wave_deflator <- wpi_data$deflator[wpi_data$wave == wave_num]
  
  # ============================
  # COMPLETE VARIABLE MAPPING
  # ============================
  
  v <- list(
    # Demographics
    hgage = paste0(wave_letter, "hgage"),
    hgsex = paste0(wave_letter, "hgsex"),
    mrcurr = paste0(wave_letter, "mrcurr"),
    
    # Employment
    esdtl = paste0(wave_letter, "esdtl"),
    wscei = paste0(wave_letter, "wscei"),
    jbhruc = paste0(wave_letter, "jbhruc"),
    
    # Financial year income (for FTB simulator)
    wsfei = paste0(wave_letter, "wsfei"),   # FY wages/salary income (all jobs)
    wsfes = paste0(wave_letter, "wsfes"),   # FY self-employment income
    ccactci = paste0(wave_letter, "ccactci"),   # annual childcare cost
    
    # Education - Main
    edhigh1 = paste0(wave_letter, "edhigh1"),
    edhists = paste0(wave_letter, "edhists"),
    edsscmp = paste0(wave_letter, "edsscmp"),  # CPQ school completion
    
    # Education - Completed (edqro*)
    edqrodc = paste0(wave_letter, "edqrodc"),  # Doctorate completed
    edqroms = paste0(wave_letter, "edqroms"),  # Masters completed
    edqrogd = paste0(wave_letter, "edqrogd"),  # Graduate diploma completed
    edqrogc = paste0(wave_letter, "edqrogc"),  # Graduate certificate completed
    edqrohd = paste0(wave_letter, "edqrohd"),  # Honours bachelor completed
    edqrobd = paste0(wave_letter, "edqrobd"),  # Bachelor completed
    edqroad = paste0(wave_letter, "edqroad"),  # Associate degree completed
    edqroav = paste0(wave_letter, "edqroav"),  # Advanced diploma completed
    edqrodp = paste0(wave_letter, "edqrodp"),  # Diploma completed
    edqrodn = paste0(wave_letter, "edqrodn"),  # Diploma nursing completed
    edqroc4 = paste0(wave_letter, "edqroc4"),  # Certificate IV completed
    edqroc3 = paste0(wave_letter, "edqroc3"),  # Certificate III completed
    edqroc2 = paste0(wave_letter, "edqroc2"),  # Certificate II completed
    edqroc1 = paste0(wave_letter, "edqroc1"),  # Certificate I completed
    edqrocd = paste0(wave_letter, "edqrocd"),  # Certificate unknown level completed
    edqrosh = paste0(wave_letter, "edqrosh"),  # Secondary school highest completed
    edqrosl = paste0(wave_letter, "edqrosl"),  # Secondary school lower completed
    
    # Education - Obtained (edqo*)
    edqodc = paste0(wave_letter, "edqodc"),   # Doctorate obtained
    edqoms = paste0(wave_letter, "edqoms"),   # Masters obtained
    edqogd = paste0(wave_letter, "edqogd"),   # Graduate diploma obtained
    edqogc = paste0(wave_letter, "edqogc"),   # Graduate certificate obtained
    edqohd = paste0(wave_letter, "edqohd"),   # Honours bachelor obtained
    edqobd = paste0(wave_letter, "edqobd"),   # Bachelor obtained
    edqoad = paste0(wave_letter, "edqoad"),   # Associate degree obtained
    edqoav = paste0(wave_letter, "edqoav"),   # Advanced diploma obtained
    edqodp = paste0(wave_letter, "edqodp"),   # Diploma obtained
    edqodn = paste0(wave_letter, "edqodn"),   # Diploma nursing obtained
    edqoc4 = paste0(wave_letter, "edqoc4"),   # Certificate IV obtained
    edqoc3 = paste0(wave_letter, "edqoc3"),   # Certificate III obtained
    edqoc2 = paste0(wave_letter, "edqoc2"),   # Certificate II obtained
    edqoc1 = paste0(wave_letter, "edqoc1"),   # Certificate I obtained
    edqocd = paste0(wave_letter, "edqocd"),   # Certificate unknown level obtained
    edqota = paste0(wave_letter, "edqota"),   # Trade certificate obtained
    edqotc = paste0(wave_letter, "edqotc"),   # Technical certificate obtained
    edqonq = paste0(wave_letter, "edqonq"),   # Nursing qualification obtained
    edqotq = paste0(wave_letter, "edqotq"),   # Teaching qualification obtained
    edqobc = paste0(wave_letter, "edqobc"),   # Business qualification obtained
    edqosh = paste0(wave_letter, "edqosh"),   # Secondary school highest obtained
    edqosl = paste0(wave_letter, "edqosl"),   # Secondary school lower obtained 
    
    # Nursing-specific qualifications
    ednrsdc = paste0(wave_letter, "ednrsdc"), # Nursing doctorate
    ednrsms = paste0(wave_letter, "ednrsms"), # Nursing masters
    ednrsgd = paste0(wave_letter, "ednrsgd"), # Nursing graduate diploma
    ednrsbd = paste0(wave_letter, "ednrsbd"), # Nursing bachelor
    ednrshd = paste0(wave_letter, "ednrshd"), # Nursing honours
    ednrsav = paste0(wave_letter, "ednrsav"), # Nursing advanced diploma
    ednrsdp = paste0(wave_letter, "ednrsdp"), # Nursing diploma
    ednrsc3 = paste0(wave_letter, "ednrsc3"), # Nursing certificate III
    ednrstd = paste0(wave_letter, "ednrstd"), # Registered Nurse qualification
    
    # Teaching-specific qualifications  
    edtchdc = paste0(wave_letter, "edtchdc"), # Teaching doctorate
    edtchms = paste0(wave_letter, "edtchms"), # Teaching masters
    edtchgd = paste0(wave_letter, "edtchgd"), # Teaching graduate diploma
    edtchgc = paste0(wave_letter, "edtchgc"), # Teaching graduate certificate
    edtchhd = paste0(wave_letter, "edtchhd"), # Teaching honours
    edtchbd = paste0(wave_letter, "edtchbd"), # Teaching bachelor
    edtchts = paste0(wave_letter, "edtchts"), # Trained Secondary Teachers Cert
    edtchtp = paste0(wave_letter, "edtchtp"), # Primary Teaching Cert
    edtchdp = paste0(wave_letter, "edtchdp"), # Teaching diploma
    edtchav = paste0(wave_letter, "edtchav"), # Teaching advanced diploma
    
    #Training
    jttdays  = paste0(wave_letter, "jttdays"),   # Number of days attended training course
    jtthrs   = paste0(wave_letter, "jtthrs"),    # Average number of training hours per course day
    jttnum   = paste0(wave_letter, "jttnum"),    # Number of different training courses attended last 12 months
    jttopot  = paste0(wave_letter, "jttopot"),   # Job-related training at some other place - in your own time
    jttopwt  = paste0(wave_letter, "jttopwt"),   # Job-related training at some other place - during paid work time
    jttpeot  = paste0(wave_letter, "jttpeot"),   # Job-related training at place of employment - in your own time
    jttpewt  = paste0(wave_letter, "jttpewt"),    # Job-related training at place of employment - during paid work time
    
    # Experience and tenure
    ehtjb = paste0(wave_letter, "ehtjb"),
    jbempt = paste0(wave_letter, "jbempt"),
    jbocct = paste0(wave_letter, "jbocct"),
    
    # Job characteristics
    jbcasab = paste0(wave_letter, "jbcasab"),
    jbmhl = paste0(wave_letter, "jbmhl"),
    jbmsl = paste0(wave_letter, "jbmsl"),
    jbmsch = paste0(wave_letter, "jbmsch"),
    jbmday = paste0(wave_letter, "jbmday"),
    jbhrcpr = paste0(wave_letter, "jbhrcpr"),
    jbpmfhr = paste0(wave_letter, "jbpmfhr"),
    jbmtuea = paste0(wave_letter, "jbmtuea"),
    jbmmply = paste0(wave_letter, "jbmmply"),
    jbmemsz = paste0(wave_letter, "jbmemsz"),
    jbmems2 = paste0(wave_letter, "jbmems2"),
    jbptrea = paste0(wave_letter, "jbptrea"),
    
    # Occupation and Industry
    jbmo61 = paste0(wave_letter, "jbmo61"),
    jbmi61 = paste0(wave_letter, "jbmi61"),
    edpsfdn = paste0(wave_letter, "edpsfdn"),
    
    # Children - household counts (original)
    hhd0_4 = paste0(wave_letter, "hhd0_4"),
    hhd5_9 = paste0(wave_letter, "hhd5_9"),
    hhd1014 = paste0(wave_letter, "hhd1014"),
    hhd1524 = paste0(wave_letter, "hhd1524"),
    ncchk = paste0(wave_letter, "ncchk"),
    ncany = paste0(wave_letter, "ncany"),
    
    # ---- NEW: Fertility history (tc* prefix) ----
    # tchad: History variable — total children ever had (bio + adopted)
    #   Available most waves as History or PQ/NPQ. Absorbing: once > 0, stays > 0.
    tchad   = paste0(wave_letter, "tchad"),     # History: Total children ever had
    tchadn  = paste0(wave_letter, "tchadn"),    # Total children ever had (NPQ version)
    tchave  = paste0(wave_letter, "tchave"),     # Number of children (CPQ version)
    tchch   = paste0(wave_letter, "tchch"),      # Check G1: has ever had children (fertility module waves)
    
    # Own children counts (DV — available ALL waves)
    tcr     = paste0(wave_letter, "tcr"),        # DV: Number of own resident children
    tcnr    = paste0(wave_letter, "tcnr"),       # DV: Number of own non-resident children
    tcyng   = paste0(wave_letter, "tcyng"),      # DV: Age youngest own child (excl step/foster/grand)
    tcdied  = paste0(wave_letter, "tcdied"),     # History: Total children since died
    
    # Own resident children by age group (DV — all waves)
    tcr04   = paste0(wave_letter, "tcr04"),      # Resident own + step/foster/grand aged 0-4
    tcr514  = paste0(wave_letter, "tcr514"),     # Same, aged 5-14
    tcr1524 = paste0(wave_letter, "tcr1524"),    # Same, aged 15-24
    tcr25   = paste0(wave_letter, "tcr25"),      # Same, aged 25+
    
    # Own non-resident children by age group (DV — all waves)
    tcnr04  = paste0(wave_letter, "tcnr04"),     # Non-resident own children aged 0-4
    tcnr514 = paste0(wave_letter, "tcnr514"),    # Same, aged 5-14
    tcnr1524= paste0(wave_letter, "tcnr1524"),   # Same, aged 15-24
    tcnr25  = paste0(wave_letter, "tcnr25"),     # Same, aged 25+
    
    # Step/foster/grandchildren
    rcstepn = paste0(wave_letter, "rcstepn"),    # N resident step/foster/grandchildren, no parent in HH
    
    # Child ages (for computing age at first birth)
    rcage1  = paste0(wave_letter, "rcage1"),     # Age of resident child 1
    rcage2  = paste0(wave_letter, "rcage2"),     # Age of resident child 2
    rcage3  = paste0(wave_letter, "rcage3"),     # Age of resident child 3
    ncage1  = paste0(wave_letter, "ncage1"),     # Age of non-resident child 1
    ncage2  = paste0(wave_letter, "ncage2"),     # Age of non-resident child 2
    
    # ---- NEW: Parenting attitudes ----
    paresp  = paste0(wave_letter, "paresp"),     # Has parenting responsibilities, children <= 17
    pahard  = paste0(wave_letter, "pahard"),      # "Being a parent is harder than I thought"
    patird  = paste0(wave_letter, "patird"),      # "I often feel tired/exhausted from children's needs"
    patrap  = paste0(wave_letter, "patrap"),      # "I feel trapped by my responsibilities as a parent"
    pashare = paste0(wave_letter, "pashare"),     # "Do fair share of looking after children"
    
    # ---- NEW: Pregnancy (fertility module waves: 15, 19, 23 in our panel) ----
    ftcpg   = paste0(wave_letter, "ftcpg"),      # Currently pregnant
    ftcppg  = paste0(wave_letter, "ftcppg"),     # Partner currently pregnant
    ftipg   = paste0(wave_letter, "ftipg"),      # Check if self or partner currently pregnant
    ftcpgm  = paste0(wave_letter, "ftcpgm"),     # Expected birth month
    ftcpgy  = paste0(wave_letter, "ftcpgy"),     # Expected birth year
    
    # ---- NEW: Other parent / partner info ----
    ncesop  = paste0(wave_letter, "ncesop"),     # Employment status of the other parent
    prlpli  = paste0(wave_letter, "prlpli"),     # Living with partner when last interviewed
    hhrhid  = paste0(wave_letter, "hhrhid"),     # Household ID (changes when people move)
    hhpxid  = paste0(wave_letter, "hhpxid"),     # Partner's cross-wave person ID (direct link)
    
    # ---- NEW: Job satisfaction & workplace quality ----
    josat   = paste0(wave_letter, "josat"),       # Overall enjoyment and satisfaction from work
    jompf   = paste0(wave_letter, "jompf"),       # "I get paid fairly for things I do"
    jomini  = paste0(wave_letter, "jomini"),      # "My job requires me to take initiative"
    jomsf   = paste0(wave_letter, "jomsf"),       # "I have a secure future in my job"
    jomwf   = paste0(wave_letter, "jomwf"),       # "I worry about the future of my job"
    jomls_say = paste0(wave_letter, "jomls"),     # "I have a lot of say about what happens on my job"
    jbnewjs = paste0(wave_letter, "jbnewjs"),     # Looked for new job in last 4 weeks
    
    # ---- NEW: Life events (asked annually — all waves) ----
    leprg   = paste0(wave_letter, "leprg"),      # Life event: pregnancy/birth of new child
    lemar   = paste0(wave_letter, "lemar"),      # Life event: got married
    lesep   = paste0(wave_letter, "lesep"),      # Life event: separated from spouse
    lejob   = paste0(wave_letter, "lejob"),      # Life event: changed jobs
    lejlf   = paste0(wave_letter, "lejlf"),      # Life event: fired/made redundant
    lemvd   = paste0(wave_letter, "lemvd"),      # Life event: moved residence
    
    # Housework and care
    lshrhw = paste0(wave_letter, "lshrhw"),
    
    # ---- NEW: Social capital / networking ----
    lsclub  = paste0(wave_letter, "lsclub"),     # Active member of sporting/community/professional org (1=Yes)
    lssocal = paste0(wave_letter, "lssocal"),    # Social gathering frequency (1=every day ... 7=never)
    
    # Flexibility
    jomflex = paste0(wave_letter, "jomflex"),
    jowpfx = paste0(wave_letter, "jowpfx"),
    
    # Commute
    jbmlkr = paste0(wave_letter, "jbmlkr"),
    
    # Parental background
    fmmsch = paste0(wave_letter, "fmmsch"),
    fmfsch = paste0(wave_letter, "fmfsch"),
    
    # Geography and weights
    hhstate = paste0(wave_letter, "hhstate"),
    hhsos = paste0(wave_letter, "hhsos"),
    hhwtrps = paste0(wave_letter, "hhwtrps"),
    hhwtrp = paste0(wave_letter, "hhwtrp"),
    lnwtrp = paste0(wave_letter, "lnwtrp")    # DV: Responding person longitudinal weight
  )
  
  # ============================
  # STEP 1: CORE VARIABLES
  # ============================
  
  cat("  Step 1: Core variables...\n")
  
  # Load data (column-restricted): reading only the mapped columns plus the
  # cross-wave id cuts read time and peak memory by roughly an order of
  # magnitude versus reading the full Combined file. any_of() tolerates
  # variables absent in a given wave; the safe_numeric helpers already
  # return NA columns for those, so downstream behaviour is unchanged.
  needed_cols <- unique(c("xwaveid", unlist(v)))
  raw <- haven::read_dta(file_path, col_select = dplyr::any_of(needed_cols))
  cat("  Loaded:", format(nrow(raw), big.mark = ","), "observations,",
      ncol(raw), "of", length(needed_cols), "mapped columns\n")
  
  wave_data <- raw %>%
    mutate(
      # Identifiers
      person_id = as.character(xwaveid),
      wave = wave_num,
      year = year_val,
      
      # Demographics
      age = safe_numeric(., v$hgage),
      age_sq = age^2,
      female = ifelse(safe_numeric(., v$hgsex) == 2, 1, 0),
      # NA-propagating (same fix family as `employed`): missing mrcurr must
      # not be coded as unmarried.
      mrcurr_raw = safe_numeric(., v$mrcurr),
      married = ifelse(is.na(mrcurr_raw), NA_real_,
                       as.numeric(mrcurr_raw %in% c(1, 2))),
      
      # Employment
      esdtl_raw = safe_numeric(., v$esdtl),
      # NA-propagating: a missing esdtl must NOT be coded as non-employed
      # (the old %in% rule silently mapped NA -> 0, while fulltime/parttime
      # below correctly stay NA). Matters for extensive-margin estimates
      # (Part H) and the training event study; the published MASTER_hh
      # results do not use `employed`, so they are unaffected.
      employed = ifelse(is.na(esdtl_raw), NA_real_,
                        as.numeric(esdtl_raw %in% c(1, 2))),
      fulltime = ifelse(esdtl_raw == 1, 1, 0),
      # numeric aliases used directly in MASTER_hh.R formulas (the
      # historical W12-W24 loader carried these; never derived in-master)
      married_num  = as.numeric(married),
      fulltime_num = as.numeric(fulltime),
      parttime = ifelse(esdtl_raw == 2, 1, 0),
      
      
      fy_wages    = safe_numeric(., v$wsfei),
      fy_self_emp = safe_numeric_zero(., v$wsfes),
      annual_income = ifelse(
        is.na(fy_wages) & is.na(fy_self_emp), NA_real_,
        coalesce(fy_wages, 0) + coalesce(fy_self_emp, 0)
      ),
      
      # Wages and hours
      weekly_earnings = safe_numeric(., v$wscei),
      hours_worked = safe_numeric(., v$jbhruc),
      hours_worked_clean = case_when(
        is.na(hours_worked) ~ NA_real_,
        hours_worked <= 0 ~ NA_real_,
        hours_worked > 100 ~ NA_real_,
        TRUE ~ hours_worked
      ),
      
      hourly_wage = ifelse(
        !is.na(weekly_earnings) & !is.na(hours_worked_clean) & hours_worked_clean > 0,
        weekly_earnings / hours_worked_clean,
        NA_real_
      ),
      
      hourly_wage_clean = case_when(
        is.na(hourly_wage) ~ NA_real_,
        hourly_wage < 10 ~ NA_real_,
        hourly_wage > 500 ~ NA_real_,
        TRUE ~ hourly_wage
      ),
      
      ln_hourly_wage = ifelse(hourly_wage_clean > 0, log(hourly_wage_clean), NA_real_),
      
      # REAL WAGE (WPI adjusted)
      deflator = wave_deflator,   # carried: MASTER_hh.R re-deflates with it
      real_hourly_wage = hourly_wage_clean * wave_deflator,
      ln_hourly_wage_real = ifelse(real_hourly_wage > 0, log(real_hourly_wage), NA_real_)
    )
  
  cat("    - Core variables created\n")
  
  # ============================
  # STEP 2: ALL EDUCATION VARIABLES EXTRACTION
  # ============================
  
  cat("  Step 2: Extracting all education variables...\n")
  
  wave_data <- wave_data %>%
    mutate(
      
      childcare_cost = safe_numeric_zero(., v$ccactci),
      
      # Main education variables
      edhists_raw = safe_numeric(., v$edhists),
      edhigh1_raw = safe_numeric(., v$edhigh1),
      edsscmp_raw = safe_numeric(., v$edsscmp),
      
      # ALL completed qualifications (edqro*)
      edqrodc_raw = safe_numeric(., v$edqrodc),
      edqroms_raw = safe_numeric(., v$edqroms),
      edqrogd_raw = safe_numeric(., v$edqrogd),
      edqrogc_raw = safe_numeric(., v$edqrogc),
      edqrohd_raw = safe_numeric(., v$edqrohd),
      edqrobd_raw = safe_numeric(., v$edqrobd),
      edqroad_raw = safe_numeric(., v$edqroad),
      edqroav_raw = safe_numeric(., v$edqroav),
      edqrodp_raw = safe_numeric(., v$edqrodp),
      edqrodn_raw = safe_numeric(., v$edqrodn),
      edqroc4_raw = safe_numeric(., v$edqroc4),
      edqroc3_raw = safe_numeric(., v$edqroc3),
      edqroc2_raw = safe_numeric(., v$edqroc2),
      edqroc1_raw = safe_numeric(., v$edqroc1),
      edqrocd_raw = safe_numeric(., v$edqrocd),
      edqrosh_raw = safe_numeric(., v$edqrosh),
      edqrosl_raw = safe_numeric(., v$edqrosl),
      
      # ALL obtained qualifications (edqo*)
      edqodc_raw = safe_numeric(., v$edqodc),
      edqoms_raw = safe_numeric(., v$edqoms),
      edqogd_raw = safe_numeric(., v$edqogd),
      edqogc_raw = safe_numeric(., v$edqogc),
      edqohd_raw = safe_numeric(., v$edqohd),
      edqobd_raw = safe_numeric(., v$edqobd),
      edqoad_raw = safe_numeric(., v$edqoad),
      edqoav_raw = safe_numeric(., v$edqoav),
      edqodp_raw = safe_numeric(., v$edqodp),
      edqodn_raw = safe_numeric(., v$edqodn),
      edqoc4_raw = safe_numeric(., v$edqoc4),
      edqoc3_raw = safe_numeric(., v$edqoc3),
      edqoc2_raw = safe_numeric(., v$edqoc2),
      edqoc1_raw = safe_numeric(., v$edqoc1),
      edqocd_raw = safe_numeric(., v$edqocd),
      edqota_raw = safe_numeric(., v$edqota),
      edqotc_raw = safe_numeric(., v$edqotc),
      edqonq_raw = safe_numeric(., v$edqonq),
      edqotq_raw = safe_numeric(., v$edqotq),
      edqobc_raw = safe_numeric(., v$edqobc),
      edqosh_raw = safe_numeric(., v$edqosh),
      edqosl_raw = safe_numeric(., v$edqosl),
      
      # Nursing-specific
      ednrsdc_raw = safe_numeric(., v$ednrsdc),
      ednrsms_raw = safe_numeric(., v$ednrsms),
      ednrsgd_raw = safe_numeric(., v$ednrsgd),
      ednrsbd_raw = safe_numeric(., v$ednrsbd),
      ednrshd_raw = safe_numeric(., v$ednrshd),
      ednrsav_raw = safe_numeric(., v$ednrsav),
      ednrsdp_raw = safe_numeric(., v$ednrsdp),
      ednrsc3_raw = safe_numeric(., v$ednrsc3),
      ednrstd_raw = safe_numeric(., v$ednrstd),
      
      # Teaching-specific
      edtchdc_raw = safe_numeric(., v$edtchdc),
      edtchms_raw = safe_numeric(., v$edtchms),
      edtchgd_raw = safe_numeric(., v$edtchgd),
      edtchgc_raw = safe_numeric(., v$edtchgc),
      edtchhd_raw = safe_numeric(., v$edtchhd),
      edtchbd_raw = safe_numeric(., v$edtchbd),
      edtchts_raw = safe_numeric(., v$edtchts),
      edtchtp_raw = safe_numeric(., v$edtchtp),
      edtchdp_raw = safe_numeric(., v$edtchdp),
      edtchav_raw = safe_numeric(., v$edtchav),
      
      # Job-related training variables
      jttdays_raw  = safe_numeric(., v$jttdays),
      jtthrs_raw   = safe_numeric(., v$jtthrs),
      jttnum_raw   = safe_numeric(., v$jttnum),
      jttopot_raw  = safe_numeric(., v$jttopot),
      jttopwt_raw  = safe_numeric(., v$jttopwt),
      jttpeot_raw  = safe_numeric(., v$jttpeot),
      jttpewt_raw  = safe_numeric(., v$jttpewt),
      
    )
  
  cat("    - All education variables extracted\n")
  
  
  
  # ============================
  # STEP 3: COMPREHENSIVE EDUCATION YEARS AND PATHWAYS
  # ============================
  
  cat("  Step 3: Creating detailed education variables...\n")
  
  wave_data <- wave_data %>%
    mutate(
      # Extract school years from history for granularity
      school_years_from_hist = case_when(
        edhists_raw == 1 ~ 12,  # Year 12
        edhists_raw == 2 ~ 11,  # Year 11
        edhists_raw == 3 ~ 10,  # Year 10
        edhists_raw == 4 ~ 9,   # Year 9
        edhists_raw == 5 ~ 8,   # Year 8
        edhists_raw == 6 ~ 7,   # Year 7
        edhists_raw == 7 ~ 6,   # Finished primary
        edhists_raw == 8 ~ 5,   # Attended but didn't finish primary
        TRUE ~ NA_real_
      ),
      
      # Create detailed education years with maximum granularity
      educ_years = case_when(
        # === DOCTORATE LEVEL (19-20 years) ===
        edqrodc_raw == 1 | edqodc_raw == 1 ~ 20,
        ednrsdc_raw == 1 ~ 19.5,  # Nursing doctorate
        edtchdc_raw == 1 ~ 19.5,  # Teaching doctorate
        
        # === MASTERS LEVEL (18-18.5 years) ===
        edqroms_raw == 1 | edqoms_raw == 1 ~ 18,
        ednrsms_raw == 1 ~ 18,    # Nursing masters
        edtchms_raw == 1 ~ 18,    # Teaching masters
        
        # === GRADUATE DIPLOMA/CERTIFICATE (17 years) ===
        edqrogd_raw == 1 | edqogd_raw == 1 ~ 17,
        ednrsgd_raw == 1 ~ 17,    # Nursing grad diploma
        edtchgd_raw == 1 ~ 17,    # Teaching grad diploma
        edqrogc_raw == 1 | edqogc_raw == 1 ~ 17,
        edtchgc_raw == 1 ~ 17,    # Teaching grad certificate
        
        # === BACHELOR LEVEL (16-16.5 years) ===
        edqrohd_raw == 1 | edqohd_raw == 1 ~ 16.5,  # Honours
        ednrshd_raw == 1 ~ 16.5,  # Nursing honours
        edtchhd_raw == 1 ~ 16.5,  # Teaching honours
        edqrobd_raw == 1 | edqobd_raw == 1 ~ 16,    # Ordinary bachelor
        ednrsbd_raw == 1 ~ 16,    # Nursing bachelor
        edtchbd_raw == 1 ~ 16,    # Teaching bachelor
        
        # === ADVANCED DIPLOMA LEVEL (15 years) ===
        edqroav_raw == 1 | edqoav_raw == 1 ~ 15,
        ednrsav_raw == 1 ~ 15,    # Nursing advanced diploma
        edtchav_raw == 1 ~ 15,    # Teaching advanced diploma
        edtchts_raw == 1 ~ 15,    # Secondary teaching cert
        
        # === DIPLOMA LEVEL (14-14.5 years) ===
        edqrodp_raw == 1 | edqodp_raw == 1 ~ 14.5,
        ednrsdp_raw == 1 ~ 14.5,  # Nursing diploma
        edtchdp_raw == 1 ~ 14.5,  # Teaching diploma
        edqroad_raw == 1 | edqoad_raw == 1 ~ 14.5,  # Associate degree
        edqrodn_raw == 1 | edqodn_raw == 1 ~ 14.5,  # Diploma nursing/unknown
        ednrstd_raw == 1 ~ 14,    # Registered Nurse diploma
        edtchtp_raw == 1 ~ 14,    # Primary teaching cert
        
        # === CERTIFICATE IV (13.5 years) ===
        edqroc4_raw == 1 | edqoc4_raw == 1 ~ 13.5,
        
        # === CERTIFICATE III/TRADE (13 years) ===
        edqroc3_raw == 1 | edqoc3_raw == 1 ~ 13,
        ednrsc3_raw == 1 ~ 13,    # Nursing cert III
        edqota_raw == 1 ~ 13,     # Trade certificate
        edqotc_raw == 1 ~ 13,     # Technical certificate
        
        # === CERTIFICATE II (12.5 years) ===
        edqroc2_raw == 1 | edqoc2_raw == 1 ~ 12.5,
        edqrocd_raw == 1 | edqocd_raw == 1 ~ 12.5,  # Cert unknown level
        
        # === CERTIFICATE I (12 years) ===
        edqroc1_raw == 1 | edqoc1_raw == 1 ~ 12,
        
        # === SECONDARY SCHOOL COMPLETION ===
        edqrosh_raw == 1 | edqosh_raw == 1 ~ 12,  # Year 12
        
        # === FROM edhigh1 (if not captured above) ===
        edhigh1_raw == 1 ~ 18,    # Postgrad
        edhigh1_raw == 2 ~ 17,    # Grad diploma/cert
        edhigh1_raw == 3 ~ 16,    # Bachelor
        edhigh1_raw == 4 ~ 15,    # Advanced diploma
        edhigh1_raw == 5 ~ 13,    # Cert III/IV
        edhigh1_raw == 8 ~ 12,    # Year 12
        
        # Use detailed school years when edhigh1 indicates below Year 12
        edhigh1_raw == 9 & !is.na(school_years_from_hist) ~ school_years_from_hist,
        edhigh1_raw == 9 & is.na(school_years_from_hist) ~ 10,  # Default to Year 10
        edhigh1_raw == 10 & !is.na(school_years_from_hist) ~ school_years_from_hist,
        
        # === FROM edsscmp (CPQ school completion) ===
        edsscmp_raw == 1 ~ 12,  # Year 12
        edsscmp_raw == 2 ~ 11,  # Year 11
        edsscmp_raw == 3 ~ 10,  # Year 10
        edsscmp_raw == 4 ~ 9,   # Year 9
        
        # === FROM edqrosl/edqosl (lower secondary) ===
        edqrosl_raw == 1 | edqosl_raw == 1 ~ 10,
        
        # === FALLBACK to school history ===
        !is.na(school_years_from_hist) ~ school_years_from_hist,
        
        # === DEFAULT ===
        TRUE ~ NA_real_
      ),
      
      # Create field indicators for enriched pathways
      has_nursing_any = coalesce(
        (ednrsdc_raw == 1 | ednrsms_raw == 1 | ednrsgd_raw == 1 | 
           ednrsbd_raw == 1 | ednrshd_raw == 1 | ednrsav_raw == 1 |
           ednrsdp_raw == 1 | ednrsc3_raw == 1 | ednrstd_raw == 1 |
           edqrodn_raw == 1 | edqodn_raw == 1 | edqonq_raw == 1), FALSE),
      
      has_teaching_any = coalesce(
        (edtchdc_raw == 1 | edtchms_raw == 1 | edtchgd_raw == 1 |
           edtchgc_raw == 1 | edtchhd_raw == 1 | edtchbd_raw == 1 |
           edtchts_raw == 1 | edtchtp_raw == 1 | edtchdp_raw == 1 |
           edtchav_raw == 1 | edqotq_raw == 1), FALSE),
      
      has_business_any = coalesce(edqobc_raw == 1, FALSE),
      has_trade_any = coalesce(edqota_raw == 1, FALSE),
      has_technical_any = coalesce(edqotc_raw == 1, FALSE),
      
      # Count qualifications for complexity
      qual_count = rowSums(across(c(starts_with("edqro"), starts_with("edqo")) & 
                                    ends_with("_raw"), ~ .x == 1), na.rm = TRUE),
      
      # Create detailed education pathway
      education_pathway = case_when(
        # === DOCTORATE PATHWAYS ===
        educ_years >= 19 & has_nursing_any ~ "Doctorate_Nursing",
        educ_years >= 19 & has_teaching_any ~ "Doctorate_Teaching",
        educ_years >= 19 & has_business_any ~ "Doctorate_Business",
        educ_years >= 19 ~ "Doctorate_General",
        
        # === MASTERS PATHWAYS ===
        educ_years >= 18 & has_nursing_any ~ "Masters_Nursing",
        educ_years >= 18 & has_teaching_any ~ "Masters_Teaching",
        educ_years >= 18 & has_business_any ~ "Masters_Business",
        educ_years >= 18 ~ "Masters_General",
        
        # === GRADUATE DIPLOMA/CERT PATHWAYS ===
        educ_years >= 17 & has_nursing_any ~ "GradDipCert_Nursing",
        educ_years >= 17 & has_teaching_any ~ "GradDipCert_Teaching",
        educ_years >= 17 & has_business_any ~ "GradDipCert_Business",
        educ_years >= 17 ~ "GradDipCert_General",
        
        # === BACHELOR PATHWAYS ===
        educ_years >= 16.5 & has_nursing_any ~ "Bachelor_Honours_Nursing",
        educ_years >= 16.5 & has_teaching_any ~ "Bachelor_Honours_Teaching",
        educ_years >= 16.5 & has_business_any ~ "Bachelor_Honours_Business",
        educ_years >= 16.5 ~ "Bachelor_Honours_General",
        
        educ_years >= 16 & has_nursing_any ~ "Bachelor_Nursing",
        educ_years >= 16 & has_teaching_any ~ "Bachelor_Teaching",
        educ_years >= 16 & has_business_any ~ "Bachelor_Business",
        educ_years >= 16 ~ "Bachelor_General",
        
        # === SPECIAL PROFESSIONAL CERTIFICATES ===
        ednrstd_raw == 1 ~ "Registered_Nurse",
        edtchts_raw == 1 ~ "Secondary_Teaching_Cert",
        edtchtp_raw == 1 ~ "Primary_Teaching_Cert",
        
        # === ADVANCED DIPLOMA PATHWAYS ===
        educ_years >= 15 & has_nursing_any ~ "AdvDiploma_Nursing",
        educ_years >= 15 & has_teaching_any ~ "AdvDiploma_Teaching",
        educ_years >= 15 & has_trade_any ~ "AdvDiploma_Trade",
        educ_years >= 15 & has_technical_any ~ "AdvDiploma_Technical",
        educ_years >= 15 ~ "Advanced_Diploma",
        
        # === DIPLOMA PATHWAYS ===
        educ_years >= 14 & has_nursing_any ~ "Diploma_Nursing",
        educ_years >= 14 & has_teaching_any ~ "Diploma_Teaching",
        educ_years >= 14 & has_trade_any ~ "Diploma_Trade",
        educ_years >= 14 & has_technical_any ~ "Diploma_Technical",
        educ_years >= 14 & has_business_any ~ "Diploma_Business",
        educ_years >= 14 ~ "Diploma",
        
        # === CERTIFICATE IV PATHWAYS ===
        educ_years >= 13.5 & has_trade_any ~ "Cert_IV_Trade",
        educ_years >= 13.5 & has_technical_any ~ "Cert_IV_Technical",
        educ_years >= 13.5 & has_nursing_any ~ "Cert_IV_Nursing",
        educ_years >= 13.5 ~ "Certificate_IV",
        
        # === CERTIFICATE III/TRADE PATHWAYS ===
        educ_years >= 13 & has_trade_any ~ "Cert_III_Trade",
        educ_years >= 13 & has_technical_any ~ "Cert_III_Technical",
        educ_years >= 13 & has_nursing_any ~ "Cert_III_Nursing",
        educ_years >= 13 ~ "Certificate_III",
        
        # === CERTIFICATE II ===
        educ_years >= 12.5 ~ "Certificate_II",
        
        # === CERTIFICATE I ===
        educ_years > 12 & educ_years < 12.5 ~ "Certificate_I",
        
        # === SECONDARY SCHOOL ===
        educ_years == 12 ~ "Year_12",
        educ_years == 11 ~ "Year_11",
        educ_years == 10 ~ "Year_10",
        educ_years == 9 ~ "Year_9",
        educ_years == 8 ~ "Year_8",
        educ_years == 7 ~ "Year_7",
        educ_years == 6 ~ "Primary_Complete",
        educ_years == 5 ~ "Primary_Incomplete",
        
        # === NO EDUCATION DATA ===
        is.na(educ_years) ~ "No_Education_Data",
        TRUE ~ "Other"
      ),
      
      # Education squared for non-linear effects
      educ_sq = educ_years^2,
      
      # Create education level categories
      education_level = case_when(
        educ_years >= 19 ~ "Doctorate",
        educ_years >= 18 ~ "Masters",
        educ_years >= 17 ~ "Graduate_Diploma_Certificate",
        educ_years >= 16.5 ~ "Bachelor_Honours",
        educ_years >= 16 ~ "Bachelor",
        educ_years >= 15 ~ "Advanced_Diploma",
        educ_years >= 14 ~ "Diploma",
        educ_years >= 13.5 ~ "Certificate_IV",
        educ_years >= 13 ~ "Certificate_III",
        educ_years >= 12.5 ~ "Certificate_II",
        educ_years > 12 ~ "Certificate_I",
        educ_years == 12 ~ "Year_12",
        educ_years == 11 ~ "Year_11",
        educ_years == 10 ~ "Year_10",
        educ_years < 10 & !is.na(educ_years) ~ "Below_Year_10",
        TRUE ~ "Missing"
      ),
      
      # Binary indicators for analysis
      has_degree = (educ_years >= 16),
      has_postsecondary = (educ_years > 12),
      has_diploma_plus = (educ_years >= 14),
      has_postgrad = (educ_years >= 17),
      
      # Data quality indicator
      education_data_source = case_when(
        edqrodc_raw == 1 | edqroms_raw == 1 | edqrogd_raw == 1 | edqrogc_raw == 1 ~ "Specific_Postgrad",
        edqrohd_raw == 1 | edqrobd_raw == 1 ~ "Specific_Bachelor",
        edqroav_raw == 1 | edqrodp_raw == 1 | edqroad_raw == 1 ~ "Specific_Diploma",
        edqroc1_raw == 1 | edqroc2_raw == 1 | edqroc3_raw == 1 | edqroc4_raw == 1 ~ "Specific_Certificate",
        !is.na(edhigh1_raw) ~ "General_Highest",
        !is.na(edhists_raw) ~ "School_History",
        !is.na(edsscmp_raw) ~ "CPQ_School",
        TRUE ~ "No_Data"
      )
    )
  
  # Education diagnostics
  cat("    - Education variables created\n")
  cat("    - Education coverage:", 
      round(100 * sum(!is.na(wave_data$educ_years)) / nrow(wave_data), 1), "%\n")
  cat("    - Mean education years:", round(mean(wave_data$educ_years, na.rm = TRUE), 2), "\n")
  cat("    - Nursing professionals:", sum(wave_data$has_nursing_any, na.rm = TRUE), "\n")
  cat("    - Teaching professionals:", sum(wave_data$has_teaching_any, na.rm = TRUE), "\n")
  cat("    - With degrees (16+):", sum(wave_data$has_degree, na.rm = TRUE), "\n")
  cat("    - With postgrad (17+):", sum(wave_data$has_postgrad, na.rm = TRUE), "\n")
  
  # ============================
  # STEP 4: EXPERIENCE AND JOB CHARACTERISTICS (CORRECTED)
  # ============================
  
  cat("  Step 4: Experience and job characteristics...\n")
  
  wave_data <- wave_data %>%
    mutate(
      # Experience
      # Use HILDA's cumulative employment history variable (ehtjb): total years in
      # paid work since leaving full-time education, constructed from the employment
      # history interview and updated each wave from the labour market activity
      # calendar. Directly captures actual experience including work-study overlap
      # and career interruptions -- superior to the Mincer potential experience proxy.
      # Fallback to potential experience (age - educ_years - 5) only when ehtjb is
      # missing. The Mincer formula ceiling is NOT applied: ehtjb legitimately exceeds
      # it for educated workers who worked during study (38-55% of post-secondary
      # respondents). Only a biological ceiling (age - 13, Australia's minimum working
      # age) is applied, removing 10 person-wave observations across 3 individuals
      # with implied labour market entry below age 10 (< 0.02% of sample).
      actual_experience    = safe_numeric_zero(., v$ehtjb),
      # NA-propagating (fixed): with the old na.rm = TRUE, a missing age or
      # educ_years collapsed to potential_experience = 0, and if ehtjb was
      # also missing, experience_years became a spurious 0 instead of NA.
      # Now missing inputs yield NA and such rows drop from experience-based
      # samples. NOTE: this can shift sample counts relative to the panel
      # used for the current manuscript; re-verify the FINAL VERIFICATION
      # block (esp. paper Sec. 3 counts) after rebuilding.
      potential_experience = pmax(0, age - educ_years - 5),
      experience_years     = pmin(
        coalesce(actual_experience, potential_experience),
        pmax(0, age - 13)            # biological ceiling: implied start age >= 13
      ),
      experience_years     = ifelse(experience_years < 0, 0, experience_years),
      experience_sq        = experience_years^2,
      
      # Total number of training days across all types
      total_training_days = coalesce(jttdays_raw, 0),
      
      # Average hours per course day (already available as jtthrs_raw)
      avg_training_hours = jtthrs_raw,
      
      # Total number of courses attended in last 12 months
      total_courses = coalesce(jttnum_raw, 0),
      
      # Indicators for training done in own time or paid work time
      # (pooled over place of employment + some other place)
      # NOTE: "mentioned" coding made robust with == 1, which is correct
      # under both 0/1 and 1/2 (yes/no) HILDA codings. The previous "> 0"
      # rule would count 2 = "not mentioned" as training if the items are
      # 1/2-coded; verify with table(jttpewt_raw) on the first run.
      training_own_time = ifelse(
        coalesce(jttopot_raw, 0) == 1 | coalesce(jttpeot_raw, 0) == 1, 1, 0),
      training_paid_time = ifelse(
        coalesce(jttopwt_raw, 0) == 1 | coalesce(jttpewt_raw, 0) == 1, 1, 0),
      
      # PLACE-OF-EMPLOYMENT-ONLY variants (for training_event_study.R):
      # jttpewt = at place of employment, during paid work time
      # jttpeot = at place of employment, in own time
      jttpewt_any = ifelse(coalesce(jttpewt_raw, 0) == 1, 1, 0),
      jttpeot_any = ifelse(coalesce(jttpeot_raw, 0) == 1, 1, 0),
      
      # Tenure
      tenure_employer = safe_numeric_zero(., v$jbempt),
      tenure_occupation = safe_numeric_zero(., v$jbocct),
      
      # Employment type 
      casual = ifelse(safe_numeric(., v$jbcasab) == 1, 1, 0),
      permanent = ifelse(safe_numeric(., v$jbcasab) == 2, 1, 0),
      
      # Leave
      has_holiday_leave = ifelse(safe_numeric(., v$jbmhl) == 1, 1, 0),
      has_sick_leave = ifelse(safe_numeric(., v$jbmsl) == 1, 1, 0),
      leave_none = ifelse(has_holiday_leave == 0 & has_sick_leave == 0, 1, 0),
      
      # Schedule
      jbmsch_raw = safe_numeric(., v$jbmsch),
      has_standard_hours = ifelse(jbmsch_raw == 1, 1, 0),
      has_shift_work = ifelse(jbmsch_raw %in% c(2,3,4,5), 1, 0),
      has_irregular_hours = ifelse(jbmsch_raw %in% c(6,7), 1, 0),
      
      # Hours preferences
      jbhrcpr_raw = safe_numeric(., v$jbhrcpr),
      underemployed = ifelse(jbhrcpr_raw == 3, 1, 0),
      overemployed = ifelse(jbhrcpr_raw == 1, 1, 0),
      
      # Union and sector
      union_member = ifelse(safe_numeric(., v$jbmtuea) == 1, 1, 0),
      public_sector = ifelse(safe_numeric(., v$jbmmply) %in% c(2,5), 1, 0),
      public_sector = ifelse(is.na(public_sector), 0, public_sector),
      
      # Part-time reason
      jbptrea_raw = safe_numeric(., v$jbptrea),
      pt_childcare = ifelse(parttime == 1 & jbptrea_raw == 2, 1, 0),
      pt_involuntary = ifelse(parttime == 1 & jbptrea_raw == 6, 1, 0)
    )
  
  # FIRM SIZE - SEPARATE MUTATE FOR CLARITY
  wave_data <- wave_data %>%
    mutate(
      # First extract the raw firm size variable
      firm_size_raw = if(wave_num <= 21) {
        safe_numeric(., v$jbmemsz)
      } else {
        safe_numeric(., v$jbmems2)
      }
    ) %>%
    mutate(
      # Then map to categories
      firm_size = case_when(
        # Old waves <=21
        wave_num <= 21 & firm_size_raw == 1 ~ "Small",
        wave_num <= 21 & firm_size_raw %in% 2:3 ~ "Medium",
        wave_num <= 21 & firm_size_raw %in% 4:9 ~ "Large",
        
        # New waves >=22
        wave_num >= 22 & firm_size_raw %in% 1:4 ~ "Small",
        wave_num >= 22 & firm_size_raw %in% 5:6 ~ "Medium",
        wave_num >= 22 & firm_size_raw %in% 7:10 ~ "Large",
        wave_num >= 22 & firm_size_raw == 11 ~ "Small",
        wave_num >= 22 & firm_size_raw %in% 12:13 ~ "Large",
        
        TRUE ~ NA_character_
      ),
      
      # Commute
      job_distance_cat = safe_numeric(., v$jbmlkr),
      job_distance_km = case_when(
        job_distance_cat == 0 ~ 0,
        job_distance_cat == 1 ~ 2.5,
        job_distance_cat == 2 ~ 7,
        job_distance_cat == 3 ~ 12,
        job_distance_cat == 4 ~ 17,
        job_distance_cat == 5 ~ 24.5,
        job_distance_cat == 6 ~ 39.5,
        job_distance_cat == 7 ~ 74.5,
        job_distance_cat == 8 ~ 124.5,
        job_distance_cat == 9 ~ 150,
        TRUE ~ NA_real_
      )
    )
  
  # Add diagnostic output for firm size
  cat("    - Firm size distribution:\n")
  firm_dist <- table(wave_data$firm_size, useNA = "ifany")
  print(firm_dist)
  
  # ============================
  # STEP 5: OCCUPATION AND INDUSTRY
  # ============================
  
  cat("  Step 5: Occupation and industry...\n")
  
  wave_data <- wave_data %>%
    mutate(
      # Occupation
      occupation_major_code = safe_numeric(., v$jbmo61),
      occupation_major = case_when(
        occupation_major_code == 1 ~ "Managers",
        occupation_major_code == 2 ~ "Professionals",
        occupation_major_code == 3 ~ "Technicians_Trades",
        occupation_major_code == 4 ~ "Community_Service",
        occupation_major_code == 5 ~ "Clerical_Admin",
        occupation_major_code == 6 ~ "Sales",
        occupation_major_code == 7 ~ "Machinery_Operators",
        occupation_major_code == 8 ~ "Labourers",
        TRUE ~ NA_character_
      ),
      
      # Industry
      industry_major_code = safe_numeric(., v$jbmi61),
      industry_major = case_when(
        industry_major_code == 1 ~ "Agriculture",
        industry_major_code == 2 ~ "Mining",
        industry_major_code == 3 ~ "Manufacturing",
        industry_major_code == 4 ~ "Utilities",
        industry_major_code == 5 ~ "Construction",
        industry_major_code == 6 ~ "Wholesale",
        industry_major_code == 7 ~ "Retail",
        industry_major_code == 8 ~ "Hospitality",
        industry_major_code == 9 ~ "Transport",
        industry_major_code == 10 ~ "Information",
        industry_major_code == 11 ~ "Finance",
        industry_major_code == 12 ~ "Real_Estate",
        industry_major_code == 13 ~ "Professional",
        industry_major_code == 14 ~ "Administrative",
        industry_major_code == 15 ~ "Public_Admin",
        industry_major_code == 16 ~ "Education",
        industry_major_code == 17 ~ "Health",
        industry_major_code == 18 ~ "Arts",
        industry_major_code == 19 ~ "Other_Services",
        TRUE ~ NA_character_
      ),
      
      # Field of education
      occupation_field_code = safe_numeric(., v$edpsfdn)
    )
  
  # ============================
  # STEP 6: CHILDREN, FERTILITY HISTORY, AND FLEXIBILITY
  # ============================
  
  cat("  Step 6: Children, fertility history, flexibility, and care...\n")
  
  wave_data <- wave_data %>%
    mutate(
      # -----------------------------------------------------------
      # 6A. HOUSEHOLD CHILDREN COUNTS (original)
      # -----------------------------------------------------------
      num_children_0to4 = safe_numeric_zero(., v$hhd0_4),
      num_children_5to9 = safe_numeric_zero(., v$hhd5_9),
      num_children_10to14 = safe_numeric_zero(., v$hhd1014),
      num_children_15to24 = safe_numeric_zero(., v$hhd1524),
      
      num_dependent_children = rowSums(
        cbind(num_children_0to4, num_children_5to9, 
              num_children_10to14, num_children_15to24),
        na.rm = TRUE
      ),
      
      num_children_under15 = rowSums(
        cbind(num_children_0to4, num_children_5to9, num_children_10to14),
        na.rm = TRUE
      ),
      
      # Non-resident children (original)
      ncchk_raw = safe_numeric(., v$ncchk),
      ncany_raw = safe_numeric(., v$ncany),
      has_nonresident_children_under18 = (ncchk_raw %in% 1:2),
      has_any_nonresident_children = (ncany_raw == 1),
      
      # -----------------------------------------------------------
      # 6B. FERTILITY HISTORY (tc* prefix — NEW)
      # -----------------------------------------------------------
      
      # tchad: "Total children ever had" — History variable, absorbing
      #   Available most waves as History or PQ/NPQ
      tchad_raw  = safe_numeric_zero(., v$tchad),
      tchadn_raw = safe_numeric_zero(., v$tchadn),
      tchave_raw = safe_numeric_zero(., v$tchave),
      tchch_raw  = safe_numeric(., v$tchch),       # Has ever had children check (fertility module waves)
      
      # Best available total children ever had (prefer tchad, fallback tchadn, tchave)
      total_children_ever_had = coalesce(tchad_raw, tchadn_raw, tchave_raw),
      
      # Own children counts — DV, available ALL waves (most reliable)
      num_own_resident_children    = safe_numeric_zero(., v$tcr),
      num_own_nonresident_children = safe_numeric_zero(., v$tcnr),
      total_own_children = rowSums(
        cbind(num_own_resident_children, num_own_nonresident_children),
        na.rm = TRUE
      ),
      
      # Age of youngest own child (DV, all waves)
      age_youngest_own_child = safe_numeric(., v$tcyng),
      
      # Children since died
      total_children_died = safe_numeric_zero(., v$tcdied),
      
      # Resident own children + step/foster/grand by age group (DV, all waves)
      own_res_children_0to4   = safe_numeric_zero(., v$tcr04),
      own_res_children_5to14  = safe_numeric_zero(., v$tcr514),
      own_res_children_15to24 = safe_numeric_zero(., v$tcr1524),
      own_res_children_25plus = safe_numeric_zero(., v$tcr25),
      
      # Non-resident own children by age group (DV, all waves)
      own_nonres_children_0to4   = safe_numeric_zero(., v$tcnr04),
      own_nonres_children_5to14  = safe_numeric_zero(., v$tcnr514),
      own_nonres_children_15to24 = safe_numeric_zero(., v$tcnr1524),
      own_nonres_children_25plus = safe_numeric_zero(., v$tcnr25),
      
      # Step/foster/grandchildren with no parent in HH
      num_step_foster_grand = safe_numeric_zero(., v$rcstepn),
      
      # Child ages (for estimating age at first birth)
      age_resident_child_1    = safe_numeric(., v$rcage1),
      age_resident_child_2    = safe_numeric(., v$rcage2),
      age_resident_child_3    = safe_numeric(., v$rcage3),
      age_nonresident_child_1 = safe_numeric(., v$ncage1),
      age_nonresident_child_2 = safe_numeric(., v$ncage2),
      
      # -----------------------------------------------------------
      # 6C. ABSORBING EVER-PARENT INDICATOR (NEW)
      # -----------------------------------------------------------
      # Three sources combined for maximum coverage:
      #   1. tchad > 0          (total children ever had — History, most waves)
      #   2. tcr + tcnr > 0     (own resident + non-resident — DV, all waves)
      #   3. hhd children > 0   (household children — original, all waves)
      #   4. ncany == 1         (any non-resident children — original, all waves)
      # Any one sufficient to classify as ever-parent.
      
      ever_parent = case_when(
        total_children_ever_had > 0 ~ 1L,
        total_own_children > 0 ~ 1L,
        num_children_under15 > 0 ~ 1L,
        has_any_nonresident_children == TRUE ~ 1L,
        # If tchad is NA but all other sources are 0, assume childless
        is.na(total_children_ever_had) & total_own_children == 0 & 
          num_children_under15 == 0 & (is.na(has_any_nonresident_children) | has_any_nonresident_children == FALSE) ~ 0L,
        # If all available and zero
        total_children_ever_had == 0 & total_own_children == 0 ~ 0L,
        TRUE ~ NA_integer_
      ),
      
      # Gender-specific ever-parent
      ever_mother = ifelse(female == 1, ever_parent, NA_integer_),
      ever_father = ifelse(female == 0, ever_parent, NA_integer_),
      
      # -----------------------------------------------------------
      # 6D. THREE-WAY MOTHERHOOD CLASSIFICATION (NEW)
      # -----------------------------------------------------------
      # never_mother:  ever_parent == 0 (never had children by any measure)
      # active_mother: ever_parent == 1 AND has children under 15 in HH now
      # post_mother:   ever_parent == 1 AND NO children under 15 in HH now
      
      motherhood_status = case_when(
        female == 0 ~ NA_character_,
        ever_parent == 1 & num_children_under15 > 0 ~ "active_mother",
        ever_parent == 1 & num_children_under15 == 0 ~ "post_mother",
        ever_parent == 0 ~ "never_mother",
        TRUE ~ NA_character_
      ),
      
      # Fatherhood equivalent
      fatherhood_status = case_when(
        female == 1 ~ NA_character_,
        ever_parent == 1 & num_children_under15 > 0 ~ "active_father",
        ever_parent == 1 & num_children_under15 == 0 ~ "post_father",
        ever_parent == 0 ~ "never_father",
        TRUE ~ NA_character_
      ),
      
      # Approximate age at first birth (from oldest known child)
      age_oldest_child = pmax(age_resident_child_1, age_resident_child_2,
                              age_resident_child_3, age_nonresident_child_1,
                              age_nonresident_child_2, na.rm = TRUE),
      approx_age_at_first_birth = ifelse(!is.na(age_oldest_child) & ever_parent == 1,
                                         age - age_oldest_child, NA_real_),
      
      # -----------------------------------------------------------
      # 6E. ORIGINAL CHILD INDICATORS (preserved for backward compatibility)
      # -----------------------------------------------------------
      has_young_children = ifelse(num_children_0to4 > 0, 1, 0),
      active_caregiver = ifelse(num_dependent_children > 0, 1, 0),
      is_mother = ifelse(female == 1 & num_dependent_children > 0, 1, 0),
      is_father = ifelse(female == 0 & num_dependent_children > 0, 1, 0),
      
      # -----------------------------------------------------------
      # 6F. PARENTING ATTITUDES (NEW)
      # -----------------------------------------------------------
      parenting_resp      = safe_numeric(., v$paresp),    # Has responsibilities, kids <= 17
      parenting_hard      = safe_numeric(., v$pahard),    # "Harder than I thought"
      parenting_tired     = safe_numeric(., v$patird),    # "Often tired/exhausted"
      parenting_trapped   = safe_numeric(., v$patrap),    # "Feel trapped"
      parenting_fairshare = safe_numeric(., v$pashare),   # "Do fair share"
      
      # -----------------------------------------------------------
      # 6G. PREGNANCY (fertility module waves: 15, 19, 23 in our panel)
      # Will be NA for non-fertility-module waves — this is expected
      # -----------------------------------------------------------
      currently_pregnant       = safe_numeric(., v$ftcpg),
      partner_pregnant         = safe_numeric(., v$ftcppg),
      self_or_partner_pregnant = safe_numeric(., v$ftipg),
      expected_birth_month     = safe_numeric(., v$ftcpgm),
      expected_birth_year      = safe_numeric(., v$ftcpgy),
      
      # -----------------------------------------------------------
      # 6H. LIFE EVENTS (asked annually — all waves) (NEW)
      # -----------------------------------------------------------
      event_birth     = safe_numeric(., v$leprg),   # Pregnancy/birth of new child
      event_married   = safe_numeric(., v$lemar),   # Got married
      event_separated = safe_numeric(., v$lesep),   # Separated from spouse
      event_job_change= safe_numeric(., v$lejob),   # Changed jobs
      event_fired     = safe_numeric(., v$lejlf),   # Fired/made redundant
      event_moved     = safe_numeric(., v$lemvd),   # Moved residence
      
      # Binary indicators for life events (HILDA coding: 1=No, 2=Yes)
      had_birth_this_wave     = ifelse(event_birth == 2, 1L, 0L),
      got_married_this_wave   = ifelse(event_married == 2, 1L, 0L),
      separated_this_wave     = ifelse(event_separated == 2, 1L, 0L),
      changed_job_this_wave   = ifelse(event_job_change == 2, 1L, 0L),
      
      # -----------------------------------------------------------
      # 6I. OTHER PARENT / PARTNER INFO (NEW)
      # -----------------------------------------------------------
      other_parent_employment = safe_numeric(., v$ncesop),
      living_with_partner_last = safe_numeric(., v$prlpli),
      # Household and partner linkage
      household_id = if (v$hhrhid %in% names(.)) as.character(.[[v$hhrhid]]) else NA_character_,
      partner_id   = if (v$hhpxid %in% names(.)) as.character(.[[v$hhpxid]]) else NA_character_,
      
      # -----------------------------------------------------------
      # 6J. JOB SATISFACTION & WORKPLACE QUALITY (NEW)
      # -----------------------------------------------------------
      job_satisfaction    = safe_numeric(., v$josat),     # Overall work enjoyment/satisfaction
      pay_fairness        = safe_numeric(., v$jompf),     # "Paid fairly"
      job_initiative      = safe_numeric(., v$jomini),    # "Job requires initiative"
      job_security        = safe_numeric(., v$jomsf),     # "Secure future in my job"
      job_worry_future    = safe_numeric(., v$jomwf),     # "Worry about future of my job"
      job_autonomy        = safe_numeric(., v$jomls_say), # "Have a lot of say"
      looking_for_new_job = safe_numeric(., v$jbnewjs),   # Looked for new job last 4 weeks
      
      # -----------------------------------------------------------
      # 6K. SOCIAL CAPITAL / NETWORKING (NEW — H4 mechanism)
      # -----------------------------------------------------------
      # lsclub: "Are you an active member of a sporting, hobby or community-based
      #          club or association?" 1=Yes, 2=No → recode to 0/1
      club_member_raw = safe_numeric(., v$lsclub),
      club_member = case_when(
        club_member_raw == 1 ~ 1L,
        club_member_raw == 2 ~ 0L,
        TRUE ~ NA_integer_
      ),
      # lssocal: "How often do you get together socially with friends/relatives
      #           not living with you?" 1=Every day, 2=Several times/week, ..., 7=Never
      #           REVERSE so higher = more social (1=never → 7=every day)
      social_freq_raw = safe_numeric(., v$lssocal),
      social_freq = case_when(
        social_freq_raw >= 1 & social_freq_raw <= 7 ~ 8L - as.integer(social_freq_raw),
        TRUE ~ NA_integer_
      ),
      social_freq_01 = ifelse(!is.na(social_freq), (social_freq - 1) / 6, NA_real_),
      high_social = ifelse(!is.na(social_freq), as.integer(social_freq >= 5), NA_integer_),
      
      # -----------------------------------------------------------
      # 6L. CARE AND FLEXIBILITY (original)
      # -----------------------------------------------------------
      
      # Care intensity (housework hours)
      housework_hours = safe_numeric_zero(., v$lshrhw),
      care_intensity = case_when(
        is.na(housework_hours) ~ 0,
        housework_hours < 0 ~ 0,
        housework_hours > 100 ~ 100,
        TRUE ~ housework_hours
      ),
      
      # Flexibility (treatment)
      jomflex_raw = safe_numeric(., v$jomflex),
      jowpfx_raw = safe_numeric(., v$jowpfx),
      has_flex_perception = ifelse(jomflex_raw >= 5, 1, 0),
      has_flex_entitlement = ifelse(jowpfx_raw == 1, 1, 0),
      flexibility_any = ifelse(has_flex_perception == 1 | has_flex_entitlement == 1, 1, 0),
      treatment = flexibility_any,
      treatment = ifelse(is.na(treatment), 0, treatment)
    )
  
  # Fertility diagnostics
  cat("    - Fertility history variables extracted\n")
  cat("    - tchad coverage:", 
      round(100 * sum(!is.na(wave_data$total_children_ever_had)) / nrow(wave_data), 1), "%\n")
  cat("    - tcr+tcnr coverage:", 
      round(100 * sum(!is.na(wave_data$num_own_resident_children)) / nrow(wave_data), 1), "%\n")
  cat("    - ever_parent == 1:", sum(wave_data$ever_parent == 1, na.rm = TRUE), 
      "(", round(100 * mean(wave_data$ever_parent == 1, na.rm = TRUE), 1), "%)\n")
  cat("    - ever_parent == 0:", sum(wave_data$ever_parent == 0, na.rm = TRUE), "\n")
  cat("    - ever_parent NA:", sum(is.na(wave_data$ever_parent)), "\n")
  
  # ============================
  # STEP 7: PRECARIOUSNESS AND FINAL VARIABLES
  # ============================
  
  cat("  Step 7: Final variables...\n")
  
  wave_data <- wave_data %>%
    mutate(
      # Precariousness score
      precarious_score = rowSums(
        cbind(casual, 1-permanent, leave_none, has_irregular_hours, underemployed),
        na.rm = TRUE
      ),
      
      precarious_category = case_when(
        precarious_score == 0 ~ "Secure",
        precarious_score == 1 ~ "Mostly_Secure",
        precarious_score == 2 ~ "Somewhat_Precarious",
        precarious_score >= 3 ~ "Highly_Precarious",
        TRUE ~ NA_character_
      ),
      
      # Parental education
      mother_educ_years = case_when(
        safe_numeric(., v$fmmsch) == 1 ~ 0,
        safe_numeric(., v$fmmsch) == 2 ~ 6,
        safe_numeric(., v$fmmsch) == 3 ~ 10,
        safe_numeric(., v$fmmsch) == 4 ~ 11,
        safe_numeric(., v$fmmsch) == 5 ~ 12,
        TRUE ~ NA_real_
      ),
      
      father_educ_years = case_when(
        safe_numeric(., v$fmfsch) == 1 ~ 0,
        safe_numeric(., v$fmfsch) == 2 ~ 6,
        safe_numeric(., v$fmfsch) == 3 ~ 10,
        safe_numeric(., v$fmfsch) == 4 ~ 11,
        safe_numeric(., v$fmfsch) == 5 ~ 12,
        TRUE ~ NA_real_
      ),
      
      # State
      state_code = safe_numeric(., v$hhstate),
      state = case_when(
        state_code == 1 ~ "NSW",
        state_code == 2 ~ "VIC",
        state_code == 3 ~ "QLD",
        state_code == 4 ~ "SA",
        state_code == 5 ~ "WA",
        state_code == 6 ~ "TAS",
        state_code == 7 ~ "NT",
        state_code == 8 ~ "ACT",
        TRUE ~ NA_character_
      ),
      
      # Weights
      analysis_weight = safe_numeric(., v$hhwtrps),
      population_weight = safe_numeric(., v$hhwtrp),
      longitudinal_weight = safe_numeric(., v$lnwtrp)  # For FE/panel analysis
    )
  
  # Networking index (composite of club_member + social_freq)
  # scale() requires full column, so computed after mutate
  wave_data <- wave_data %>%
    mutate(
      networking_index = {
        cm_z <- if (sd(club_member, na.rm = TRUE) > 0) as.numeric(scale(club_member)) else 0
        sf_z <- if (sd(social_freq, na.rm = TRUE) > 0) as.numeric(scale(social_freq)) else 0
        ifelse(!is.na(club_member) & !is.na(social_freq), cm_z + sf_z, NA_real_)
      }
    )
  
  cat("    - Networking: club_member non-NA =",
      sum(!is.na(wave_data$club_member)), 
      ", social_freq non-NA =", 
      sum(!is.na(wave_data$social_freq)), "\n")
  
  # Clean up raw data
  rm(raw)
  gc()
  
  # ============================
  # APPLY FILTERS — TWO-TIER APPROACH
  # ============================
  # 
  # Tier 1 (broad_sample): Everyone in age range with valid demographics.
  #   Includes unemployed, out-of-labor-force, missing wages.
  #   Used for: Heckman selection models, participation analysis,
  #   motherhood transition tracking, attrition diagnostics.
  #
  # Tier 2 (in_wage_sample): Employed with valid wages.
  #   Used for: Mincer equations, FE wage regressions.
  #   Flagged as in_wage_sample == 1 within the broad sample.
  
  cat("  Applying filters...\n")
  
  # Calculate wage quantiles (from employed only)
  wage_q005 <- quantile(wave_data$hourly_wage_clean, 0.005, na.rm = TRUE)
  wage_q995 <- quantile(wave_data$hourly_wage_clean, 0.995, na.rm = TRUE)
  
  # ---- TIER 1: Broad sample (includes non-employed) ----
  wave_broad <- wave_data %>%
    filter(
      # Demographics available
      !is.na(age), 
      age >= 21, 
      age <= 52,
      !is.na(female),
      !is.na(educ_years),
      !is.na(state)
    )
  
  cat("  Broad sample (Tier 1):", format(nrow(wave_broad), big.mark = ","), "observations\n")
  
  # ---- TIER 2: Wage sample flag within broad ----
  wave_broad <- wave_broad %>%
    mutate(
      in_wage_sample = case_when(
        # Must be employed
        !(esdtl_raw %in% c(1, 2)) ~ 0L,
        # Must have valid wages
        is.na(ln_hourly_wage_real) ~ 0L,
        hourly_wage_clean < wage_q005 | hourly_wage_clean > wage_q995 ~ 0L,
        # Must have valid hours
        is.na(hours_worked_clean) | hours_worked_clean <= 0 | hours_worked_clean > 100 ~ 0L,
        # Must have valid experience
        is.na(experience_years) | experience_years < 0 ~ 0L,
        # Must have weights
        is.na(analysis_weight) | analysis_weight <= 0 ~ 0L,
        # All conditions met
        TRUE ~ 1L
      )
    )
  
  n_wage <- sum(wave_broad$in_wage_sample == 1)
  n_nonwage <- sum(wave_broad$in_wage_sample == 0)
  cat("  Wage sample (Tier 2):", format(n_wage, big.mark = ","), "observations\n")
  cat("  Non-wage (available for selection):", format(n_nonwage, big.mark = ","), "observations\n")
  cat("  Employment rate in broad sample:", 
      round(100 * mean(wave_broad$esdtl_raw %in% c(1,2), na.rm = TRUE), 1), "%\n")
  cat("  Retention (wage/total):", round(100 * n_wage / nrow(wave_data), 1), "%\n")
  
  # ============================
  # EDUCATION DIAGNOSTICS (WAGE SAMPLE)
  # ============================
  
  cat("\n  Education Diagnostics (wage sample):\n")
  
  wave_wage <- wave_broad %>% filter(in_wage_sample == 1)
  
  # Education level distribution
  educ_dist <- wave_wage %>%
    count(education_level) %>%
    mutate(pct = round(100 * n / sum(n), 1))
  
  cat("    Education level distribution:\n")
  print(educ_dist, n = 20)
  
  # Education pathway distribution (top 20)
  pathway_dist <- wave_wage %>%
    count(education_pathway, sort = TRUE) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    head(20)
  
  cat("\n    Top 20 education pathways:\n")
  print(pathway_dist)
  
  # Field-specific education
  cat("\n    Field-specific education:\n")
  cat("      Nursing:", sum(wave_wage$has_nursing_any, na.rm = TRUE), 
      "(", round(100 * mean(wave_wage$has_nursing_any, na.rm = TRUE), 1), "%)\n")
  cat("      Teaching:", sum(wave_wage$has_teaching_any, na.rm = TRUE),
      "(", round(100 * mean(wave_wage$has_teaching_any, na.rm = TRUE), 1), "%)\n")
  cat("      Business:", sum(wave_wage$has_business_any, na.rm = TRUE),
      "(", round(100 * mean(wave_wage$has_business_any, na.rm = TRUE), 1), "%)\n")
  cat("      Trade:", sum(wave_wage$has_trade_any, na.rm = TRUE),
      "(", round(100 * mean(wave_wage$has_trade_any, na.rm = TRUE), 1), "%)\n")
  cat("      Technical:", sum(wave_wage$has_technical_any, na.rm = TRUE),
      "(", round(100 * mean(wave_wage$has_technical_any, na.rm = TRUE), 1), "%)\n")
  
  # NEW: Fertility history diagnostics (wage sample)
  cat("\n  Fertility History Diagnostics (wage sample):\n")
  cat("    ever_parent == 1:", sum(wave_wage$ever_parent == 1, na.rm = TRUE), "\n")
  cat("    ever_parent == 0:", sum(wave_wage$ever_parent == 0, na.rm = TRUE), "\n")
  cat("    ever_parent NA:", sum(is.na(wave_wage$ever_parent)), "\n")
  if(sum(wave_wage$female == 1) > 0) {
    cat("    Women - motherhood_status:\n")
    print(table(wave_wage$motherhood_status, useNA = "ifany"))
  }
  cat("    Birth events this wave:", sum(wave_wage$had_birth_this_wave == 1, na.rm = TRUE), "\n")
  
  # NEW: Selection diagnostics (broad vs wage sample)
  cat("\n  Selection Diagnostics (broad sample):\n")
  selection_by_group <- wave_broad %>%
    filter(female == 1, !is.na(ever_parent)) %>%
    group_by(motherhood_status) %>%
    summarise(
      N_broad = n(),
      N_wage = sum(in_wage_sample == 1),
      Employment_rate = round(100 * mean(esdtl_raw %in% c(1,2), na.rm = TRUE), 1),
      Wage_sample_rate = round(100 * mean(in_wage_sample == 1), 1),
      .groups = "drop"
    )
  cat("    Women by motherhood status (broad vs wage):\n")
  print(selection_by_group)
  
  # Clean up
  rm(wave_data, wave_wage)
  gc()
  
  # ---- DROP intermediate columns to save memory ----
  # Keep esdtl_raw for selection diagnostics, drop all other _raw columns
  raw_cols <- grep("_raw$", names(wave_broad), value = TRUE)
  raw_cols <- setdiff(raw_cols, "esdtl_raw")  # keep this one
  wave_broad <- wave_broad %>% dplyr::select(-any_of(raw_cols))
  
  cat("  Columns retained:", ncol(wave_broad), "\n")
  cat("  Memory:", format(object.size(wave_broad), units = "MB"), "\n")
  
  return(wave_broad)
}

# =====================================================================
# PROCESS ALL WAVES
# =====================================================================

process_all_waves <- function(data_dir, wpi_data) {
  
  cat("\n############################################\n")
  cat("# HILDA DATA PROCESSING - WAVES 15-24     #\n")
  cat("############################################\n")
  
  # ------------------------------------------------------------------
  # CHECKPOINT VERSION: bump this string whenever the wave processing
  # logic changes (variable construction, recoding, new variables).
  # Old checkpoints with a different version tag are ignored and the
  # wave is reprocessed from the raw .dta file.
  # Current change (v5): potential_experience fallback is NA-propagating
  # (missing age/educ_years -> NA, not 0); employed is NA-propagating;
  # raw reads are column-restricted. Earlier (v4): experience_years uses
  # ehtjb with biological ceiling (age-13) only; Mincer cap removed.
  # ------------------------------------------------------------------
  # v6: unified loader, waves 12-24, deflator column carried, all
  #     outputs written from one pass. (v5: NA-propagating experience
  #     and employed; column-restricted reads.)
  CHECKPOINT_VERSION <- "v7_wpi_fixed_w12_24"
  
  temp_dir <- file.path(data_dir, "wave_checkpoints")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  cat("  Checkpoint directory:", temp_dir, "\n")
  cat("  Checkpoint version:  ", CHECKPOINT_VERSION, "\n")
  
  # Warn if stale (unversioned) checkpoints exist from a previous run
  stale <- list.files(temp_dir, pattern = "^wave_\\d+\\.rds$", full.names = TRUE)
  if (length(stale) > 0) {
    cat("  WARNING:", length(stale), "unversioned checkpoint(s) found and will be ignored.\n")
    cat("  Delete them manually if you no longer need them:\n")
    cat("    unlink(list.files(\'" , temp_dir, "\', pattern=\'wave_\\\\d+\\.rds$\', full.names=TRUE))\n")
  }
  
  wave_files <- c()
  
  for(w in WAVES) {
    temp_file <- file.path(temp_dir, paste0("wave_", w, "_", CHECKPOINT_VERSION, ".rds"))
    
    # RESUME: skip only if a checkpoint with the CURRENT version tag exists
    if (file.exists(temp_file)) {
      cat(sprintf("  Wave %d: checkpoint found (%s), skipping processing\n", w, CHECKPOINT_VERSION))
      wave_files <- c(wave_files, temp_file)
      next
    }
    
    wave_result <- process_wave_complete(w, data_dir, wpi_data)
    
    if(!is.null(wave_result)) {
      # Save to disk immediately, free memory
      saveRDS(wave_result, temp_file)
      wave_files <- c(wave_files, temp_file)
      rm(wave_result)
    }
    
    gc()
    cat("  Memory freed. Files saved:", length(wave_files), "\n\n")
  }
  
  # Combine by reading back one at a time
  cat("\n===========================================\n")
  cat("Combining all waves from disk...\n")
  cat("===========================================\n")
  
  combined_data <- NULL
  for(f in wave_files) {
    cat("  Reading:", basename(f), "...")
    wave_chunk <- readRDS(f)
    if(is.null(combined_data)) {
      combined_data <- wave_chunk
    } else {
      combined_data <- dplyr::bind_rows(combined_data, wave_chunk)
    }
    rm(wave_chunk)
    gc()
    cat(" Total rows:", format(nrow(combined_data), big.mark = ","), "\n")
  }
  
  # Versioned checkpoints are retained for resume capability.
  # To force full reprocessing, either bump CHECKPOINT_VERSION above
  # or delete the checkpoint folder:
  #   unlink(file.path(data_dir, "wave_checkpoints"), recursive = TRUE)
  cat("  Versioned checkpoints retained in:", temp_dir, "\n")
  cat("  (Bump CHECKPOINT_VERSION or delete folder to force reprocessing)\n")
  
  # ================================================================
  # CUMMAX FIX: Enforce monotonicity in ever_parent
  # ================================================================
  # Back-switching detected: 88 individuals (0.54%) had ever_parent
  # go from 1→0 due to tchad reporting inconsistencies across waves.
  # Fix: within each person, once ever_parent = 1, it stays 1.
  # Uses base R ave() to avoid dplyr memory overhead.
  
  cat("\n--- Enforcing ever_parent monotonicity with cummax() ---\n")
  n_before <- sum(combined_data$ever_parent == 1, na.rm = TRUE)
  
  # Sort by person + wave (required for cummax to work forward in time)
  combined_data <- combined_data[order(combined_data$person_id, combined_data$wave), ]
  
  # Base R cummax via ave() — no tibble copies, minimal memory
  combined_data$ever_parent <- as.integer(ave(
    as.numeric(coalesce(combined_data$ever_parent, 0L)),
    combined_data$person_id,
    FUN = cummax
  ))
  gc()
  
  n_after <- sum(combined_data$ever_parent == 1, na.rm = TRUE)
  cat(sprintf("  ever_parent == 1: %s -> %s (gained %d from cummax)\n",
              format(n_before, big.mark = ","),
              format(n_after, big.mark = ","),
              n_after - n_before))
  
  # Re-derive all downstream parenthood variables after cummax fix
  combined_data$ever_mother <- ifelse(combined_data$female == 1, combined_data$ever_parent, NA_integer_)
  combined_data$ever_father <- ifelse(combined_data$female == 0, combined_data$ever_parent, NA_integer_)
  
  combined_data$motherhood_status <- with(combined_data, case_when(
    female == 0 ~ NA_character_,
    ever_parent == 1 & num_children_under15 > 0 ~ "active_mother",
    ever_parent == 1 & num_children_under15 == 0 ~ "post_mother",
    ever_parent == 0 ~ "never_mother",
    TRUE ~ NA_character_
  ))
  
  combined_data$fatherhood_status <- with(combined_data, case_when(
    female == 1 ~ NA_character_,
    ever_parent == 1 & num_children_under15 > 0 ~ "active_father",
    ever_parent == 1 & num_children_under15 == 0 ~ "post_father",
    ever_parent == 0 ~ "never_father",
    TRUE ~ NA_character_
  ))
  
  # Verify: no back-switching should remain
  n_bad <- sum(ave(
    as.numeric(combined_data$ever_parent),
    combined_data$person_id,
    FUN = function(x) any(diff(x) < 0, na.rm = TRUE)
  ) > 0) / max(table(combined_data$person_id))  # rough count
  # Cleaner check:
  ep_check <- tapply(combined_data$ever_parent, combined_data$person_id,
                     function(x) any(diff(x) < 0, na.rm = TRUE))
  n_bad_persons <- sum(ep_check, na.rm = TRUE)
  cat(sprintf("  Back-switching after fix: %d individuals (should be 0)\n", n_bad_persons))
  if (n_bad_persons > 0) warning("cummax() did not eliminate all back-switching!")
  
  # Report updated group counts
  cat("\n  Updated group counts after cummax:\n")
  cat("  ever_parent == 1:", sum(combined_data$ever_parent == 1, na.rm = TRUE), "\n")
  cat("  ever_parent == 0:", sum(combined_data$ever_parent == 0, na.rm = TRUE), "\n")
  cat("  Mothers ever:", sum(combined_data$ever_mother == 1, na.rm = TRUE), "\n")
  cat("  Never-mothers:", sum(combined_data$ever_mother == 0, na.rm = TRUE), "\n")
  cat("  Post-mothers:", sum(combined_data$motherhood_status == "post_mother", na.rm = TRUE), "\n")
  gc()
  
  # Add final variables
  combined_data <- combined_data %>%
    mutate(
      # Period indicators
      period = case_when(
        year %in% 2015:2019 ~ "Pre_COVID",
        year %in% 2020:2021 ~ "COVID",
        year %in% 2022:2024 ~ "Post_COVID"
      ),
      
      # Interaction terms (original)
      educ_female = educ_years * female,
      educ_mother = educ_years * is_mother,
      educ_treatment = educ_years * treatment,
      educ_precarious = educ_years * (precarious_score >= 2),
      
      # NEW: Interaction terms with absorbing motherhood
      educ_ever_mother = educ_years * coalesce(ever_mother, 0L),
      educ_ever_father = educ_years * coalesce(ever_father, 0L),
      
      # Centered education
      educ_centered = educ_years - mean(educ_years, na.rm = TRUE),
      educ_centered_sq = educ_centered^2,
      
      # ================================================================
      # GROUP CLASSIFICATIONS
      # ================================================================
      
      # PRIMARY (new): Ever-parent based — absorbing, no back-switching
      group_ever = case_when(
        female == 0 & ever_father == 1  ~ "fathers_ever",
        female == 0 & ever_father == 0  ~ "childless_men_ever",
        female == 1 & ever_mother == 1  ~ "mothers_ever",
        female == 1 & ever_mother == 0  ~ "never_mothers",
        TRUE ~ NA_character_
      ),
      
      # SENSITIVITY (current): Original children-under-15 based — time-varying
      group_current = case_when(
        female == 0 & num_children_under15 > 0 ~ "fathers",
        female == 0 & num_children_under15 == 0 ~ "childless_men",
        female == 1 & num_children_under15 > 0 ~ "mothers",
        female == 1 & num_children_under15 == 0 ~ "childless_women",
        TRUE ~ NA_character_
      ),
      
      # THREE-WAY WOMEN SPLIT (for decomposing contamination)
      group_3way_women = case_when(
        female == 0 ~ NA_character_,
        motherhood_status == "never_mother"  ~ "never_mothers",
        motherhood_status == "active_mother" ~ "active_mothers",
        motherhood_status == "post_mother"   ~ "post_mothers",
        TRUE ~ NA_character_
      )
    )
  
  # ============================
  # COMPREHENSIVE EDUCATION DIAGNOSTICS
  # ============================
  
  cat("\n############################################\n")
  cat("# COMPREHENSIVE EDUCATION DIAGNOSTICS      #\n")
  cat("############################################\n")
  
  # Overall education statistics
  cat("\n=== EDUCATION STATISTICS ===\n")
  educ_stats <- combined_data %>%
    summarise(
      N = n(),
      Mean = mean(educ_years, na.rm = TRUE),
      SD = sd(educ_years, na.rm = TRUE),
      Min = min(educ_years, na.rm = TRUE),
      Q1 = quantile(educ_years, 0.25, na.rm = TRUE),
      Median = median(educ_years, na.rm = TRUE),
      Q3 = quantile(educ_years, 0.75, na.rm = TRUE),
      Max = max(educ_years, na.rm = TRUE),
      Missing = sum(is.na(educ_years))
    )
  print(educ_stats)
  
  # Education level distribution
  cat("\n=== EDUCATION LEVEL DISTRIBUTION ===\n")
  educ_level_dist <- combined_data %>%
    count(education_level) %>%
    mutate(pct = round(100 * n / nrow(combined_data), 2)) %>%
    arrange(desc(n))
  print(educ_level_dist, n = 20)
  
  # Education pathway distribution
  cat("\n=== TOP 30 EDUCATION PATHWAYS ===\n")
  pathway_dist <- combined_data %>%
    count(education_pathway) %>%
    mutate(pct = round(100 * n / nrow(combined_data), 2)) %>%
    arrange(desc(n)) %>%
    head(30)
  print(pathway_dist)
  
  # Field-specific breakdown
  cat("\n=== FIELD-SPECIFIC EDUCATION ===\n")
  field_stats <- combined_data %>%
    summarise(
      Nursing_N = sum(has_nursing_any, na.rm = TRUE),
      Nursing_Pct = round(100 * mean(has_nursing_any, na.rm = TRUE), 2),
      Teaching_N = sum(has_teaching_any, na.rm = TRUE),
      Teaching_Pct = round(100 * mean(has_teaching_any, na.rm = TRUE), 2),
      Business_N = sum(has_business_any, na.rm = TRUE),
      Business_Pct = round(100 * mean(has_business_any, na.rm = TRUE), 2),
      Trade_N = sum(has_trade_any, na.rm = TRUE),
      Trade_Pct = round(100 * mean(has_trade_any, na.rm = TRUE), 2),
      Technical_N = sum(has_technical_any, na.rm = TRUE),
      Technical_Pct = round(100 * mean(has_technical_any, na.rm = TRUE), 2)
    )
  print(t(field_stats))
  
  # Education by gender
  cat("\n=== EDUCATION BY GENDER ===\n")
  gender_educ <- combined_data %>%
    group_by(female) %>%
    summarise(
      N = n(),
      Mean_Years = round(mean(educ_years, na.rm = TRUE), 2),
      Pct_Degree = round(100 * mean(has_degree, na.rm = TRUE), 1),
      Pct_Postgrad = round(100 * mean(has_postgrad, na.rm = TRUE), 1),
      Pct_Nursing = round(100 * mean(has_nursing_any, na.rm = TRUE), 1),
      Pct_Teaching = round(100 * mean(has_teaching_any, na.rm = TRUE), 1)
    )
  print(gender_educ)
  
  # Education data source quality
  cat("\n=== EDUCATION DATA SOURCE ===\n")
  source_dist <- combined_data %>%
    count(education_data_source) %>%
    mutate(pct = round(100 * n / nrow(combined_data), 2)) %>%
    arrange(desc(n))
  print(source_dist)
  
  # ============================
  # NEW: FERTILITY HISTORY DIAGNOSTICS
  # ============================
  
  cat("\n############################################\n")
  cat("# FERTILITY HISTORY DIAGNOSTICS            #\n")
  cat("############################################\n")
  
  cat("\n=== EVER-PARENT DISTRIBUTION ===\n")
  print(table(combined_data$ever_parent, combined_data$female, useNA = "ifany",
              dnn = c("ever_parent", "female")))
  
  cat("\n=== GROUP CLASSIFICATIONS ===\n")
  cat("\nAbsorbing (group_ever):\n")
  print(table(combined_data$group_ever, useNA = "ifany"))
  cat("\nTime-varying (group_current):\n")
  print(table(combined_data$group_current, useNA = "ifany"))
  cat("\nThree-way women (group_3way_women):\n")
  print(table(combined_data$group_3way_women, useNA = "ifany"))
  
  cat("\n=== CONTAMINATION CHECK ===\n")
  cat("Women classified 'childless' under current definition who are ever_parent==1:\n")
  contam <- combined_data %>%
    filter(female == 1, group_current == "childless_women") %>%
    summarise(
      N_total = n(),
      N_ever_parent = sum(ever_parent == 1, na.rm = TRUE),
      Pct_contaminated = round(100 * mean(ever_parent == 1, na.rm = TRUE), 1)
    )
  print(contam)
  
  cat("\n=== BIRTH EVENTS ACROSS PANEL ===\n")
  birth_events <- combined_data %>%
    group_by(year) %>%
    summarise(
      N_births = sum(had_birth_this_wave == 1, na.rm = TRUE),
      N_women_births = sum(had_birth_this_wave == 1 & female == 1, na.rm = TRUE),
      N_men_births = sum(had_birth_this_wave == 1 & female == 0, na.rm = TRUE)
    )
  print(birth_events)
  
  # Summary
  cat("\n############################################\n")
  cat("# PROCESSING COMPLETE                      #\n")
  cat("############################################\n")
  cat("Total observations (broad sample):", format(nrow(combined_data), big.mark = ","), "\n")
  cat("Wage sample observations:", format(sum(combined_data$in_wage_sample == 1), big.mark = ","), "\n")
  cat("Non-wage observations:", format(sum(combined_data$in_wage_sample == 0), big.mark = ","), "\n")
  cat("Unique individuals:", format(n_distinct(combined_data$person_id), big.mark = ","), "\n")
  cat("Unique individuals in wage sample:", 
      format(n_distinct(combined_data$person_id[combined_data$in_wage_sample == 1]), big.mark = ","), "\n")
  
  # Selection into employment by motherhood status
  cat("\n=== EMPLOYMENT SELECTION BY MOTHERHOOD STATUS ===\n")
  selection_overall <- combined_data %>%
    filter(female == 1, !is.na(motherhood_status)) %>%
    group_by(motherhood_status) %>%
    summarise(
      N_broad = n(),
      N_employed = sum(esdtl_raw %in% c(1,2), na.rm = TRUE),
      N_wage_sample = sum(in_wage_sample == 1),
      Employment_rate = round(100 * N_employed / N_broad, 1),
      Wage_sample_rate = round(100 * N_wage_sample / N_broad, 1),
      .groups = "drop"
    )
  print(selection_overall)
  
  return(combined_data)
}

# =====================================================================
# EXECUTE DATA PROCESSING
# =====================================================================

# Set directory
# data_dir is defined by config.R (sourced at the top); do not hardcode it here.

# Process all waves
cat("\nProcessing HILDA waves 15-24...\n")
analysis_sample <- process_all_waves(data_dir, wpi_data)

# =====================================================================
# SAVE OUTPUTS
# =====================================================================

# Full all-variables frame: large and rarely needed. Off by default;
# set SAVE_FULL <- TRUE to archive it (xz: smallest, slow to write).
SAVE_FULL <- FALSE
if (SAVE_FULL) saveRDS(analysis_sample, "hilda_analysis_full.rds", compress = "xz")

panel_full <- analysis_sample   # no disk round trip


# =====================================================================
# NAME-BASED VARIABLE SELECTION
# =====================================================================
# Using dplyr::select with any_of() for variables that may not exist
# in all configurations (e.g., mechanism variables created in analysis script)

analysis_sample <- panel_full %>%
  dplyr::select(
    # === Identifiers & panel structure ===
    person_id, any_of("xwaveid"), wave, year,
    
    # === Demographics ===
    age, age_sq, female, married,
    
    # === Wages (outcome) ===
    ln_hourly_wage, ln_hourly_wage_real,
    hourly_wage_clean, real_hourly_wage, deflator,
    weekly_earnings, hours_worked_clean,
    annual_income, childcare_cost,
    
    # === Education ===
    educ_years, educ_sq, education_level, education_pathway,
    has_degree, has_postsecondary, has_diploma_plus, has_postgrad,
    has_nursing_any, has_teaching_any, has_business_any,
    has_trade_any, has_technical_any,
    education_data_source,
    
    # === Experience & tenure ===
    experience_years, experience_sq, actual_experience,
    tenure_employer, tenure_occupation,
    
    # === Employment characteristics ===
    employed, fulltime, parttime, married_num, fulltime_num,
    any_of("esdtl_raw"),
    in_wage_sample,
    casual, permanent,
    public_sector, union_member,
    firm_size,
    
    # === Job quality ===
    has_holiday_leave, has_sick_leave, leave_none,
    has_standard_hours, has_shift_work, has_irregular_hours,
    underemployed, overemployed,
    precarious_score, precarious_category,
    flexibility_any, treatment,
    
    # === Job satisfaction & workplace quality (NEW) ===
    job_satisfaction, pay_fairness, job_initiative,
    job_security, job_worry_future, job_autonomy,
    looking_for_new_job,
    
    # === Occupation & industry ===
    occupation_major, occupation_major_code,
    industry_major, industry_major_code,
    
    # === Training ===
    any_of(c("any_training", "total_training_hours",
             "promotion_training", "maintenance_training")),
    total_training_days, avg_training_hours, total_courses,
    training_own_time, training_paid_time,
    jttpewt_any, jttpeot_any,   # place-of-employment-only indicators
    
    # === Networking / social capital (created in analysis script) ===
    any_of(c("club_member", "social_freq", "social_freq_01",
             "high_social", "networking_index")),
    
    # === Children & care (original — time-varying) ===
    num_dependent_children, num_children_under15,
    num_children_0to4, num_children_5to9, num_children_10to14,
    num_children_15to24,
    has_young_children, active_caregiver, is_mother, is_father,
    housework_hours, care_intensity,
    
    # === Non-resident children (original) ===
    has_nonresident_children_under18, has_any_nonresident_children,
    
    # === Fertility history (NEW — absorbing) ===
    total_children_ever_had,
    num_own_resident_children, num_own_nonresident_children, total_own_children,
    age_youngest_own_child, total_children_died,
    num_step_foster_grand,
    
    # Own children by age group
    own_res_children_0to4, own_res_children_5to14,
    own_res_children_15to24, own_res_children_25plus,
    own_nonres_children_0to4, own_nonres_children_5to14,
    own_nonres_children_15to24, own_nonres_children_25plus,
    
    # Child ages (for age at first birth)
    age_resident_child_1, age_resident_child_2, age_resident_child_3,
    age_nonresident_child_1, age_nonresident_child_2,
    
    # === Absorbing parenthood indicators (NEW) ===
    ever_parent, ever_mother, ever_father,
    motherhood_status, fatherhood_status,
    approx_age_at_first_birth, age_oldest_child,
    
    # === Group classifications (NEW) ===
    group_ever, group_current, group_3way_women,
    
    # === Parenting attitudes (NEW) ===
    parenting_resp, parenting_hard, parenting_tired,
    parenting_trapped, parenting_fairshare,
    
    # === Pregnancy (NEW — fertility module waves only) ===
    currently_pregnant, partner_pregnant, self_or_partner_pregnant,
    expected_birth_month, expected_birth_year,
    
    # === Life events (NEW — all waves) ===
    event_birth, event_married, event_separated,
    event_job_change, event_fired, event_moved,
    had_birth_this_wave, got_married_this_wave,
    separated_this_wave, changed_job_this_wave,
    
    # === Other parent / partner (NEW) ===
    other_parent_employment, living_with_partner_last,
    any_of(c("household_id", "partner_id")),
    
    # === Parental background (instruments) ===
    mother_educ_years, father_educ_years,
    
    # === Geography & weights ===
    state, analysis_weight, population_weight, longitudinal_weight,
    
    # === Derived interactions ===
    educ_centered, educ_centered_sq,
    educ_female, educ_mother, educ_treatment, educ_precarious,
    educ_ever_mother, educ_ever_father,
    
    # === Period ===
    period
  )



# 
# Fix education_level to have unique categories
analysis_sample <- analysis_sample %>%
  mutate(
    education_level = case_when(
      educ_years >= 19.5 ~ "Doctorate",
      educ_years >= 18 & educ_years < 19.5 ~ "Masters",
      educ_years >= 17 & educ_years < 18 ~ "Graduate_Diploma_Certificate",
      educ_years >= 16.5 & educ_years < 17 ~ "Bachelor_Honours",
      educ_years >= 16 & educ_years < 16.5 ~ "Bachelor",
      educ_years >= 15 & educ_years < 16 ~ "Advanced_Diploma",
      educ_years >= 14.5 & educ_years < 15 ~ "Diploma",
      educ_years >= 14 & educ_years < 14.5 ~ "Diploma_Other",  # Or "Registered_Nurse" if that's what it represents
      educ_years >= 13.5 & educ_years < 14 ~ "Certificate_IV",
      educ_years >= 13 & educ_years < 13.5 ~ "Certificate_III",
      educ_years >= 12.5 & educ_years < 13 ~ "Certificate_II",
      educ_years == 12 ~ "Year_12",                # MUST precede the Cert_I band:
      educ_years > 12 & educ_years < 12.5 ~ "Certificate_I",  # first-match would
      # otherwise swallow all 12s
      educ_years == 11 ~ "Year_11",
      educ_years == 10 ~ "Year_10",
      educ_years == 9 ~ "Year_9",
      educ_years == 8 ~ "Year_8",
      educ_years == 7 ~ "Year_7",
      educ_years == 6 ~ "Primary_Complete",
      educ_years == 5 ~ "Primary_Incomplete",
      TRUE ~ "Other"
    )
  )

# Verify no duplicates
check_duplicates <- analysis_sample %>%
  group_by(education_level) %>%
  summarise(
    mean_years = mean(educ_years),
    min_years = min(educ_years),
    max_years = max(educ_years),
    n = n()
  ) %>%
  arrange(desc(mean_years))

print(check_duplicates)


education_order <- c(
  "Other",
  "Primary_Incomplete", "Primary_Complete",
  "Year_7", "Year_8", "Year_9",
  "Year_10", "Year_11", "Year_12",
  "Certificate_I", "Certificate_II", "Certificate_III", "Certificate_IV",
  "Diploma_Other", "Diploma", "Advanced_Diploma",
  "Bachelor", "Bachelor_Honours",
  "Graduate_Diploma_Certificate", "Masters", "Doctorate"
)
# guard: every produced category must be a level, or factor() silently NAs it
stopifnot(all(unique(analysis_sample$education_level) %in% education_order))

analysis_sample$education_level <- factor(
  analysis_sample$education_level,
  levels = education_order,
  ordered = TRUE
)

education_pathway_order <- c(
  "Primary_Incomplete", "Primary_Complete", "Year_7", "Year_8", "Year_9",
  "Below_Year_10", "Year_10", "Year_11", "Year_12",
  "Certificate_II", "Cert_III_Trade", "Cert_III_Technical", "Cert_III_Nursing", "Certificate_III",
  "Cert_IV_Trade", "Cert_IV_Technical", "Certificate_IV",
  "Diploma_Trade", "Diploma_Teaching", "Diploma_Nursing", "Diploma", "Advanced_Diploma", "AdvDiploma_Trade",
  "Bachelor_General", "Bachelor_Honours_General", "Bachelor_Nursing", "Bachelor_Teaching",
  "GradDipCert_General", "GradDipCert_Nursing", "GradDipCert_Teaching",
  "Masters_General", "Masters_Teaching",
  "Doctorate_General", "Doctorate_Teaching",
  "Registered_Nurse"
)


analysis_sample$education_pathway <- factor(
  analysis_sample$education_pathway,
  levels = education_pathway_order,
  ordered = TRUE
)

# =====================================================================
# WRITE ALL ANALYSIS DATASETS (single source of truth)
# =====================================================================
stopifnot(min(analysis_sample$wave) <= 12, max(analysis_sample$wave) >= 24)

# (1) Full-window panel for MASTER_hh.R Parts A-F (name kept for
#     drop-in compatibility with the master's readRDS call)
saveRDS(analysis_sample, "hilda_panel_data_W12_W24_slim.rds")
cat(sprintf("Saved hilda_panel_data_W12_W24_slim.rds: %s rows (waves %d-%d)\n",
            format(nrow(analysis_sample), big.mark = ","),
            min(analysis_sample$wave), max(analysis_sample$wave)))

# (2) W15-24 subset for PART G / PART H (clean Reform 1 identification)
ext_sample <- dplyr::filter(analysis_sample, wave >= 15)
saveRDS(ext_sample, "hilda_panel_data_extended.rds")
cat(sprintf("Saved hilda_panel_data_extended.rds: %s rows (waves %d-%d)\n",
            format(nrow(ext_sample), big.mark = ","),
            min(ext_sample$wave), max(ext_sample$wave)))
rm(ext_sample)

# =====================================================================
# SLIM TRAINING PANEL (input for training_event_study.R)
# =====================================================================
# Only the columns the training event study needs: identifiers, cohort
# inputs (female, ever_parent), employment, training outcomes, weights,
# and a couple of context variables. Indicators/counts stored as integer
# (halves in-memory size and compresses better); xz for smallest file.
training_vars <- c(
  "person_id", "xwaveid", "wave", "year", "age", "female", "married",
  "ever_parent", "employed",
  "hours_worked_clean", "ln_hourly_wage_real", "experience_years",
  "analysis_weight", "longitudinal_weight",
  "total_training_days", "avg_training_hours", "total_courses",
  "training_own_time", "training_paid_time",
  "jttpewt_any", "jttpeot_any"
)
int_vars <- c("wave", "year", "female", "married", "ever_parent", "employed",
              "total_courses", "training_own_time", "training_paid_time",
              "jttpewt_any", "jttpeot_any")

training_panel <- analysis_sample %>%
  dplyr::select(dplyr::any_of(training_vars)) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(int_vars),
                              ~ as.integer(round(as.numeric(.x)))))

missing_tv <- setdiff(training_vars, names(training_panel))
if (length(missing_tv))
  cat("NOTE: training panel missing columns:",
      paste(missing_tv, collapse = ", "), "\n")

saveRDS(training_panel, "hilda_panel_data_training.rds", compress = "xz")
cat(sprintf("Saved hilda_panel_data_training.rds: %d rows x %d cols\n",
            nrow(training_panel), ncol(training_panel)))