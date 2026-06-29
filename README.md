# Replication Package

## Lost Experience and the Motherhood Wage Penalty: A Dynamic Account of Household Specialisation in Australia

Ainura Tursunalieva (CSIRO) and Ratbek Dzhumashev (Monash University)

This README follows the Social Science Data Editors template. It documents the
data, the code, the computing environment, and the steps an independent
researcher needs to reproduce every table and figure in the paper. The analysis
data are confidential and are not included in this package; the Data
availability and provenance section below explains how to obtain them.

---

## 1. Overview

The code in this repository constructs the analysis sample from the raw HILDA
release and reproduces all exhibits in the paper. The pipeline is written in R.
A single top-level script, `run_all.R`, executes the full sequence in order:
sample construction, the main analysis (Parts A to F), the FTB-B reform
instrumental-variables analysis, the policy-economics extension, the
child-penalty and training event studies, and all publication figures.

The reference run completed in approximately 4 hours 15 minutes on the machine
described in Section 4.

---

## 2. Data availability and provenance statements

### 2.1 Statement of rights

The authors certify that they have legitimate access to, and permission to use,
the data used in this manuscript, obtained under a confidentiality deed with the
Australian Data Archive. The authors did not obtain permission to redistribute
the data, and therefore the data are not included in this package.

### 2.2 License for the data

The HILDA Survey data are confidential and are governed by the Australian Data
Archive confidentiality deed and the conditions of the Department of Social
Services. They cannot be redistributed. No part of the source data is included
in this repository.

### 2.3 Summary of data availability

The data cannot be made publicly available. They are available free of charge to
approved researchers through a restricted-access application process described
below.

### 2.4 Data Availability Statement

The analysis uses unit-record data from the HILDA Survey, General Release 24
(DOI 10.26193/6M1BMR), managed by the Melbourne Institute and funded by the
Australian Government Department of Social Services (DSS). The data are
confidential and cannot be redistributed. They are available free of charge to
approved researchers by application to the DSS Longitudinal Studies Dataverse
within the Australian Data Archive (https://dataverse.ada.edu.au/dataverse/hilda),
subject to a signed confidentiality deed; the General Release is available to
researchers outside Australia. Approval typically takes a few weeks. All code
required to construct the analysis sample from the raw release and to reproduce
every table and figure is provided in this repository; the raw data are not
included.

Note to the editor and Data Editor: access is restricted under a confidentiality
deed, so the authors cannot provide the microdata even privately. The Data Editor
can obtain identical data by signing the same confidentiality deed with DSS
through the Australian Data Archive. The authors commit to preserving the code
for no less than five years and to assisting with reasonable clarification
requests.

### 2.5 Data citation

Department of Social Services; Melbourne Institute of Applied Economic and Social
Research. 2024. *The Household, Income and Labour Dynamics in Australia (HILDA)
Survey, General Release 24.* ADA Dataverse. https://doi.org/10.26193/6M1BMR.

### 2.6 Details on the data source

| Item | Detail |
| --- | --- |
| Dataset | HILDA Survey, General Release 24, Waves 12 to 24 (2012 to 2024) |
| Provider | Australian Data Archive (ADA), DSS Longitudinal Studies Dataverse |
| DOI | 10.26193/6M1BMR |
| Access | Restricted; free to approved researchers under a confidentiality deed |
| Format expected | Stata `.dta` wave files from the raw release |
| Location in package | Place the raw release in `data/raw/`; it is gitignored and not distributed |

The HILDA release is the only source dataset. No other external data are used.

---

## 3. Dataset list

The pipeline builds three derived analysis panels from the raw release. These are
intermediate files, are gitignored, and are not distributed. They are written to
`data/analysis/` by the loader.

| File | Built by | Used by |
| --- | --- | --- |
| `hilda_panel_data_W12_W24_slim.rds` | `00_load_hilda.R` | Parts A to F (`01_master_hh.R`) |
| `hilda_panel_data_extended.rds` | `00_load_hilda.R` | Parts G and H |
| `hilda_panel_data_training.rds` | `00_load_hilda.R` | training event study |

---

## 4. Computational requirements

### 4.1 Software

The reference results were produced with R 4.5.2 (2025-10-31) on Windows 11 x64.
The full environment is recorded in `sessionInfo.txt`, and pinned versions are in
`renv.lock`. Restore the environment with:

```r
renv::restore()
```

If `renv` is not used, install dependencies with `source("PACKAGES.R")`, which
installs and version-checks the CRAN packages: tidyverse, fixest, lmtest,
sandwich, patchwork, DoubleML, mlr3, mlr3learners, ranger, glmnet, grf,
car, AER, MASS, and future.

### 4.2 Controlled randomness

The random seed is set to 42 in `01_master_hh.R`, before the double machine
learning and causal-forest steps. The causal forest and the cross-fitting folds
are the only sources of randomness.

### 4.3 Memory and runtime

The reference run completed in approximately 4 hours 15 minutes (start
2026-06-15 19:49, finish 2026-06-16 00:05). Peak memory was not separately
recorded; a workstation with at least 16 GB of RAM is recommended because the
causal forest and double machine learning steps are memory intensive. The thread
count is detected automatically: on a SLURM batch job the pipeline uses
`SLURM_CPUS_PER_TASK`, and interactively it caps at four threads. Override with
the `MASTER_HH_THREADS` environment variable.

---

## 5. Description of programs and code

Configuration and orchestration live at the repository root. The numbered
scripts in `code/` run in sequence. The ordering matters: Part H requires
in-memory objects from Part G and the results file written by the child-penalty
event study, so the event studies run before the policy extension.

```
repo-root/
├── README.md                  this file
├── LICENSE.txt                BSD-3-Clause for code; CC-BY-4.0 for documentation
├── .gitignore
├── renv.lock                  pinned R and package versions
├── sessionInfo.txt            environment record from the reference run
├── config.R                   sets data_dir and output_dir in one place
├── run_all.R                  master runner; sources code/ scripts in order
├── code/
│   ├── 00_load_hilda.R        raw HILDA waves to the three analysis panels
│   ├── 01_master_hh.R         Parts A to F: main couple analysis, housework
│   │                          mediation, panel IV (diagnostic), couple DML and
│   │                          causal forest, robustness battery, Oster bounds
│   ├── 02_cs_event_study.R    child-penalty event study (Sun and Abraham 2021)
│   ├── 03_training_event_study.R  training event study around first birth
│   ├── 04_part_g_ftbb_iv.R    FTB-B 2015 reform simulated-instrument IV
│   ├── 05_part_h_policy.R     policy economics of the FTB-B reform
│   └── 06_figures.R           all publication figures
├── data/
│   ├── raw/                   empty; replicator places HILDA Release 24 here
│   │   └── README.md          which waves, how to obtain, the DOI
│   └── analysis/              derived panels written by 00_load_hilda.R (gitignored)
├── output/
│   ├── tables/
│   └── figures/
└── docs/
    └── variables.md           HILDA variables used and their derivations
```

This layout renames the working scripts to the numbered scheme above. The
correspondence to the original filenames is: `00_load_hilda.R` is the HILDA
loader; `01_master_hh.R` is `MASTER_hh.R`; `02_cs_event_study.R` is
`CS_event_study.R`; `03_training_event_study.R` is `training_event_study.R`;
`04_part_g_ftbb_iv.R` is `MASTER_PART_G_FTBB_IV.R`; `05_part_h_policy.R` is
`Master_part_h_policy.R`; and `06_figures.R` is `RUN_FIGURES_hh.R`.

---

## 6. Instructions to replicators

1. Obtain HILDA General Release 24 from the Australian Data Archive, following
   Section 2. Place the raw wave files in `data/raw/`.
2. Restore the environment with `renv::restore()`, or run `source("PACKAGES.R")`.
3. Set the data and output paths in `config.R` if they differ from the defaults.
4. From the repository root, run the full pipeline:

   ```sh
   Rscript run_all.R
   ```

   `run_all.R` builds the analysis panels, runs the full analysis, and writes all
   tables and figures to `output/`. It runs in a single R session because the
   figure script reuses objects built earlier in the run.

To run a single stage, source `config.R` first, then the numbered script. Stages
02 through 06 assume the analysis panels from stage 00 are on disk, and stage 05
assumes stages 02 and 04 have already run in the same session.

---

## 7. List of exhibits and the code that produces them

Table numbers in the paper follow order of appearance; the labels below are the
manuscript's internal labels. Confirm the final numbering against the compiled
manuscript. All figures are written by `06_figures.R`, which reads result files
produced by the upstream stages noted in the dependency column.

### Tables

| Manuscript label | Exhibit | Produced by |
| --- | --- | --- |
| `tab:panel_descriptives` | Descriptive statistics by motherhood group | `01_master_hh.R` (descriptives) |
| `tab:couple_descriptives` | Within-couple labour-market outcomes | `01_master_hh.R` (Part A) |
| `tab:hypotheses` | Empirical hypotheses and estimands | Set in the manuscript; not computed |
| `tab:housework_4grp` | Housework hours by gender and parenthood | `01_master_hh.R` (Part B) |
| `tab:hw_within_fe` | Within-person hours change at parenthood | `01_master_hh.R` (Part B) |
| `tab:specialisation` | Specialisation test, cross-partner elasticity | `01_master_hh.R` (Part A) |
| `tab:event_study` | Child-penalty event study, women minus men | `02_cs_event_study.R` |
| `tab:iv_main` | FTB-B first stage and contemporaneous-wage IV | `04_part_g_ftbb_iv.R` |
| (training participation table) | Training around first birth, women minus men | `03_training_event_study.R` |
| `tab:within_couple` | Within-couple DML gender gap in complementarity | `01_master_hh.R` (Part D) |
| `tab:main_beta2` | Education-experience complementarity by group | `01_master_hh.R` (Parts D, E) |
| `tab:heckman` | Heckman selection correction | `01_master_hh.R` (Part E) |
| `tab:transitioner_bound` | Within-person complementarity by transition type | `01_master_hh.R` (Part E) |
| `tab:definition_sensitivity` | Motherhood-definition sensitivity | `01_master_hh.R` (Part E) |
| `tab:robustness_master` | Master robustness summary | `01_master_hh.R` (Part E) |
| `tab:ftbb_first_stage` | FTB-B first stages and implied elasticity | `04_part_g_ftbb_iv.R`, `05_part_h_policy.R` |

### Figures

| Figure file | Exhibit | Upstream dependency |
| --- | --- | --- |
| `fig2_event_study_beta2.pdf` | Complementarity event study | `event_study_beta2.rds` (Part E) |
| `fig3_employment_birth.pdf` | Employment around birth | `01_master_hh.R` |
| `fig5C1_gender_gap_vanishes.pdf` | Within-couple gender gap | Part D objects |
| `fig7_robustness_forest.pdf` | Robustness forest | `robustness_forest_data.rds` (Part E) |
| `fig_event_study.pdf` | Child-penalty event study | `cs_event_study_results.rds` (stage 02) |
| `fig_training_event_study.pdf` | Training event study | `training_event_study_results.rds` (stage 03) |
| `fig_iv_coefplot.pdf` | FTB-B IV coefficients | `ftbb_reform1_iv_results.rds` (stage 04) |
| `fig_hw_within_fe.pdf` | Within-couple housework allocation | `01_master_hh.R` (Part B) |

---

## 8. License

The code is released under the BSD-3-Clause license and the documentation under
CC-BY-4.0; see `LICENSE.txt`. The license covers this repository only and does
not extend to the HILDA data, which remain subject to the ADA confidentiality
deed.

---

## 9. References

Sun, Liyang, and Sarah Abraham. 2021. "Estimating Dynamic Treatment Effects in
Event Studies with Heterogeneous Treatment Effects." *Journal of Econometrics*
225 (2): 175-199.

Oster, Emily. 2019. "Unobservable Selection and Coefficient Stability: Theory and
Evidence." *Journal of Business and Economic Statistics* 37 (2): 187-204.

Department of Social Services; Melbourne Institute of Applied Economic and Social
Research. 2024. *The Household, Income and Labour Dynamics in Australia (HILDA)
Survey, General Release 24.* ADA Dataverse. https://doi.org/10.26193/6M1BMR.
