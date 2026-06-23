# Replication Package

## Lost Experience and the Motherhood Wage Penalty: A Dynamic Account of Household Specialisation in Australia

Ainura Tursunalieva (CSIRO) and Ratbek Dzhumashev (Monash University)

This README follows the Social Science Data Editors template. It documents the
data, the code, the computing environment, and the steps an independent
researcher needs to reproduce every table and figure in the paper. The analysis
data are confidential and are not included in this package; Section 2 explains
how to obtain them. The package also ships a synthetic dataset so the entire
pipeline can be run and inspected without the confidential data, as described in
Section 1a.

---

## 1. Overview

The code in this repository constructs the analysis sample from the raw HILDA
release and reproduces all exhibits in the paper. The pipeline is written in R.
A single top-level script, `run_all.R`, executes the full sequence in order:
sample construction, the main analysis (Parts A to F), the FTB-B reform
instrumental-variables analysis, the policy-economics extension, the
child-penalty and training event studies, and all publication figures.

The reference run on the confidential data completed in approximately 4 hours
15 minutes on the machine described in Section 4. In synthetic mode (Section 1a)
the pipeline completes in a few minutes.

---

## 1a. Running the pipeline (real data and synthetic mode)

The pipeline runs in two modes, selected by the `DATA_MODE` setting in
`config.R`, or by the `DATA_MODE` environment variable, which overrides the file
without editing it. The default is synthetic, so a fresh clone runs end to end
with no data and no manual setup.

`DATA_MODE = "synthetic"` (default). The three analysis panels are built on the
fly by `00b_make_synthetic_hilda.R` from a fixed seed, and the full analysis and
all figures run against them. This requires no access to HILDA and reaches every
code path: the couple analysis, the housework mediation, the panel and FTB-B
instrumental-variables steps, the double machine learning and causal forest, the
robustness battery, the event studies, and the figures. Use this mode to verify
that the code installs and runs, to read the pipeline, or to develop against it.

`DATA_MODE = "real"`. The pipeline uses the panels built from the restricted
HILDA release by `00_load_hilda.R`. Only the authors and approved researchers
who have obtained the data under the confidentiality deed (Section 2) can use
this mode. Build the panels once with `Rscript 00_load_hilda.R`, then run the
rest.

To run the full pipeline from the repository root:

```sh
Rscript run_all.R                      # synthetic by default
DATA_MODE=real Rscript run_all.R       # uses the restricted panels (authors)
```

In synthetic mode `run_all.R` builds the panels if they are not already on disk,
then runs `01_master.R` (which sources the event studies and the FTB-B IV) and
`06_figures.R` in one session. In real mode it stops with an instruction if the
panels are missing, rather than overwriting them.

### Synthetic data are illustrative, not the paper's results

The synthetic panels are random data generated to match the schema of the real
analysis panels: the same variables, types, and factor levels, the same panel
and couple structure, and a built-in motherhood interruption and
education-by-experience wage gradient so the estimates carry sensible signs.
They exist so the code can be run and inspected without the confidential data.

Any number produced in synthetic mode is an artefact of the random generator.
These numbers are not estimates of the paper's quantities, will change with the
seed or sample size, and must not be reported or cited. Every result, table, and
figure in the manuscript comes from `DATA_MODE = "real"` on HILDA General
Release 24. The synthetic panels and all generated outputs are gitignored and
are not part of the distributed package.

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
| Location in package | Place the raw release where `config.R`'s `data_dir` points (set `HILDA_DATA_DIR`); it is gitignored and not distributed |

The HILDA release is the only source dataset. No other external data are used.

---

## 3. Dataset list

The pipeline builds three derived analysis panels. In real mode they are written
by `00_load_hilda.R` from the raw release; in synthetic mode they are written by
`00b_make_synthetic_hilda.R` from a fixed seed. Both modes write the same three
filenames with the same schema, so every downstream stage is identical across
modes. The panels, and the intermediate result caches the stages write
(`cs_event_study_results.rds`, `ftbb_reform1_iv_results.rds`,
`event_study_beta2.rds`, `robustness_forest_data.rds`, `oster_bounds.rds`, and
others), are all `.rds`, are gitignored, and are not distributed.

| File | Built by (real / synthetic) | Used by |
| --- | --- | --- |
| `hilda_panel_data_W12_W24_slim.rds` | `00_load_hilda.R` / `00b_make_synthetic_hilda.R` | Parts A to F (`01_master.R`) |
| `hilda_panel_data_extended.rds` | `00_load_hilda.R` / `00b_make_synthetic_hilda.R` | Parts G and H |
| `hilda_panel_data_training.rds` | `00_load_hilda.R` / `00b_make_synthetic_hilda.R` | training event study |

---

## 4. Computational requirements

### 4.1 Software

The reference results were produced with R 4.5.2 (2025-10-31) on Windows 11 x64.
The full environment is recorded in `sessionInfo.txt`. Install and version-check
dependencies with:

```r
source("PACKAGES.R")
```

`PACKAGES.R` installs and version-checks the CRAN packages: tidyverse, fixest,
lmtest, sandwich, patchwork, DoubleML, mlr3, mlr3learners, ranger, glmnet, grf,
ivreg, car, AER, MASS, future, and haven. If you maintain a `renv` lockfile,
`renv::restore()` is an equivalent route.

One dependency is not on CRAN. The Contextual Robustness test in `01_master.R`
uses `ivcrtest::iv_cr_test()` from the custom package `ivcrtest` (reference
version 0.1.0), which depends on `ivreg`. Install `ivcrtest` from source before
running; the relevant section stops with an instruction if it is absent, and it
is the only section that needs it.

Note on namespace masking: `config.R` rebinds the dplyr verbs that other attached
packages mask (`select`, `filter`, `mutate`, `between`, and others) to their
dplyr definitions. `run_all.R` sources `config.R` after `PACKAGES.R` so these
bindings take effect for the whole run.

### 4.2 Controlled randomness

In real mode the seed is set to 42 in `01_master.R`, before the double machine
learning and causal-forest steps; the causal forest and the cross-fitting folds
are the only sources of randomness. In synthetic mode the data generator uses its
own fixed seed in `00b_make_synthetic_hilda.R`, so the synthetic panels are
reproducible across runs.

### 4.3 Memory and runtime

The reference run on the confidential data completed in approximately 4 hours
15 minutes. A workstation with at least 16 GB of RAM is recommended because the
causal forest and double machine learning steps are memory intensive. The thread
count is detected automatically: on a SLURM batch job the pipeline uses
`SLURM_CPUS_PER_TASK`, and interactively it caps at four threads. Override with
the `MASTER_HH_THREADS` environment variable. In synthetic mode, with a smaller
default sample, the pipeline completes in a few minutes; set `N_PERSONS` to
adjust the synthetic sample size.

---

## 5. Description of programs and code

Configuration and orchestration live alongside the analysis scripts in the
repository root; the scripts source each other by name, so they must stay in one
folder. The numbered scripts run in sequence. The ordering matters: Part H
requires in-memory objects from Part G and the results file written by the
child-penalty event study, so the event studies run before the policy extension.

```
repo-root/
├── README.md                     this file
├── LICENSE                       BSD-3-Clause for code
├── .gitignore
├── sessionInfo.txt               environment record from the reference run
├── config.R                      DATA_MODE, paths, and dplyr namespace guards
├── run_all.R                     master runner: PACKAGES, config, data, analysis, figures
├── PACKAGES.R                    install and version-check dependencies
├── 00_load_hilda.R               real mode: raw HILDA waves to the three analysis panels
├── 00b_make_synthetic_hilda.R    synthetic mode: seed-based stand-in panels (same schema)
├── 01_master.R                   Parts A to F: main couple analysis, housework
│                                 mediation, panel IV and CR test, couple DML and
│                                 causal forest, robustness battery, Oster bounds.
│                                 Sources 02, 04, 05, 03 in that order.
├── 02_cs_event_study.R           child-penalty event study (Sun and Abraham 2021)
├── 03_training_event_study.R     training event study around first birth
├── 04_part_g_ftbb_iv.R           FTB-B 2015 reform simulated-instrument IV (Part G)
├── 05_part_h_policy.R            policy economics of the FTB-B reform (Part H)
└── 06_figures.R                  all publication figures
```

Outputs are written to `output/figures/` and `output/tables/`, and per-run console
logs to `run_logs/`. These directories, the analysis panels, and the intermediate
`.rds` caches are all gitignored.

Correspondence to the original working filenames, for anyone holding an earlier
copy: `01_master.R` was `MASTER_hh.R`; `02_cs_event_study.R` was
`CS_event_study.R`; `03_training_event_study.R` was `training_event_study.R`;
`04_part_g_ftbb_iv.R` was `MASTER_PART_G_FTBB_IV.R`; `05_part_h_policy.R` was
`Master_part_h_policy.R`; and `06_figures.R` was `RUN_FIGURES_hh.R`.

---

## 6. Instructions to replicators

To run on the confidential data and reproduce the paper:

1. Obtain HILDA General Release 24 from the Australian Data Archive, following
   Section 2. Place the raw wave files where `config.R`'s `data_dir` points (set
   the `HILDA_DATA_DIR` environment variable, or edit `config.R`).
2. Install dependencies with `source("PACKAGES.R")`, and install the custom
   `ivcrtest` package from source.
3. From the repository root, build the panels and run the full pipeline:

   ```sh
   Rscript 00_load_hilda.R            # build the three analysis panels
   DATA_MODE=real Rscript run_all.R   # run the analysis and write all exhibits
   ```

   `run_all.R` runs the full analysis and writes all tables and figures to
   `output/`. It runs in a single R session because the figure script reuses
   objects built earlier in the run.

To run without the confidential data, using the synthetic panels, simply run
`Rscript run_all.R`; synthetic is the default mode and the panels are built
automatically. See Section 1a for the caveat on synthetic results.

To run a single stage, source `config.R` first, then the numbered script. Stages
02 through 06 assume the analysis panels are on disk, and stage 05 assumes stages
02 and 04 have already run in the same session.

---

## 7. List of exhibits and the code that produces them

Table numbers in the paper follow order of appearance; the labels below are the
manuscript's internal labels. Confirm the final numbering against the compiled
manuscript. All figures are written by `06_figures.R`, which reads result files
produced by the upstream stages noted in the dependency column.

### Tables

| Manuscript label | Exhibit | Produced by |
| --- | --- | --- |
| `tab:panel_descriptives` | Descriptive statistics by motherhood group | `01_master.R` (descriptives) |
| `tab:couple_descriptives` | Within-couple labour-market outcomes | `01_master.R` (Part A) |
| `tab:hypotheses` | Empirical hypotheses and estimands | Set in the manuscript; not computed |
| `tab:housework_4grp` | Housework hours by gender and parenthood | `01_master.R` (Part B) |
| `tab:hw_within_fe` | Within-person hours change at parenthood | `01_master.R` (Part B) |
| `tab:specialisation` | Specialisation test, cross-partner elasticity | `01_master.R` (Part A) |
| `tab:event_study` | Child-penalty event study, women minus men | `02_cs_event_study.R` |
| `tab:iv_main` | FTB-B first stage and contemporaneous-wage IV | `04_part_g_ftbb_iv.R` |
| (training participation table) | Training around first birth, women minus men | `03_training_event_study.R` |
| `tab:within_couple` | Within-couple DML gender gap in complementarity | `01_master.R` (Part D) |
| `tab:main_beta2` | Education-experience complementarity by group | `01_master.R` (Parts D, E) |
| `tab:heckman` | Heckman selection correction | `01_master.R` (Part E) |
| `tab:transitioner_bound` | Within-person complementarity by transition type | `01_master.R` (Part E) |
| `tab:definition_sensitivity` | Motherhood-definition sensitivity | `01_master.R` (Part E) |
| `tab:robustness_master` | Master robustness summary | `01_master.R` (Part E) |
| `tab:ftbb_first_stage` | FTB-B first stages and implied elasticity | `04_part_g_ftbb_iv.R`, `05_part_h_policy.R` |

### Figures

| Figure file | Exhibit | Upstream dependency |
| --- | --- | --- |
| `fig2_event_study_beta2.pdf` | Complementarity event study | `event_study_beta2.rds` (Part E) |
| `fig3_employment_birth.pdf` | Employment around birth | `01_master.R` |
| `fig5C1_gender_gap_vanishes.pdf` | Within-couple gender gap | Part D objects |
| `fig7_robustness_forest.pdf` | Robustness forest | `robustness_forest_data.rds` (Part E) |
| `fig_event_study.pdf` | Child-penalty event study | `cs_event_study_results.rds` (stage 02) |
| `fig_training_event_study.pdf` | Training event study | `training_event_study_results.rds` (stage 03) |
| `fig_iv_coefplot.pdf` | FTB-B IV coefficients | `ftbb_reform1_iv_results.rds` (stage 04) |
| `fig_hw_within_fe.pdf` | Within-couple housework allocation | `01_master.R` (Part B) |

---

## 8. License

The code is released under the BSD-3-Clause license; see `LICENSE`. The license
covers this repository only and does not extend to the HILDA data, which remain
subject to the ADA confidentiality deed.

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
