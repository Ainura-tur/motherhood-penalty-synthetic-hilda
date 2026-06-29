# Replication Package

## Lost Experience and the Motherhood Wage Penalty: A Dynamic Account of Household Specialisation in Australia

Ainura Tursunalieva (CSIRO) and Ratbek Dzhumashev (Monash University)

This README follows the Social Science Data Editors template. It documents the
data, the code, the computing environment, and the steps an independent
researcher needs to reproduce every table and figure in the paper. The analysis
data are confidential and are not included in this package; the Data
availability and provenance section below explains how to obtain them. The
package also ships a synthetic HILDA panel so the entire pipeline can be run and
inspected without the confidential data, as described in Section 1.

---

## 1. Overview

The code in this repository constructs the analysis sample from the raw HILDA
release and reproduces all exhibits in the paper. The pipeline is written in R.
A single top-level script, `run_all.R`, executes the full sequence in order:
sample construction (or synthetic-panel generation), the main analysis (Parts A
to F), the estimator-aligned complementarity contrast (Part O), the
child-penalty event study, the implied dynamic cost and the Gelbach
experience-channel decomposition, the FTB-B reform instrumental-variables
analysis, the policy-economics extension, the training event study, and all
publication figures.

By default the package runs in **synthetic mode**: `run_all.R` builds synthetic
stand-in panels via `00b_make_synthetic_hilda.R`, so a fresh clone runs end to
end without the restricted HILDA data. Results on the synthetic panel are
illustrative and do **not** reproduce the paper's estimates. To reproduce the
published results, set `DATA_MODE=real` and build the panels from the restricted
release with `00_load_hilda.R` (Sections 2 and 6).

The reference run (real data) completed in approximately 4 hours 15 minutes on
the machine described in Section 4.

### 1a. Running the pipeline (real data and synthetic mode)

The pipeline runs in two modes, selected by the `DATA_MODE` setting in
`config.R` (or by the `DATA_MODE` environment variable, which overrides the
file without editing it). The default is synthetic, so a fresh clone runs end
to end with no data and no manual setup.

`DATA_MODE = "synthetic"` (default). The three analysis panels are built on the
fly by `00b_make_synthetic_hilda.R` from a fixed seed, and the full analysis
and all figures run against them. This requires no access to HILDA and reaches
every code path: the couple analysis, the housework mediation, the panel and
FTB-B instrumental-variables steps, the double machine learning and causal
forest, the robustness battery, the event studies, the implied-cost and Gelbach
decompositions, and the figures. Use this mode to verify that the code installs
and runs, to read the pipeline, or to develop against it.

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
then runs `01_master.R` (which sources the event studies, the implied-cost and
Gelbach decompositions, and the FTB-B IV) and `06_figures.R` in one session. In
real mode it stops with an instruction if the panels are missing, rather than
overwriting them.

#### Synthetic data are illustrative, not the paper's results

The synthetic panels are random data generated to match the schema of the real
analysis panels: the same variables, types, and factor levels, the same panel
and couple structure, and a built-in motherhood interruption and
education-by-experience wage gradient so the estimates carry sensible signs.
They exist so the code can be run and inspected without the confidential data.

Any number produced in synthetic mode is an artefact of the random generator.
These numbers are not estimates of the paper's quantities, will change with the
seed or sample size, and must not be reported or cited. Every result, table,
and figure in the manuscript comes from `DATA_MODE = "real"` on HILDA General
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
included. A synthetic stand-in panel is included so the pipeline can be executed
without the confidential data.

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
| Format expected | Stata `.dta` wave files (`Combined_<wave><release>0c.dta`) from the raw release |
| Location in package | Place the raw release in `data/raw/` or set `HILDA_DATA_DIR`; it is gitignored and not distributed |

The HILDA release is the only source dataset. No other external data are used.

---

## 3. Dataset list

In **real mode** the loader builds three derived analysis panels from the raw
release; in **synthetic mode** the same three panels are generated as synthetic
stand-ins. Either way they are written by bare name to the working directory (the
repository root) and are gitignored, so they are not distributed.

| File | Built by (real / synthetic) | Used by |
| --- | --- | --- |
| `hilda_panel_data_W12_W24_slim.rds` | `00_load_hilda.R` / `00b_make_synthetic_hilda.R` | Parts A to F (`01_master.R`), Part O, event study, implied cost, Gelbach |
| `hilda_panel_data_extended.rds` | `00_load_hilda.R` / `00b_make_synthetic_hilda.R` | Parts G and H |
| `hilda_panel_data_training.rds` | `00_load_hilda.R` / `00b_make_synthetic_hilda.R` | training event study |

---

## 4. Computational requirements

### 4.1 Software

The reference results were produced with R 4.5.2 (2025-10-31) on Windows 11 x64.
The full environment is recorded in `sessionInfo.txt`. Install dependencies with:

```r
source("PACKAGES.R")
```

`PACKAGES.R` installs and version-checks the CRAN packages used by the pipeline
and prints, for each, the installed version against the reference-run pin:
tidyverse, fixest, lmtest, sandwich, patchwork, DoubleML, mlr3, mlr3learners,
ranger, glmnet, grf, car, AER, MASS, future, and haven. `run_all.R` sources
`PACKAGES.R` as its first step, so a separate install call is not required when
running the full pipeline.

### 4.2 Controlled randomness

The random seed is set to 42 in `01_master.R`, before the double machine learning
and causal-forest steps, and is re-set locally before each stochastic draw. The
causal forest and the cross-fitting folds are the only sources of randomness in
the analysis. In synthetic mode, `00b_make_synthetic_hilda.R` uses its own fixed
seed so the synthetic panel is itself reproducible.

### 4.3 Memory and runtime

The real-data reference run completed in approximately 4 hours 15 minutes. Peak
memory was not separately recorded; a workstation with at least 16 GB of RAM is
recommended because the causal forest and double machine learning steps are
memory intensive. The synthetic-mode run is substantially faster. The thread
count is detected automatically: on a SLURM batch job the pipeline uses
`SLURM_CPUS_PER_TASK`, and interactively it caps at four threads. Override with
the `MASTER_HH_THREADS` environment variable.

---

## 5. Description of programs and code

Configuration and orchestration live at the repository root, and the numbered
scripts sit alongside them at the root (the scripts read and write by bare
filename, so they must run from the directory that contains them; `run_all.R`
sets that working directory automatically). The ordering matters: the Gelbach
decomposition reads the event-study results and the implied-cost output, and
Part H requires in-memory objects from Part G and the event-study results file,
so the event studies and the implied cost run before the decomposition and the
policy extension.

```
repo-root/
├── README.md                          this file
├── LICENSE                            BSD-3-Clause for code; CC-BY-4.0 for documentation
├── .gitignore
├── sessionInfo.txt                   environment record from the reference run
├── PACKAGES.R                        installs and version-checks dependencies
├── config.R                          DATA_MODE, data_dir, and output paths in one place
├── run_all.R                         master runner; sources the scripts below in order
├── 00_load_hilda.R                   raw HILDA waves to the three analysis panels (real mode)
├── 00b_make_synthetic_hilda.R        synthetic stand-in panels (synthetic mode, default)
├── 01_master.R                       Parts A to F: main couple analysis, housework
│                                     mediation, panel IV (diagnostic), couple DML and
│                                     causal forest, robustness battery, Oster bounds;
│                                     then sources the modules below in dependency order
├── 02b_part_o_aligned.R              estimator-aligned betaX contrast + Lee bounds (Table A4 / Table 13)
├── 02_cs_event_study.R               child-penalty event study (Sun and Abraham 2021)
├── 07b_implied_dynamic_cost_concave.R  late-horizon experience share (return x gap)
├── 08_gelbach_decomposition.R        experience-channel decomposition (Table 10 / Appendix A6)
├── 04_part_g_ftbb_iv.R               FTB-B 2015 reform simulated-instrument IV
├── 05_part_h_policy.R                policy economics of the FTB-B reform
├── 03_training_event_study.R         training event study around first birth
├── 06_figures.R                      all publication figures
├── data/
│   └── raw/                          empty; replicator places HILDA Release 24 here (real mode)
└── output/
    ├── tables/
    └── figures/
```

`run_all.R` sources `PACKAGES.R`, then `config.R`, then ensures the three panels
exist (building synthetic panels in the default mode), then sources `01_master.R`
and `06_figures.R` in the same R session. `01_master.R` internally sources the
analysis modules in this dependency order: `02b_part_o_aligned.R` (independent;
runs early), `02_cs_event_study.R`, `07b_implied_dynamic_cost_concave.R`,
`08_gelbach_decomposition.R`, `04_part_g_ftbb_iv.R`, `05_part_h_policy.R`, and
`03_training_event_study.R`. The 02 → 07b → 08 sequence is required because the
Gelbach output reads the implied-cost late-horizon share, and the policy
extension (05) reads Part G's (04) in-memory objects.

For reference, the correspondence to the authors' original working filenames is:
`01_master.R` was `MASTER_hh.R`; `02_cs_event_study.R` was `CS_event_study.R`;
`03_training_event_study.R` was `training_event_study.R`; `04_part_g_ftbb_iv.R`
was `MASTER_PART_G_FTBB_IV.R`; `05_part_h_policy.R` was `Master_part_h_policy.R`;
`06_figures.R` was `RUN_FIGURES_hh.R`; `02b_part_o_aligned.R` was
`Master_part_o.R`; `07b_implied_dynamic_cost_concave.R` was
`implied_dynamic_cost.R`; and `08_gelbach_decomposition.R` was
`gelbach_decomposition.R`.

---

## 6. Instructions to replicators

### Synthetic mode (no restricted data; default)

From the repository root:

```sh
Rscript run_all.R
```

This installs dependencies, builds the synthetic panels, runs the full analysis,
and writes all tables and figures to `output/`. Results are illustrative and do
not match the paper's estimates. At the end, `run_all.R` prints a
result-completeness check listing each expected result file as PASS or MISSING.

### Real mode (reproduces the published results)

1. Obtain HILDA General Release 24 from the Australian Data Archive, following
   Section 2. Place the raw wave files in `data/raw/`, or set `HILDA_DATA_DIR` to
   the folder that holds them.
2. Install dependencies with `source("PACKAGES.R")` (or let `run_all.R` do it).
3. Set any non-default paths via `config.R`, `HILDA_DATA_DIR`, or
   `MASTER_HH_THREADS`.
4. Run the full pipeline in real mode:

   ```sh
   DATA_MODE=real Rscript run_all.R
   ```

   `run_all.R` builds the analysis panels with `00_load_hilda.R`, runs the full
   analysis, and writes all tables and figures to `output/`. It runs in a single
   R session because the figure script reuses objects built earlier in the run.

To run a single stage, source `config.R` first, then the numbered script. The
analysis modules assume the panels from stage 00/00b are on disk;
`08_gelbach_decomposition.R` assumes `02_cs_event_study.R` and
`07b_implied_dynamic_cost_concave.R` have run; and `05_part_h_policy.R` assumes
`02_cs_event_study.R` and `04_part_g_ftbb_iv.R` have run in the same session.

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
| `tab:gelbach` | Gelbach experience-channel decomposition (Table 10) | `08_gelbach_decomposition.R` |
| `tab:iv_main` | FTB-B first stage and contemporaneous-wage IV | `04_part_g_ftbb_iv.R` |
| (training participation table) | Training around first birth, women minus men | `03_training_event_study.R` |
| `tab:within_couple` | Within-couple DML gender gap in complementarity | `01_master.R` (Part D) |
| `tab:main_beta2` | Education-experience complementarity by group | `01_master.R` (Parts D, E) |
| `tab:heckman` | Heckman selection correction | `01_master.R` (Part E) |
| `tab:transitioner_bound` | Within-person complementarity by transition type (Table A3) | `01_master.R` (Part E) |
| `tab:estimator_aligned` | Estimator-aligned gender contrast (Table A4) | `02b_part_o_aligned.R` |
| `tab:betaX_selection` | Within-couple betaX under selection corrections (Table 13) | `02b_part_o_aligned.R` |
| `tab:definition_sensitivity` | Motherhood-definition sensitivity | `01_master.R` (Part E) |
| `tab:robustness_master` | Master robustness summary | `01_master.R` (Part E) |
| `tab:ftbb_first_stage` | FTB-B first stages and implied elasticity | `04_part_g_ftbb_iv.R`, `05_part_h_policy.R` |

The reliable late-horizon (k = 5 to 10) experience share reported in Section V.D
and Appendix A6 (the return-times-gap figure) is computed by
`07b_implied_dynamic_cost_concave.R` and read into the Gelbach output by
`08_gelbach_decomposition.R`.

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

The code is released under the BSD-3-Clause license and the documentation under
CC-BY-4.0; see `LICENSE`. The license covers this repository only and does
not extend to the HILDA data, which remain subject to the ADA confidentiality
deed.

---

## 9. References

Gelbach, Jonah B. 2016. "When Do Covariates Matter? And Which Ones, and How
Much?" *Journal of Labor Economics* 34 (2): 509-543.

Lee, David S. 2009. "Training, Wages, and Sample Selection: Estimating Sharp
Bounds on Treatment Effects." *Review of Economic Studies* 76 (3): 1071-1102.

Oster, Emily. 2019. "Unobservable Selection and Coefficient Stability: Theory and
Evidence." *Journal of Business and Economic Statistics* 37 (2): 187-204.

Sun, Liyang, and Sarah Abraham. 2021. "Estimating Dynamic Treatment Effects in
Event Studies with Heterogeneous Treatment Effects." *Journal of Econometrics*
225 (2): 175-199.

Department of Social Services; Melbourne Institute of Applied Economic and Social
Research. 2024. *The Household, Income and Labour Dynamics in Australia (HILDA)
Survey, General Release 24.* ADA Dataverse. https://doi.org/10.26193/6M1BMR.
