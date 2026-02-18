<div align="center">

```
 ██████╗██╗     ██╗███╗   ██╗    ██████╗ ██╗  ██╗ █████╗ ██████╗ ███╗   ███╗
██╔════╝██║     ██║████╗  ██║    ██╔══██╗██║  ██║██╔══██╗██╔══██╗████╗ ████║
██║     ██║     ██║██╔██╗ ██║    ██████╔╝███████║███████║██████╔╝██╔████╔██║
██║     ██║     ██║██║╚██╗██║    ██╔═══╝ ██╔══██║██╔══██║██╔══██╗██║╚██╔╝██║
╚██████╗███████╗██║██║ ╚████║    ██║     ██║  ██║██║  ██║██║  ██║██║ ╚═╝ ██║
 ╚═════╝╚══════╝╚═╝╚═╝  ╚═══╝   ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝
```

# Clinical Pharmacology Decision Engine

**A production-grade R Shiny application for end-to-end drug development decision support**

[![Launch App](https://img.shields.io/badge/▶%20Launch%20App-4A90D9?style=for-the-badge&logoColor=white)](https://your-shinyapps-link-here)
&nbsp;
[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r)](https://www.r-project.org/)
&nbsp;
[![Shiny](https://img.shields.io/badge/Shiny-1.7+-blue?style=flat-square)](https://shiny.posit.co/)
&nbsp;
[![rxode2](https://img.shields.io/badge/rxode2-ODE%20engine-orange?style=flat-square)](https://nlmixr2.github.io/rxode2/)
&nbsp;
[![License: MIT](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)

</div>

---

## What This Is

A fully-integrated PK/PD simulation and analysis platform covering the complete drug development workflow — from population PK simulation and covariate analysis through TMDD mechanistic modelling, Bayesian dose individualisation, survival analysis, and clinical trial power — all in a single Shiny interface with no NONMEM licence required.


---

## Modules

### PK/PD + Safety
Standard one-compartment PK simulation with Emax pharmacodynamics. Simulates multiple dose regimens simultaneously across user-defined subject populations, computing Cmax, Cmin, AUC, and probability of target attainment for each arm. Logistic toxicity modelling overlaid on the same simulation.

### TMDD Analysis — Target-Mediated Drug Disposition
Three mechanistic models selectable at runtime:
- **Full TMDD** — tracks free drug, free target, and drug-target complex with kon/koff/kint/ksyn/kdeg
- **QSS approximation** — quasi-steady-state simplification for fast equilibrating targets
- **Michaelis-Menten** — collapsed approximation for saturated binding

Outputs free/total drug concentration, target occupancy over time, and trough occupancy attainment rates.

### Indirect Response PD
Types I–IV indirect response models (inhibition/stimulation of production or loss) driven by the PK simulation output. Individual variability on kin/kout via log-normal IIV. Can be layered on top of either standard PK or TMDD — concentration column auto-detected.

### Covariate Analysis
Generates realistic virtual patient populations with allometric weight scaling, renal function (eGFR-based), hepatic impairment, sex effects, and CYP2D6 genotype. Per-subject CL/V adjustments propagated through all downstream simulations. Visualised as:
- Population characteristic histograms (age, weight, BMI, eGFR, sex, renal category)
- Exposure scatter plots: AUC vs weight, AUC vs eGFR, Cmax vs weight, AUC by sex with regression overlays

### Survival Analysis
Exposure-driven time-to-event simulation. AUC from the best-performing regimen drives a Weibull hazard model. Outputs:
- Kaplan-Meier curves stratified by high/low AUC exposure with at-risk counts
- Cox proportional hazards table (HR, 95% CI, p-value)

### Model Validation (VPC)
Visual Predictive Check against user-uploaded observed PK data:
- Auto-detects and renames common column name variants (ID/SUBJECTID, TIME/TAD, DV/CONC/cp)
- Strips dose records (MDV rows) before quantile computation
- Matches simulation time range to observed data automatically
- 50–500 simulation replicates, sequential or parallel (multi-core)
- Shaded 90% prediction interval ribbons with observed quantile overlays, individual data points, and automatic misfit annotation

### Therapeutic Drug Monitoring
Bayesian MAP estimation for individual dose adjustment. Takes observed concentration + time, updates log-normal priors on CL/V/KA, and recommends an adjusted dose to hit target AUC or Cmin with prediction plot.

### Bayesian PV Signal Detection
Beta-binomial pharmacovigilance signal detection with reporting odds ratio (with Haldane-Anscombe continuity correction), posterior probability of association, and Bayes factor.

### Postmarketing Signal
Adverse event reporting odds ratio calculator with continuity-corrected 95% CI.

### Clinical Trial Simulator
Binary endpoint trial simulation with sample size and power computation. Monte Carlo replication of trial outcomes.

### Trial + Power
Full power curves across sample sizes for binary and continuous endpoints. Monte Carlo with configurable simulation count.

---



## Key Technical Decisions

**rxode2 throughout** — all ODE systems (PK, TMDD, IDR) solved via rxode2 for speed and reproducibility. Event tables drive dosing schedules; covariate data frames passed via `covsDf` argument to avoid silent column-ignore bugs.

**Modular reactive graph** — each input panel is a Shiny module returning a reactive list. The central `observeEvent(run_btn)` collects all module outputs, runs the engine pipeline in sequence (covariate → PK → TMDD/IDR → metrics → decision), and stores results in `reactiveVal`s for all downstream plots and tables.

**Parallel-safe design** — `future::multisession` used throughout (never `multicore`, which is blocked inside Shiny forked processes). All `furrr` calls use `furrr_options(seed=TRUE)` for reproducibility. Plans are always restored with `on.exit()` to avoid leaking parallel state into subsequent reactive evaluations.

**Column normalisation layer** — TMDD outputs `Cfree`/`Ctotal`, standard PK outputs `cp`. All downstream consumers (IDR engine, VPC, metrics, plots) detect and alias whichever column is present rather than hardcoding `$cp`.

---



---

## Dependencies

```r
# Core
shiny, dplyr, ggplot2, readr, rlang

# PK/PD simulation
rxode2

# Survival analysis
survival

# Parallel computation
future, furrr, progressr

# Bayesian
rstan

# Data
jsonlite, tools
```

---


## Uploading Observed Data for VPC

The VPC tab accepts standard PK CSV files. Required columns (case-insensitive, common aliases auto-detected):

| Standard | Also accepted |
|----------|---------------|
| `id` | `ID`, `subject`, `SubjectID`, `patient` |
| `time` | `TIME`, `TAD`, `hours` |
| `DV` | `conc`, `concentration`, `cp`, `observed` |

Dose records (`DV=NA` or `DV=0` at time=0, standard in NONMEM format) are automatically excluded from observed quantile computation.

---

<div align="center">

Built with R · rxode2 · Shiny · ggplot2

*For research and educational use in clinical pharmacology and drug development.*

</div>

