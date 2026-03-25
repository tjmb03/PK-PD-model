# 🧪 PopPK / PKPD Workflow Shiny App

### An end-to-end population pharmacokinetic and pharmacodynamic modelling workflow built on the **nlmixr2 + rxode2** ecosystem, delivered as an interactive R Shiny application.


[![Launch App](https://img.shields.io/badge/Launch-Live%20Shiny%20App-blue?style=for-the-badge)](https://tjmb03.shinyapps.io/app26/)


------------------------------------------------------------------------

## Overview

This app provides a guided, modular workflow for fitting population PK and PK/PD models to clinical data — from raw CSV upload through parameter estimation, model comparison, diagnostics, and report export — without writing a single line of modelling code. It is designed around the philosophy: **simulate before you estimate**, then let the data refine your starting values.

---

## Features

### Data Handling
- Upload PK and PKPD datasets as CSV files
- Interactive column mapping — works with any column naming convention (e.g. `CONC`, `CP`, `Y` all map to `DV` automatically)
- Auto-detection of `EVID`, `MDV`, and `AMT` from data structure if columns are absent
- Support for multi-endpoint PKPD datasets via `DVID`

### Starting Value Estimation
Rather than guessing `ini()` parameters blind, the app offers two data-driven helpers:

**Guess PK from data (NCA)**
Runs non-compartmental analysis on the uploaded PK dataset to estimate:
- Typical **CL** from `Dose / AUC` (linear trapezoid)
- Typical **V** from `Dose / Cmax`
- Typical **KA / KTR** from `1.5 × ln(2) / Tmax`
- **Omega²** (BSV) from `var(log(individual estimates))`
- **Proportional residual error** from CV of Cmax across subjects

**Guess PD from data**
Estimates PD starting values from paired concentration-effect observations:
- **E0** from pre-dose baseline observations per subject
- **C50** by interpolating the concentration paired with ~50% of the maximum observed effect change
- **kout** from the effect recovery half-life after peak (turnover models)
- **Omega² E0 and C50** from log-variance across subjects
- Automatically detects **inhibitory vs stimulatory** effect direction
- **PD additive error** from within-subject baseline SD

All fields are optional — leave blank to use sensible built-in defaults. A "Starting Values" tab shows a live preview of exactly what will enter the `ini()` block.

### PK Model Library
| Model | Description |
|---|---|
| 1-cpt FO absorption (solved `linCmt`) | Fastest; uses analytical solution |
| 1-cpt FO absorption (ODE) | Explicit depot–central ODE system |
| 1-cpt + 1 transit compartment (ODE) | Captures modest absorption lag |
| 1-cpt + 5 transit compartments (ODE) | Full sigmoid absorption curve |

All models use log-normal parameterisation (`exp(lx + eta.x)`) and a combined proportional + additive residual error structure.

### PKPD Model Library
| Model | Description |
|---|---|
| Simultaneous Emax (Emax=1) with 1 transit PK | Direct-effect inhibitory/stimulatory model |
| Sequential PD-only (turnover, Emax=1) | Uses PK EBEs as fixed inputs; faster than simultaneous |

### Estimation
- **SAEM** (Stochastic Approximation EM) — robust global search; default for PK
- **FOCEi** (First-Order Conditional Estimation with Interaction) — faster; default for PD and PKPD
- Configurable burn-in / EM iterations for SAEM
- Full EBE (empirical Bayes estimates) extraction for sequential PKPD workflow

### Model Comparison
Stores every run automatically with a timestamped key. Comparison table shows:
- OBJF, AIC, BIC
- Parameter count (`Npar`)
- Mean ETA shrinkage
- Trophy 🏆 icon on the winner (lowest AIC)

### Diagnostics
- **ggPMX**: DV vs IPRED, NPDE vs Time
- **xpose.nlmixr2**: DV vs PRED
- **SAEM traceplot**: Convergence monitoring
- **VPC** (Visual Predictive Check): 300+ simulated trials with observed percentile overlay, linear and log-Y scale

### Export
One-click **self-contained HTML report** embedding all results:
- Data summary and preview
- Exploratory concentration-time spaghetti plot
- Model comparison table
- Fixed-effect parameter estimates
- GOF diagnostic plots
- VPC
- PKPD fit summary
- Starting values used in `ini()`

No internet connection required to view the exported file — all plots are base64-encoded inline.

---

## Workflow

```
Upload CSV  →  Map columns  →  Guess starting values (NCA)
     ↓
Run PK model (SAEM/FOCEi)  →  Compare models  →  VPC  →  Diagnostics
     ↓
[Optional] Extract EBEs  →  Merge into PD data  →  Run PD-only fit (sequential)
     ↓
[Optional] Run simultaneous PKPD model
     ↓
Export HTML report
```

A workflow status panel in the sidebar shows ✅ / ⬜ for each step so you always know where you are.

---

## Input Data Format

The app accepts standard NONMEM-style PK datasets. Required columns (names are flexible — the column mapper handles alternatives):

| Column | Standard name | Accepted alternatives |
|---|---|---|
| Subject ID | `ID` | `SUBJ`, `SUBJECT`, `USUBJID` |
| Time | `TIME` | `T`, `NOMTIME`, `NTIME` |
| Dependent variable | `DV` | `CONC`, `CP`, `Y`, `VALUE`, `OBS` |
| Dose amount | `AMT` | `DOSE`, `AMOUNT` |
| Event ID | `EVID` | auto-derived from AMT if absent |
| Missing DV flag | `MDV` | auto-derived from DV if absent |
| Compartment | `CMT` | optional |
| Endpoint ID | `DVID` | `ENDPOINT` — required for PKPD multi-endpoint data |

Optional covariates: `WT`, `AGE`, `SEX`, `SPARSE`

---

Both the variable-assignment approach (`eta.ka ~ eta_ka_sv`) and `bquote()`/`eval()` fail because nlmixr2 re-serialises the function body via `body()` → deparse → reparse before passing it to `lotri`, stripping any R-level substitutions. The `sprintf` approach is immune to this.

### BSV and Shrinkage Interpretation

High shrinkage (> ~30%) on an omega means the individual EBEs for that parameter should be interpreted with caution — the data per subject is insufficient to move individual estimates far from the population mean. High BSV paired with high shrinkage means the variability is real but unidentifiable from the available sampling scheme (typically a sparse early absorption phase for KA/KTR parameters with transit compartment models). Population-level estimates remain valid.

---

## Example Output

Running the warfarin PK dataset from the PopSim 2024 course yields the following model ranking:

| Model | OBJF | AIC | Winner |
|---|---|---|---|
| 1cpt + 5 transit (ODE) | 248.9 | 270.9 | 🏆 |
| 1cpt + 1 transit (ODE) | 341.1 | 363.1 | |
| 1cpt FO absorption (ODE) | 450.1 | 472.1 | |
| 1cpt FO absorption (solved) | 452.2 | 474.2 | |

Winner parameters: **KTR = 3.80 1/h** (BSV 50%, shrinkage 45%), **CL = 0.133 L/h** (BSV 28%, shrinkage 3%), **V = 7.99 L** (BSV 21%, shrinkage 12%). The large ΔOBJF gaps (>90 units) confirm a biologically meaningful absorption delay requiring multiple transit compartments to describe adequately.

---

## Dependencies

| Package | Role |
|---|---|
| `nlmixr2` | Population PK/PD parameter estimation (SAEM, FOCEi) |
| `rxode2` | ODE simulation engine; powers VPC |
| `vpc` | Visual Predictive Check plotting |
| `ggPMX` | GOF diagnostics (DV vs IPRED, NPDE) |
| `xpose.nlmixr2` | xpose-style diagnostic plots |
| `data.table` | Fast data manipulation |
| `ggplot2` + `xgxr` | Exploratory plots with log-scale support |
| `patchwork` | Plot composition |
| `htmltools` | HTML escaping for report export |
| `base64enc` | Plot embedding in HTML report |
| `nonmem2rx` *(optional)* | Import NONMEM `.ctl`/`.lst` runs |
| `babelmixr2` *(optional)* | Monolix interoperability |

---

## References

- Schoemaker R. *rxode2 and nlmixr2 — PopSim Course Copenhagen 2024*
- Fidler M et al. (2019). Nonlinear mixed-effects model development and simulation using nlmixr and related R open-source packages. *CPT Pharmacometrics Syst Pharmacol* 8(9):621–633.
- Keizer RJ et al. (2013). Modeling and Simulation Workbench for NONMEM: Tutorial on Pirana, PsN, and Xpose. *CPT Pharmacometrics Syst Pharmacol* 2:e50.

---


*Built with the nlmixr2 + rxode2 open-source pharmacometrics ecosystem.*
> © 2026 tjmb03. Source code for the interactive dashboards is **available on request** for academic and research use.
------------------------------------------------------------------------
