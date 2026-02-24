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


[![Launch App](https://img.shields.io/badge/▶%20Launch%20App-4A90D9?style=for-the-badge&logoColor=white)](https://tjmb03.shinyapps.io/MIDDadvanced/)
&nbsp;
[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r)](https://www.r-project.org/)
&nbsp;
[![Shiny](https://img.shields.io/badge/Shiny-1.7+-blue?style=flat-square)](https://shiny.posit.co/)
&nbsp;
[![rxode2](https://img.shields.io/badge/rxode2-ODE%20engine-orange?style=flat-square)](https://nlmixr2.github.io/rxode2/)

> A model-informed drug development (MIDD) platform built in R Shiny for forward simulation, multi-regimen optimisation, and clinical pharmacology decision support — from dose selection through post-marketing safety surveillance.

---

## Overview

This app implements the **forward simulation** half of the pharmacometrics workflow: given known or assumed PK/PD parameters, it simulates virtual patient populations across dosing scenarios and translates pharmacokinetic exposure into actionable clinical decisions. It is designed to answer the questions a clinical pharmacologist faces at each stage of drug development:

| Stage | Question answered |
|---|---|
| Phase I design | What dose range should we explore? |
| Phase II/III design | What dose achieves the target in ≥90% of patients? |
| Special populations | Does the dose need adjustment for renal impairment? |
| NDA/BLA submission | What are the exposure-response relationships? |
| Post-marketing | Is there a safety signal emerging? |
| TDM / clinical practice | What dose should this individual patient receive? |

> **Tip:** This app pairs naturally with a population PK estimation tool. Run estimation first to get CL, V, KA, and omega estimates, then use those values here as inputs for dose optimisation and trial design.

---

## Features

### Structural PK Models
Four compartment structures selectable from the sidebar — the app simulates the correct ODE for your drug, not just a default 1-compartment approximation:

| Model | Description | When to use |
|---|---|---|
| 1-compartment FO | Standard depot–central system | Most small molecules, rapid distribution |
| **2-compartment FO** | Adds peripheral tissue compartment | Biphasic decline, high-volume drugs, most biologics without TMDD |
| 1-cpt + 1 transit | One absorption transit compartment | Moderate lag to Cmax |
| 1-cpt + 5 transit | Five absorption transit compartments | Prolonged sigmoid rise, enteric-coated, complex formulations |
| **TMDD** | Full target-mediated disposition (full / QSS / Michaelis-Menten) | Monoclonal antibodies, high-affinity drugs |

Extra parameters (Q, V2 for 2-cpt; KTR for transit models) appear in the sidebar only when the relevant model is selected.

> **Why 2-compartment matters for dosing:** A 1-cpt model places all drug in the central compartment at trough. A 2-cpt model redistributes from peripheral tissue during the dosing interval, raising predicted Cmin — sometimes substantially. PTA calculations based on Cmin can be materially wrong if the wrong structural model is used.

---

### Dose Optimisation — Four Modes

The core of the app. Select a mode from the **Dosing Strategy** panel in the sidebar:

#### Manual
Simulate specific doses you enter (e.g. `50,100,200`). Computes PK/PD metrics and Pr(Success) for each. The original grid-search approach — use this when you have a short list of doses to evaluate.

#### Auto Sweep (PTA Curve)
Automatically sweeps a user-defined dose range (e.g. 10–500 mg in 10 mg steps) and plots three objectives simultaneously:
- **PTA** (probability of target attainment for the selected PK/PD index)
- **PD attainment** (Pr(Emax ≥ target))
- **Pr(Safe Cmax)** (logistic exposure–toxicity model)
- **Pr(Success)** = all three simultaneously

The recommended dose (maximum Pr(Success)) is highlighted, along with the full response surface so you can see whether you are in a flat region (forgiving) or a sharp peak (narrow therapeutic window requiring precision).

#### 2D Grid (Dose × Interval)
Co-optimises dose amount **and** dosing interval simultaneously. Produces:
- A **Pr(Success) heatmap** — darker = higher probability; winner marked with a star
- A **Cmin heatmap** — useful for trough-dependent drugs
- A best-dose-per-interval summary table

Use this when the dosing interval is a free design choice — the heatmap immediately shows whether a lower dose given more frequently outperforms a higher dose given less frequently for your specific drug.

#### Utility Maximisation
Uses a weighted utility function:

```
U = w_eff × PTA_PK  +  w_safe × Pr(Safe)  −  w_penalty × Pr(Toxic)
```

You control the three weights, making the tradeoff between efficacy and safety explicit and transparent. The utility curve is plotted alongside PTA curves so you can see exactly where the peak lies and how flat the optimum is.

---

### PK/PD Index Selection

Rather than a fixed Cmin target, the app supports the full pharmacodynamic index library:

| Index | Formula | Drug class |
|---|---|---|
| **Cmin** | Trough concentration | Beta-lactams, anticoagulants, many orals |
| **AUC** | Total area under curve | Vancomycin, most oral drugs |
| **Cmax/MIC** | Peak ÷ MIC | Aminoglycosides |
| **AUC/MIC** | AUC ÷ MIC | Fluoroquinolones, vancomycin (newer guidance) |
| **%T > threshold** | Fraction of interval above concentration | Time-dependent bacterial killing |
| **Cmax** | Peak concentration | Cytotoxics, peak-dependent PD |

MIC and threshold inputs appear conditionally when the relevant index is selected. The same index flows through manual simulation, PTA curves, the 2D grid, and subgroup analyses — everything is consistent.

---

### Subgroup PTA

Evaluates whether a single dose covers patients across the spectrum of renal function — critical for any renally eliminated drug. Four groups are simulated simultaneously:

| Group | eGFR | CL scaling |
|---|---|---|
| Normal | > 90 mL/min | Baseline |
| Mild RI | 60–90 mL/min | (75/90)^0.5 |
| Moderate RI | 30–60 mL/min | (45/90)^0.5 |
| Severe RI | < 30 mL/min | (20/90)^0.5 |

Bar charts show PTA and median Cmin per group with the 90% target line. If any group falls below 90%, the labelling may need a dose adjustment recommendation.

---

### Loading Dose

A checkbox in the Dosing Strategy panel enables a loading dose. Enter the loading/maintenance ratio (e.g. 2× for a double first dose). Pharmacokinetically, this front-loads exposure to hit therapeutic concentrations without waiting for accumulation over multiple half-lives — essential for long-half-life drugs (biologics, antivirals, anticoagulants).

---

### PD Modelling

| Model | Description |
|---|---|
| Direct Emax | Instantaneous Hill equation: E = Emax × C^n / (EC50^n + C^n) |
| **Indirect Response (IDR)** | All four IDR types: inhibit/stimulate kin or kout. Captures delayed response and turnover. With BSV on kin and kout. |
| **TMDD occupancy** | Target receptor occupancy from TMDD simulation |

---

### Covariate Analysis

Generate a virtual population with realistic covariate distributions (age, weight, BMI, eGFR, sex, optional genotype) and apply allometric/empirical covariate effects to CL and V. Effect magnitudes are fully user-configurable from the sidebar — no source code editing required:

- Weight on CL (power exponent, e.g. 0.75 for standard allometric scaling)
- eGFR on CL (power exponent)
- Weight on V (power exponent, typically 1.0)
- Sex on CL (female/male ratio)

Visualises: age/weight/BMI/eGFR histograms, AUC vs weight, AUC vs eGFR, AUC by sex.

---

### Therapeutic Drug Monitoring (TDM)

Bayesian dose individualisation from a single observed concentration:

1. Enter the dosing history (dose amount, dosing interval, number of prior doses — these are now wired to user inputs, not hardcoded)
2. Enter the observed concentration and the time since the last dose
3. Set AUC and Cmin targets
4. The app back-calculates individual CL, V, KA and recommends the adjusted dose

Plots current vs recommended dose predictions, with the observed concentration marked.

---

### Adaptive TDM Concept

As observations accumulate, each new concentration can be entered to progressively refine the individual PK estimate. The updated individual parameters drive a new forward simulation to recalculate the recommended dose — this is model predictive control applied to drug dosing.

---

### Survival Analysis

Exposure-driven time-to-event modelling from the simulated population:
- AUC per subject drives a parametric survival model (Weibull baseline hazard + beta_AUC log-hazard coefficient)
- Kaplan-Meier curves stratified by high vs low AUC exposure
- Cox proportional hazards model with hazard ratio, 95% CI, and p-value
- Quantifies the exposure-efficacy or exposure-toxicity relationship at the population level

---

### Bayesian Pharmacovigilance Signal Detection

Beta-Binomial model for post-marketing safety signal evaluation:
- Inputs: number exposed, AE reports in drug arm; comparator N and AE count
- Prior: Beta(a0, b0) — weakly informative by default
- Outputs: posterior mean RR, 95% credible interval, Pr(RR > threshold), signal/no-signal call
- 50,000 Monte Carlo draws for credible interval estimation

---

### Clinical Trial Design

**Trial Simulator** — simulate a binary endpoint trial (e.g. responder/non-responder) and read off observed response rates.

**Monte Carlo Power** — compute power for a binary outcome trial:
- Optionally use the simulated drug response rate directly from the PK/PD simulation
- Custom N grid for a power vs sample size curve
- Fisher's exact test or proportions test
- Produces a power curve plot with the 80% reference line

---

### VPC (Visual Predictive Check)

Upload observed PK data and compare it against simulations using the current parameter set. Handles:
- Auto-detection and renaming of non-standard column names
- Removal of dose records (MDV=1) before computing quantiles
- Parallel simulation using `future`/`furrr` for speed
- Time range extension to match the observed data range

---

### Postmarketing AE Signal Simulation

Simulate a post-marketing pharmacovigilance dataset under a specified background AE rate and drug relative risk, to understand power and time-to-signal characteristics.

---

### Export

**Download Run Package (ZIP)** — a reproducibility archive containing:
- `inputs.json` — all simulation parameters (doses, PK, PD, DDI, targets, safety thresholds)
- `subject_metrics.csv` — per-subject Cmax, Cmin, AUC, Emax, p_tox for the best regimen
- `regimen_summary.csv` — PTA table for all simulated regimens
- `sessionInfo.txt` — R and package version snapshot

A reviewer or colleague can reconstruct the exact simulation from the JSON inputs file.

---

## Workflow

```
Enter drug PK/PD parameters
    ↓
Select structural model (1-cpt / 2-cpt / transit / TMDD)
    ↓
Set dosing strategy (PK/PD index, targets, toxicity model)
    ↓
Choose optimisation mode:
  Manual → evaluate specific doses
  Sweep  → PTA curves over dose range
  2D Grid → dose × interval heatmap
  Utility → utility-weighted optimum
    ↓
Review Dose Optimisation tab:
  - PTA curves (efficacy + safety + success overlaid)
  - Pareto front (tradeoff surface)
  - Utility curve (recommended dose)
    ↓
Check Subgroup PTA tab → does the dose work in renal impairment?
    ↓
Refine with Covariate Analysis → weight/eGFR/sex effects on exposure
    ↓
Design the trial in Trial + Power tab
    ↓
Monitor patients with Bayesian TDM
    ↓
Watch for safety signals in Bayesian PV Signal tab
    ↓
Export Run Package (ZIP) for reproducibility / submission
```

---

## Input Parameters

### PK Parameters (via `mod_pk.R`)

| Parameter | Description | Default |
|---|---|---|
| CL | Typical population clearance (L/h) | 10 |
| V | Typical central volume (L) | 50 |
| KA | First-order absorption rate constant (1/h) | 1.2 |
| omega_CL | Between-subject variability on CL (SD of log-normal) | 0.30 |
| omega_V | Between-subject variability on V | 0.25 |
| ii | Dosing interval (h) | 24 |
| addl | Number of additional doses | 6 |
| n | Virtual subjects per simulation | 200 |

**2-compartment extras:** Q (inter-compartmental CL, L/h), V2 (peripheral volume, L), omega_V2

**Transit extras:** KTR (transit rate constant, 1/h), omega_KTR. Mean transit time = (n_compartments + 1) / KTR.

### PD Parameters (via `mod_pd.R`)

| Parameter | Description |
|---|---|
| Emax | Maximum drug effect (0–1) |
| EC50 | Concentration producing 50% Emax (mg/L) |
| Hill | Sigmoidicity factor |

### TMDD Parameters (via `mod_tmdd.R`)

kon, koff, ksyn, kdeg, kint, R0, target_occupancy, model_type (full / qss / mm)

### Indirect Response PD Parameters (via `mod_pd_idr.R`)

kin, kout, Imax, IC50, R0, IDR type (inhibit/stimulate kin or kout), omega_kin, omega_kout

### DDI Parameters (via `mod_ddi.R`)

Inhibitor concentration and Ki — computes CL fold-change as `1 + [I]/Ki` and scales CL accordingly.

---

## Dosing Strategy Panel

All optimisation controls are in the purple **Dosing Strategy** panel in the sidebar:

| Control | Purpose |
|---|---|
| Optimisation mode | Manual / Auto sweep / 2D grid / Utility |
| Min / Max / Step dose | Dose range for sweep and grid modes |
| Min / Max / Step interval | Interval range for 2D grid |
| Loading dose toggle | Adds a loading dose at a user-specified multiple of the maintenance dose |
| PK/PD index | Cmin, AUC, Cmax/MIC, AUC/MIC, %T>threshold, Cmax |
| MIC / threshold | Appear conditionally when relevant index is selected |
| Target index value | PTA = Pr(index ≥ this value) |
| Target Emax | PD attainment threshold |
| Logistic alpha / beta | Exposure–toxicity model parameters |
| Utility weights | w_eff, w_safe, w_penalty (utility mode only) |

---

## Technical Design Notes

### Reactive Architecture

All simulation results are stored in `reactiveVal` stores, cleanly separated from the UI:

| Store | Contents |
|---|---|
| `sim_r` | Full concentration-time data for all subjects and regimens |
| `met_r` | Per-subject metrics (Cmax, Cmin, AUC, Emax, p_tox) |
| `regsum_r` | Decision summary table (PTA, Pr_Success) per regimen |
| `chosen_reg_r` | Best regimen label |
| `opt_curve_r` | PTA curve data from sweep/utility optimisation |
| `grid2d_r` | Full dose × interval grid results |
| `subgroup_r` | Renal subgroup PTA results |
| `opt_winner_r` | Recommended dose, interval, and all probability metrics |

The main simulation (`run_btn_r`) and the optimisation sweep (`run_optimisation`) are decoupled — they write to different stores and render into different tabs. Running optimisation does not overwrite a manual simulation and vice versa.

### Simulation Cap

`compute_pta_curve()` caps the per-dose subject count at `min(n_subj, 200)` for sweep/grid modes and validates that the total simulation count stays ≤ 2000. This keeps sweep runs (50 doses × 200 subjects = 10,000 ODE solves) fast enough for interactive use. Manual simulation uses the full `n_subj` from the PK module.

### TMDD Column Aliasing

The TMDD engine returns `Cmax_free` / `Cmin_free` / `AUC_free` instead of `Cmax` / `Cmin` / `AUC`. These are aliased to generic names before `summarize_decision()` is called, so the decision layer is engine-agnostic. The PK plot similarly detects whichever concentration column is available (`cp`, `Cfree`, or `Ctotal`).

### VPC Simulation Matching

VPC simulations use the same structural model, dosing scheme, and PK parameters as the main simulation — including the structural model type and 2-cpt/transit extra parameters. A mismatch between the model used to generate observed data and the model used for VPC would invalidate the diagnostic.

---

## Relationship to Estimation Workflows

This app is a **forward simulation** tool — it projects outcomes from assumed parameters. It does not fit data or estimate parameters. The natural upstream tool is a population PK estimation platform (e.g. nlmixr2 / NONMEM / Monolix) which recovers the parameters from observed concentration data.

The recommended integrated workflow:

```
Phase I data
     ↓
PopPK estimation  →  CL = 8.4 L/h, V = 42 L, omega_CL = 0.31 ...
                              ↓
                   This app  →  Dose optimisation, PTA, subgroup, trial design
```

---

*Built on the rxode2 ODE simulation engine. Designed for model-informed drug development across the full clinical pharmacology workflow.*





