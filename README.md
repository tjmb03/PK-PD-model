# PK/PD Modeling & Simulation

A collection of pharmacokinetic/pharmacodynamic (PK/PD) modeling tools, interactive Shiny applications, and example scripts for Model-Informed Drug Development (MIDD).

---

## Repository Structure

```
├── MIDD_advanced_shiny/          # Advanced MIDD interactive Shiny app
├── PopPKPD_workflow_shiny/       # Population PK/PD workflow Shiny app
├── rxode2_PKPDsimulator_shiny/   # RxODE2-based PK/PD simulator Shiny app
├── .github/workflows/            # CI/CD workflows
│
├── A simple indirect response PK:PD model.R
├── Convert RxODE2 QD:BID PK:PD simulation.R
├── PBPK_modeling_example.Rmd
├── TMDD_basic.Rmd
├── Steady-state PTA
└── basic example-nlmixr2.R
```

---

## Components

### 🖥️ Shiny Applications

| App | Description |
|-----|-------------|
| [MIDD_advanced_shiny](MIDD_advanced_shiny/) | Advanced model-informed drug development workflows |
| [PopPKPD_workflow_shiny](PopPKPD_workflow_shiny/) | Interactive population PK/PD analysis and simulation |
| [rxode2_PKPDsimulator_shiny](rxode2_PKPDsimulator_shiny/) | ODE-based PK/PD simulation using RxODE2 |

---

### 📄 Scripts & Notebooks

| File | Description |
|------|-------------|
| `A simple indirect response PK:PD model.R` | Indirect response model implementation in R |
| `Convert RxODE2 QD:BID PK:PD simulation.R` | Convert between once- and twice-daily dosing simulations |
| `PBPK_modeling_example.Rmd` | Physiologically-based PK (PBPK) modeling walkthrough |
| `TMDD_basic.Rmd` | Target-mediated drug disposition (TMDD) model example |
| `Steady-state PTA` | Probability of target attainment (PTA) at steady state |
| `basic example-nlmixr2.R` | Population PK fitting with nlmixr2 |

---

## Getting Started

See folders for scripts and usage details.

---

## Key Modelling Frameworks

| Framework | Purpose |
|-----------|---------|
| **RxODE2** | ODE-based PK/PD simulation in R |
| **nlmixr2** | Nonlinear mixed-effects population PK/PD fitting |
| **PBPK** | Physiologically-based pharmacokinetic modelling |
| **TMDD** | Target-mediated drug disposition modelling |
| **Indirect response** | PD models for delayed drug effects |
| **PTA** | Probability of target attainment for dose optimisation |

---

------------------------------------------------------------------------

Designed for translational scientists, pharmacometricians, and biomarker
strategists.
