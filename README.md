# PK–PD Modeling Toolkit (R / RxODE2 / nlmixr2)

A modular pharmacokinetic–pharmacodynamic (PK–PD) modeling framework in R.

This repository provides reproducible workflows for:

- Population PK simulation  
- Mechanistic PK/PD modeling  
- Target-mediated drug disposition (TMDD)  
- PBPK modeling  
- Exposure–response analysis  
- Steady-state Probability of Target Attainment (PTA)  

Structured as a professional pharmacometrics toolkit rather than standalone scripts.

---

# Repository Structure

<img width="443" height="263" alt="Screenshot 2026-02-14 at 12 05 17 PM" src="https://github.com/user-attachments/assets/72fc21b2-a98a-47d8-a8fd-b3264c89e3d2" />


---

# Core Modules

## 1️⃣ Population Simulation (PopSim2024/)

Monte Carlo simulation framework for population PK.

Focus:
- Inter-individual variability (IIV)  
- Exposure metrics (Cmin, Cmax, AUC)  
- Steady-state simulation  
- Dose comparison  

---

## 2️⃣ Indirect Response PK–PD

File:
`A simple indirect response PK:PD model.R`

Implements:
- Inhibitory/stimulatory indirect response models  
- Exposure–biomarker linkage  
- Time-course simulation  

---

## 3️⃣ Target-Mediated Drug Disposition (TMDD/)

Mechanistic nonlinear PK modeling.

Focus:
- Binding kinetics  
- Nonlinear clearance  
- Sensitivity exploration  
- Template expansion for biologics  

---

## 4️⃣ PBPK Modeling

Files:
- `PBPK_modeling_example.Rmd`
- `PBPK_modelling_minimal.Rmd`

Implements:
- Physiological compartment modeling  
- Structured parameterization  
- Simulation and visualization  

---

## 5️⃣ Nonlinear Mixed Effects Modeling

File:
`basic example-nlmixr2.R`

Demonstrates:
- Model specification  
- Random effects  
- Parameter estimation  
- Diagnostics  

---

## 6️⃣ Steady-State Probability of Target Attainment (PTA)

File:
`Steady-state PTA`

Implements a decision-support workflow:

1. Simulate population PK  
2. Extract steady-state exposure  
3. Apply PD target threshold  
4. Compute PTA across doses  

Designed for dose optimization and development-stage risk assessment.

---

# Capabilities Demonstrated

- Population PK simulation  
- Mechanistic PK/PD modeling  
- TMDD nonlinear dynamics  
- PBPK structural modeling  
- Mixed effects modeling  
- Exposure–response integration  
- Decision-oriented PTA computation  

---

# Technical Stack

- R ≥ 4.2  
- rxode2  
- nlmixr2  
- dplyr  
- ggplot2  
- rmarkdown  

Install core packages:

```r
install.packages(c(
  "rxode2",
  "nlmixr2",
  "dplyr",
  "ggplot2",
  "rmarkdown"
))
