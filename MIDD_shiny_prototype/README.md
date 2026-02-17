# Clinical Pharmacology Decision Engine

A Model-Informed Drug Development (MIDD) prototype for PopPK/PD simulation, safety evaluation, pharmacovigilance signal updating, and clinical trial power estimation.

🔗 **Launch the Application**
👉 [Click here to open the live app](https://tjmb03.shinyapps.io/MIDDengine/)

---

## Overview

The Clinical Pharmacology Decision Engine integrates:

- Population PK simulation (rxode2-based)
- Exposure–response (Emax / sigmoid Emax)
- Safety threshold evaluation (Cmax-based)
- Drug–drug interaction modeling (CYP inhibition framework)
- Bayesian pharmacovigilance signal updating
- Trial simulation + Monte Carlo power estimation
- FDA-style Clinical Pharmacology Summary generation
- Reproducible Run Package export (ZIP)

This tool demonstrates how quantitative pharmacology can be connected end-to-end:
Dose → Exposure → Effect → Safety → Signal → Trial Design

---

## Core Capabilities

### 1. Population PK Simulation
- One-compartment model with first-order absorption
- Log-normal interindividual variability
- DDI inhibition effect on clearance
- Multi-regimen comparison

### 2. Pharmacodynamic Modeling
- Emax / sigmoid Emax
- Target attainment probability
- Composite success scoring (PK + PD + Safety)

### 3. Safety Layer
- Toxicity threshold evaluation
- Exposure–toxicity logistic modeling
- Regimen ranking based on composite success probability

### 4. Bayesian Signal Detection
- Posterior updating for AE rates
- Risk ratio threshold evaluation
- Signal probability estimation

### 5. Trial Simulation
- Binary endpoint simulation
- Monte Carlo power estimation
- Optional integration of simulated PK/PD response probability

### 6. Reproducibility Package
Each run generates a downloadable ZIP containing:
- Input parameters
- Subject-level metrics
- Regimen summary table
- Session information

---

## Intended Audience

- Clinical Pharmacologists
- Pharmacometrics Scientists
- Translational Medicine Teams
- MIDD / Regulatory Strategy Professionals
- Drug Development Decision Makers

---

## Technical Architecture

- R + Shiny front-end
- rxode2 simulation engine
- Modular engine design
- Monte Carlo simulation framework
- Bayesian posterior updating
- Reproducibility-first output structure

Source code is proprietary and not publicly distributed.

---

## Disclaimer

This application is a research and portfolio demonstration tool.
It is not intended for clinical or regulatory decision-making without proper validation and review.

---

## Contact

For collaboration, demonstration access, or technical discussion:
📩 [tjmb03@gmail.com]

