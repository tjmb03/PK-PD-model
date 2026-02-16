# PK-PD Modeling Toolkit (R)

A focused portfolio of pharmacokinetic--pharmacodynamic (PK--PD)
modeling workflows built in R for translational decision-making.

## Capabilities

-   Population PK simulation with IIV
-   Indirect response (Imax / turnover) models
-   TMDD modeling
-   PBPK examples
-   Steady-state exposure metrics (Cmin, Cmax, AUC)
-   Probability of Target Attainment (PTA)
-   Interactive Shiny simulators

## Repository Structure

PopPKPD_workflow_shiny/ → Interactive PopPK workflow\
rxode2_PKPD_Shiny_simulator/ → PK/PD simulation app\
TMDD/ → Target-mediated drug disposition models

A simple indirect response PK:PD model.R\
basic example-nlmixr2.R\
PBPK_modeling_example.Rmd\
TMDD_basic.Rmd\
Steady-state PTA

## Example Workflow

Dose\
↓\
Population PK\
↓\
Exposure Metric\
↓\
Biomarker Model\
↓\
Probability of Effect\
↓\
PTA / Dose Decision

## Core Packages

rxode2 · nlmixr2 · dplyr · tibble · purrr · ggplot2 · shiny

## Modeling Philosophy

-   Mechanistic clarity\
-   Transparent assumptions\
-   Reproducible simulation\
-   Decision-oriented outputs\
-   No black-box modeling

------------------------------------------------------------------------

Designed for translational scientists, pharmacometricians, and biomarker
strategists.
