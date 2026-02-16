# 🧪 PopPK / PKPD Workflow Shiny App

### End-to-End Interactive Population Modeling Demo (nlmixr2 + rxode2)

[![Launch App](https://img.shields.io/badge/Launch-Live%20Shiny%20App-blue?style=for-the-badge)](https://tjmb03.shinyapps.io/app26/)


------------------------------------------------------------------------

## 🚀 Overview

A workflow-oriented Shiny application implementing a structured
**Population PK / PKPD modeling pipeline**:

-   Data exploration\
-   Base PK model building\
-   Diagnostics & GOF\
-   VPC generation\
-   Model comparison\
-   Covariate scaffolding\
-   Sequential PK → PD\
-   Simultaneous PKPD

Built on the modern **nlmixr2 + rxode2 ecosystem**.

------------------------------------------------------------------------

## 🧠 Modeling Workflow

``` mermaid
flowchart LR
    A[Upload PK Data] --> B[Explore & QC]
    B --> C[Base PK Model]
    C --> D[Diagnostics]
    D --> E[VPC]
    C --> F[Model Comparison]
    F --> G[Covariate Scaffold]

    C --> H[Extract EBEs]
    H --> I[Sequential PD Model]

    C --> J[Simultaneous PKPD Model]
```

------------------------------------------------------------------------

## 🔬 Implemented Features

### ✅ PK Structural Models

-   1-compartment (linCmt)
-   1-compartment ODE
-   1-transit absorption
-   5-transit absorption

### ✅ Estimation Methods

-   SAEM
-   FOCEi (v5.0.0 compatible)

### ✅ Diagnostics

-   DV vs IPRED
-   NPDE vs TIME
-   ggPMX integration
-   Optional xpose support

### ✅ Model Comparison Panel

  Metric      Description
  ----------- --------------------------------
  OBJF        Objective function value
  Npar        Parameter count
  Shrinkage   Mean ETA shrinkage
  AIC         Akaike Information Criterion
  BIC         Bayesian Information Criterion
  🏆 Winner   Best model flag

------------------------------------------------------------------------

### ✅ VPC Module

-   Custom binning
-   Prediction intervals
-   Observed median overlay
-   Log-scale support

------------------------------------------------------------------------

### ✅ Sequential PK → PD

Workflow:

    PK Fit → Extract EBEs → Merge into PD → PD-only Turnover Model

Includes:

-   Character-safe ID merging
-   EBE validation
-   FOCEi stabilization (bobyqa + central gradients)
-   ODE tolerance control (rxControl)

------------------------------------------------------------------------

### ✅ Simultaneous PKPD

Implements:

-   Immediate Emax (Emax=1)
-   Transit PK
-   Multi-endpoint residual structure
-   PKPD VPC

------------------------------------------------------------------------

## 🏗 Architecture

### Tech Stack

-   nlmixr2 (v5.0.0)
-   rxode2
-   ggPMX
-   xpose.nlmixr2
-   vpc
-   xgxr
-   data.table
-   shiny


------------------------------------------------------------------------

## 🎯 Intended Use Cases

-   Transit absorption screening
-   Structural PK model selection
-   Exposure--response linking
-   Biomarker modeling
-   Sequential vs simultaneous comparison
-   Teaching NLME modeling concepts

------------------------------------------------------------------------

## 👤 Author

**tjmb03**\
Population PK/PD & Translational Modeling\
Boston, MA

------------------------------------------------------------------------

⭐ It is intended as a portfolio demonstration of pharmacometric workflow design, not a source-code distribution repository. Full source code available upon request.
