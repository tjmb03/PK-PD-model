# rxode2 PK/PD Shiny Simulator

Interactive PK/PD simulation tool built with **rxode2 + Shiny**.

## 🔬 Features

- Switch ODE systems (base / transit / effect / custom)
- Modify PK parameters (ka, CL, V, ktr, ke0, EC50, Hill…)
- Bolus + repeated infusion block
- Multi-scenario overlay:
  - Dose levels
  - Parameter sweep grid
- Automatic summary metrics:
  - Cmax, Tmax, AUC₀–t
  - Cmaxτ, Cminτ, AUCτ (τ = 12h / 24h / II)
  - Emaxτ (if PD model)
- Export:
  - Time-series CSV
  - Summary CSV

---

## 🚀 Live Demo

👉 **[Launch App](https://tjmb03.shinyapps.io/popsim/)**

---

## 📸 Preview

![App Preview](assets/<img width="2274" height="678" alt="Screenshot 2026-02-15 at 9 26 55 PM" src="https://github.com/user-attachments/assets/25ace38e-3f19-49af-9cdb-5e09b317663e" />
preview.png)

---

## 💻 Run Locally

⭐ It is intended as a portfolio demonstration of pharmacometric workflow design, not a source-code distribution repository. Full source code available upon request.

