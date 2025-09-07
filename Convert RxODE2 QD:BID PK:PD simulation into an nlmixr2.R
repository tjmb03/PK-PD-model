#Define the population PK/PD model
library(rxode2)
library(nlmixr2)
library(dplyr)

pkpd_model <- function() {
  ini({
    # Fixed effects (log-transformed)
    tKA   <- log(0.294); eta.KA ~ 0.1
    tCL   <- log(18.6); eta.CL ~ 0.1
    tV2   <- log(40.2); eta.V2 ~ 0.1
    Q     <- 10.5
    V3    <- 297
    Kin   <- 1
    Kout  <- 1
    EC50  <- 200
    add.sd <- 0.5
  })
  
  model({
    # Individual parameters
    KA <- exp(tKA + eta.KA)
    CL <- exp(tCL + eta.CL)
    V2 <- exp(tV2 + eta.V2)
    
    # Concentrations
    C2 <- centr / V2
    C3 <- peri / V3
    
    # ODEs
    d/dt(depot)  <- -KA * depot
    d/dt(centr)  <- KA * depot - CL*C2 - Q*C2 + Q*C3
    d/dt(peri)   <- Q*C2 - Q*C3
    d/dt(eff)    <- Kin - Kout*(1 - C2/(EC50 + C2)) * eff
    
    # Observations
    cp <- C2
    cp ~ add(add.sd)
  })
}

#Create QD/BID dosing schedules

# QD dosing rows
qd_dose <- data.frame(
  id   = 1,
  time = seq(0, by = 24, length.out = 5),  # dosing at 0,24,48,72,96
  amt  = 10000,
  cmt  = 1,
  evid = 1,   # dosing record
  dv   = NA
)

# QD observation rows
qd_obs <- data.frame(
  id   = 1,
  time = c(0:24, seq(28, 96, by = 4)),  # hourly first 24h, then q4h
  amt  = 0,
  cmt  = 2,   # central compartment (observed concentration)
  evid = 0,   # observation
  dv   = NA
)

# BID dosing rows
bid_dose <- data.frame(
  id   = 2,
  time = seq(0, by = 12, length.out = 10),  # dosing every 12h
  amt  = 10000,
  cmt  = 1,
  evid = 1,
  dv   = NA
)

# BID observation rows
bid_obs <- data.frame(
  id   = 2,
  time = c(0:24, 96 + 0:24),  # first 24h, then again at 96–120h
  amt  = 0,
  cmt  = 2,
  evid = 0,
  dv   = NA
)

pkpd_data <- dplyr::bind_rows(qd_dose, qd_obs, bid_dose, bid_obs) %>%
  dplyr::arrange(id, time, desc(evid))

head(pkpd_data)

fit <- nlmixr2(pkpd_model, data = pkpd_data,
               est = "saem", saemControl(print=0))
print(fit)
