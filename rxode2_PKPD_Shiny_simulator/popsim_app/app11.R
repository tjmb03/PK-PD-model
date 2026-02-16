# app.R
library(shiny)
library(rxode2)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---------- helpers ----------
make_base_ode <- function() {
  "
  d/dt(depot)   = -ka*depot;
  d/dt(central) =  ka*depot - (cl/v)*central;
  C1 = central/v;
  "
}

make_1trans_ode <- function() {
  "
  d/dt(depot)   = -ka*depot;
  d/dt(trans1)  =  ka*depot - ktr*trans1;
  d/dt(central) =  ktr*trans1 - (cl/v)*central;
  C1 = central/v;
  "
}

make_ntrans_ode <- function(n = 5) {
  stopifnot(n >= 1)
  lines <- c("d/dt(depot) = -ka*depot;")
  lines <- c(lines, "d/dt(trans1) = ka*depot - ktr*trans1;")
  if (n > 1) {
    for (i in 2:n) {
      lines <- c(lines, sprintf("d/dt(trans%d) = ktr*trans%d - ktr*trans%d;", i, i - 1, i))
    }
  }
  lines <- c(lines, sprintf("d/dt(central) = ktr*trans%d - (cl/v)*central;", n))
  lines <- c(lines, "C1 = central/v;")
  paste(lines, collapse = "\n")
}

make_effect_ode <- function(base = c("base", "1trans", "ntrans"), ntrans = 5) {
  base <- match.arg(base)
  ode <- switch(
    base,
    base   = make_base_ode(),
    `1trans` = make_1trans_ode(),
    ntrans = make_ntrans_ode(ntrans)
  )
  paste0(
    ode,
    "\n",
    "
    d/dt(Ce) = ke0*(C1 - Ce);
    E = E0 + (Emax * (Ce^Hill)) / (EC50^Hill + Ce^Hill);
    "
  )
}

compile_model_safe <- function(ode_text) {
  tryCatch(rxode2::rxode2(ode_text), error = function(e) e)
}

build_event_table <- function(input, bolus_amt_override = NULL) {
  ev <- rxode2::eventTable()
  ev$add.sampling(seq(0, input$t_end, input$t_step))
  
  if (isTRUE(input$bolus_on)) {
    bol_amt <- if (is.null(bolus_amt_override)) input$bolus_amt else bolus_amt_override
    ev$add.dosing(
      dose = bol_amt,
      start.time = input$bolus_time,
      dosing.to = input$bolus_cmt
    )
  }
  
  if (isTRUE(input$inf_on)) {
    ev$add.dosing(
      dose = input$inf_amt,
      nbr.doses = input$inf_n,
      dosing.to = input$inf_cmt,
      dosing.interval = input$inf_ii,
      rate = input$inf_rate,
      start.time = input$inf_start
    )
  }
  
  ev
}

# Convert comma-separated or range string to numeric vector
parse_num_vec <- function(x) {
  x <- trimws(x)
  if (x == "") return(numeric(0))
  
  # allow "0, 100, 250, 500"
  if (grepl(",", x, fixed = TRUE)) {
    vals <- strsplit(x, ",", fixed = TRUE)[[1]]
    vals <- suppressWarnings(as.numeric(trimws(vals)))
    return(vals[!is.na(vals)])
  }
  
  # allow "start:step:end" e.g. "0:50:500"
  if (grepl(":", x, fixed = TRUE)) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    parts <- suppressWarnings(as.numeric(trimws(parts)))
    if (length(parts) == 3 && all(!is.na(parts))) {
      return(seq(parts[1], parts[3], by = parts[2]))
    }
  }
  
  # single number
  v <- suppressWarnings(as.numeric(x))
  if (!is.na(v)) return(v)
  numeric(0)
}

summarize_scenarios <- function(df, tau = 24) {
  df <- df %>% arrange(scenario, time)
  
  has_c1 <- "C1" %in% names(df)
  has_e  <- "E"  %in% names(df)
  
  # ---- full-horizon summary ----
  full <- df %>%
    group_by(scenario) %>%
    summarize(
      n_time = n(),
      t_end  = max(time, na.rm = TRUE),
      
      Cmax = if (has_c1) max(C1, na.rm = TRUE) else NA_real_,
      Tmax = if (has_c1) time[which.max(C1)][1] else NA_real_,
      AUC_0_t = if (has_c1) {
        sum(diff(time) * (head(C1, -1) + tail(C1, -1)) / 2, na.rm = TRUE)
      } else NA_real_,
      Ct_last = if (has_c1) C1[which.max(time)][1] else NA_real_,
      Cmin = if (has_c1) min(C1, na.rm = TRUE) else NA_real_,
      
      Emax = if (has_e) max(E, na.rm = TRUE) else NA_real_,
      t_Emax = if (has_e) time[which.max(E)][1] else NA_real_,
      Emin = if (has_e) min(E, na.rm = TRUE) else NA_real_,
      
      .groups = "drop"
    )
  
  # ---- last-interval (tau) summary ----
  # Compute over [t_end - tau, t_end] per scenario
  last_tau <- df %>%
    group_by(scenario) %>%
    mutate(
      t_end = max(time, na.rm = TRUE),
      t0 = t_end - tau
    ) %>%
    ungroup() %>%
    filter(time >= t0) %>%
    group_by(scenario) %>%
    summarize(
      tau_used = tau,
      t_start_tau = min(time, na.rm = TRUE),
      t_end_tau   = max(time, na.rm = TRUE),
      n_time_tau  = n(),
      
      Cmax_tau = if (has_c1) max(C1, na.rm = TRUE) else NA_real_,
      Tmax_tau = if (has_c1) time[which.max(C1)][1] else NA_real_,
      Cmin_tau = if (has_c1) min(C1, na.rm = TRUE) else NA_real_,
      AUC_tau  = if (has_c1) {
        sum(diff(time) * (head(C1, -1) + tail(C1, -1)) / 2, na.rm = TRUE)
      } else NA_real_,
      Cavg_tau = if (has_c1) {
        # Use actual covered duration in case grid/sampling is shorter than tau
        dur <- max(time, na.rm = TRUE) - min(time, na.rm = TRUE)
        if (isTRUE(dur > 0)) AUC_tau / dur else NA_real_
      } else NA_real_,
      
      Emax_tau = if (has_e) max(E, na.rm = TRUE) else NA_real_,
      t_Emax_tau = if (has_e) time[which.max(E)][1] else NA_real_,
      Emin_tau = if (has_e) min(E, na.rm = TRUE) else NA_real_,
      
      .groups = "drop"
    )
  
  # Join full + tau summaries
  out <- full %>%
    left_join(last_tau, by = "scenario")
  
  out
}



# ---------- UI ----------
ui <- fluidPage(
  titlePanel("rxode2 PK/PD Simulator (single + scenario overlay + CSV export)"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      h4("1) Choose model"),
      radioButtons(
        "model_type", "ODE system",
        choices = c(
          "Base (Depot -> Central)" = "base",
          "1-Transit (Depot -> Trans1 -> Central)" = "1trans",
          "N-Transit (Depot -> Trans1..N -> Central)" = "ntrans",
          "Base + Effect compartment" = "base_eff",
          "1-Transit + Effect compartment" = "1trans_eff",
          "N-Transit + Effect compartment" = "ntrans_eff",
          "Custom ODE (edit below)" = "custom"
        ),
        selected = "base"
      ),
      
      conditionalPanel(
        condition = "input.model_type == 'ntrans' || input.model_type == 'ntrans_eff'",
        sliderInput("n_trans", "Number of transit compartments (N)", min = 1, max = 10, value = 5, step = 1)
      ),
      
      conditionalPanel(
        condition = "input.model_type == 'custom'",
        tags$div(
          style = "margin-top:10px;",
          textAreaInput("ode_custom", "Paste/edit ODE text here", rows = 14, width = "100%",
                        value = make_base_ode()),
          helpText("Tip: define C1=central/v if you want concentration plots."),
          actionButton("compile_btn", "Compile custom ODE", class = "btn-primary")
        )
      ),
      
      hr(),
      h4("2) Parameters"),
      sliderInput("ka", "ka (1/h)", min = 0.01, max = 5, value = log(2)/0.5, step = 0.01),
      sliderInput("cl", "CL (L/h)", min = 0.01, max = 2, value = 0.135, step = 0.005),
      sliderInput("v",  "V (L)",    min = 1, max = 100, value = 8, step = 0.5),
      
      conditionalPanel(
        condition = "input.model_type == '1trans' || input.model_type == 'ntrans' || input.model_type == '1trans_eff' || input.model_type == 'ntrans_eff'",
        sliderInput("ktr", "ktr (1/h)", min = 0.01, max = 2, value = log(2)/5, step = 0.01)
      ),
      
      conditionalPanel(
        condition = "input.model_type == 'base_eff' || input.model_type == '1trans_eff' || input.model_type == 'ntrans_eff'",
        tags$div(
          h5("Effect compartment / PD"),
          sliderInput("ke0",  "ke0 (1/h)", min = 0.001, max = 2, value = 0.15, step = 0.001),
          numericInput("E0",   "E0 (baseline effect)", value = 0),
          numericInput("Emax", "Emax", value = 1),
          numericInput("EC50", "EC50 (same units as Ce)", value = 0.5),
          sliderInput("Hill",  "Hill", min = 0.5, max = 5, value = 1, step = 0.1)
        )
      ),
      
      hr(),
      h4("3) Dosing"),
      checkboxInput("bolus_on", "Include bolus", value = TRUE),
      fluidRow(
        column(6, numericInput("bolus_amt", "Bolus amt (mg)", value = 500)),
        column(6, numericInput("bolus_time", "Bolus time (h)", value = 0))
      ),
      numericInput("bolus_cmt", "Bolus compartment (cmt)", value = 1, min = 1, step = 1),
      
      checkboxInput("inf_on", "Include repeated infusion block", value = TRUE),
      fluidRow(
        column(6, numericInput("inf_amt", "Inf amt (mg)", value = 250)),
        column(6, numericInput("inf_rate", "Rate (mg/h)", value = 125))
      ),
      fluidRow(
        column(4, numericInput("inf_n", "Nbr doses", value = 3, min = 1)),
        column(4, numericInput("inf_ii", "II (h)", value = 12, min = 0)),
        column(4, numericInput("inf_start", "Start (h)", value = 36))
      ),
      numericInput("inf_cmt", "Infusion compartment (cmt)", value = 2, min = 1, step = 1),
      
      hr(),
      h4("4) Sampling grid"),
      fluidRow(
        column(6, numericInput("t_end", "End time (h)", value = 120, min = 0.1)),
        column(6, numericInput("t_step", "Step (h)", value = 0.1, min = 0.001))
      ),
      
      hr(),
      actionButton("run_btn", "Run single scenario", class = "btn-success"),
      br(), br(),
      verbatimTextOutput("status"),
      
      hr(),
      h4("5) Scenario overlay"),
      radioButtons(
        "scenario_mode", "Overlay mode",
        choices = c("Dose levels (override bolus dose)" = "dose",
                    "Parameter sweep grid" = "sweep"),
        selected = "dose"
      ),
      
      hr(),
      h4("Steady-state interval metrics"),
      radioButtons(
        "tau_mode", "Tau (τ) for last-interval metrics",
        choices = c("Use infusion II (if infusion ON)" = "use_ii",
                    "12 hours" = "12",
                    "24 hours" = "24"),
        selected = "use_ii"
      ),
      helpText("Last-interval metrics are computed over [t_end-τ, t_end]. If 'Use infusion II' but infusion is OFF, τ falls back to 24h."),
      
      
      conditionalPanel(
        condition = "input.scenario_mode == 'dose'",
        tags$div(
          helpText("Enter doses as: 0,100,250,500  OR  0:50:500 (start:step:end)."),
          textInput("dose_list", "Bolus dose levels (mg)", value = "0, 100, 250, 500"),
          actionButton("run_scen_btn", "Run overlay", class = "btn-info")
        )
      ),
      
      conditionalPanel(
        condition = "input.scenario_mode == 'sweep'",
        tags$div(
          helpText("Choose 1 or 2 parameters to sweep. Each uses values like: 0.1,0.2,0.5  or  0.1:0.1:1.0"),
          selectInput("sweep_p1", "Parameter 1", choices = c("ka","cl","v","ktr","ke0","EC50","Emax","Hill"), selected = "ka"),
          textInput("sweep_v1", "Values for Parameter 1", value = "0.5, 1, 2"),
          selectInput("sweep_p2", "Parameter 2 (optional)", choices = c("None","ka","cl","v","ktr","ke0","EC50","Emax","Hill"), selected = "None"),
          textInput("sweep_v2", "Values for Parameter 2 (optional)", value = ""),
          actionButton("run_scen_btn", "Run overlay", class = "btn-info")
        )
      )
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Single: Concentration / Effect",
                 plotOutput("pk_plot", height = 320),
                 conditionalPanel(
                   condition = "input.model_type == 'base_eff' || input.model_type == '1trans_eff' || input.model_type == 'ntrans_eff'",
                   plotOutput("pd_plot", height = 320)
                 )
        ),
        tabPanel("Single: Compartments",
                 plotOutput("amt_plot", height = 420)
        ),
        tabPanel("Scenarios (overlay)",
                 plotOutput("overlay_pk_plot", height = 330),
                 plotOutput("overlay_pd_plot", height = 330),
                 
                 fluidRow(
                   column(6, downloadButton("download_csv", "Download time-series CSV")),
                   column(6, downloadButton("download_summary_csv", "Download summary CSV"))
                 ),
                 br(),
                 
                 h4("Scenario summary (per curve)"),
                 tableOutput("overlay_summary_tbl"),
                 br(),
                 
                 h4("Time-series preview (head)"),
                 tableOutput("overlay_head")
        ),
        
        tabPanel("Data (single head)",
                 tableOutput("head_tbl")
        ),
        tabPanel("Current ODE",
                 verbatimTextOutput("ode_show")
        )
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  tau_ss <- reactive({
    if (input$tau_mode == "use_ii") {
      if (isTRUE(input$inf_on) && !is.null(input$inf_ii) && is.finite(input$inf_ii) && input$inf_ii > 0) {
        return(as.numeric(input$inf_ii))
      } else {
        return(24) # fallback
      }
    }
    as.numeric(input$tau_mode) # "12" or "24"
  })
  
  
  custom_mod <- reactiveVal(NULL)
  custom_msg <- reactiveVal("Custom ODE not compiled yet.")
  
  observeEvent(input$compile_btn, {
    m <- compile_model_safe(input$ode_custom)
    if (inherits(m, "error")) {
      custom_mod(NULL)
      custom_msg(paste("Compile error:", m$message))
    } else {
      custom_mod(m)
      custom_msg("Custom ODE compiled successfully.")
    }
  })
  
  ode_text <- reactive({
    switch(
      input$model_type,
      base      = make_base_ode(),
      `1trans`  = make_1trans_ode(),
      ntrans    = make_ntrans_ode(input$n_trans),
      base_eff  = make_effect_ode("base"),
      `1trans_eff` = make_effect_ode("1trans"),
      ntrans_eff   = make_effect_ode("ntrans", ntrans = input$n_trans),
      custom    = input$ode_custom
    )
  })
  
  output$ode_show <- renderText(ode_text())
  
  current_model <- reactive({
    if (input$model_type == "custom") {
      custom_mod()
    } else {
      m <- compile_model_safe(ode_text())
      if (inherits(m, "error")) NULL else m
    }
  })
  
  base_params <- reactive({
    p <- c(ka = input$ka, cl = input$cl, v = input$v)
    
    if (input$model_type %in% c("1trans", "ntrans", "1trans_eff", "ntrans_eff")) {
      p <- c(p, ktr = input$ktr)
    }
    
    if (input$model_type %in% c("base_eff", "1trans_eff", "ntrans_eff")) {
      p <- c(p,
             ke0 = input$ke0,
             E0 = input$E0,
             Emax = input$Emax,
             EC50 = input$EC50,
             Hill = input$Hill)
    }
    
    p
  })
  
  output$status <- renderText({
    if (input$model_type == "custom") paste("Custom status:", custom_msg()) else "Ready."
  })
  
  # ---- single scenario ----
  sim_data <- eventReactive(input$run_btn, {
    m <- current_model()
    if (is.null(m)) return(list(err = "Model is not available (compile failed or not compiled).", df = NULL))
    
    ev <- build_event_table(input)
    res <- tryCatch(as.data.frame(rxode2::rxSolve(m, base_params(), ev)), error = function(e) e)
    
    if (inherits(res, "error")) return(list(err = paste("Simulation error:", res$message), df = NULL))
    list(err = NULL, df = res)
  })
  
  output$head_tbl <- renderTable({
    sd <- sim_data()
    if (!is.null(sd$err)) return(data.frame(error = sd$err))
    head(sd$df, 12)
  })
  
  output$pk_plot <- renderPlot({
    sd <- sim_data()
    if (!is.null(sd$err)) { plot.new(); text(0.5, 0.5, sd$err); return() }
    df <- sd$df
    if (!("C1" %in% names(df))) { plot.new(); text(0.5, 0.5, "C1 not found. Define C1=central/v in ODE."); return() }
    
    ggplot(df, aes(x = time, y = C1)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (h)", y = "C1", title = "Single scenario: C1") +
      theme_minimal()
  })
  
  output$pd_plot <- renderPlot({
    sd <- sim_data()
    if (!is.null(sd$err)) { plot.new(); text(0.5, 0.5, sd$err); return() }
    df <- sd$df
    if (!all(c("Ce","E") %in% names(df))) { plot.new(); text(0.5, 0.5, "Ce/E not found. Use an effect-compartment model."); return() }
    
    ggplot(df, aes(x = time, y = E)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (h)", y = "E", title = "Single scenario: Effect (E)") +
      theme_minimal()
  })
  
  output$amt_plot <- renderPlot({
    sd <- sim_data()
    if (!is.null(sd$err)) { plot.new(); text(0.5, 0.5, sd$err); return() }
    df <- sd$df
    
    state_cols <- setdiff(names(df), "time")
    state_cols <- state_cols[sapply(df[state_cols], is.numeric)]
    
    df_long <- df %>%
      select(time, all_of(state_cols)) %>%
      pivot_longer(-time, names_to = "state", values_to = "value")
    
    ggplot(df_long, aes(x = time, y = value)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ state, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = "Value", title = "Single scenario: states") +
      theme_minimal()
  })
  
  output$overlay_summary_tbl <- renderTable({
    ss <- scen_summary()
    if (is.null(ss)) return(data.frame(error = "No scenario data yet. Click 'Run overlay'."))
    ss
  })
  
  output$download_summary_csv <- downloadHandler(
    filename = function() {
      paste0("rxode2_scenario_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      ss <- scen_summary()
      if (is.null(ss)) {
        write.csv(data.frame(error = "No summary available"), file, row.names = FALSE)
      } else {
        write.csv(ss, file, row.names = FALSE)
      }
    }
  )
  
  
  # ---- scenarios overlay ----
  scen_data <- eventReactive(input$run_scen_btn, {
    m <- current_model()
    if (is.null(m)) return(list(err = "Model not available (compile failed or not compiled).", df = NULL))
    
    p0 <- base_params()
    
    # Build scenario table (each row is a scenario with param overrides and/or bolus override)
    scenarios <- NULL
    
    if (input$scenario_mode == "dose") {
      doses <- parse_num_vec(input$dose_list)
      doses <- doses[!is.na(doses)]
      if (length(doses) == 0) return(list(err = "No valid doses parsed. Example: 0,100,250,500 or 0:50:500", df = NULL))
      
      scenarios <- tibble(scenario = paste0("Dose_", doses), bolus_amt = doses)
      
      # If bolus is turned off, dose overlay won’t change anything -> warn as error for clarity
      if (!isTRUE(input$bolus_on)) {
        return(list(err = "Bolus is OFF. Turn on 'Include bolus' to use dose overlay mode.", df = NULL))
      }
    } else {
      p1 <- input$sweep_p1
      v1 <- parse_num_vec(input$sweep_v1)
      if (length(v1) == 0) return(list(err = "No valid values for Parameter 1.", df = NULL))
      
      p2 <- input$sweep_p2
      if (p2 == "None") {
        scenarios <- expand_grid(v1 = v1) %>%
          mutate(scenario = paste0(p1, "=", v1))
      } else {
        v2 <- parse_num_vec(input$sweep_v2)
        if (length(v2) == 0) return(list(err = "Parameter 2 selected but no valid values provided.", df = NULL))
        scenarios <- expand_grid(v1 = v1, v2 = v2) %>%
          mutate(scenario = paste0(p1, "=", v1, "; ", p2, "=", v2))
      }
      
      scenarios <- scenarios %>%
        mutate(p1 = p1, p2 = ifelse(p2 == "None", NA_character_, p2))
    }
    
    
    
    
    # Run all scenarios
    out <- vector("list", nrow(scenarios))
    
    for (i in seq_len(nrow(scenarios))) {
      pi <- p0
      
      # dose override
      ev <- NULL
      if (input$scenario_mode == "dose") {
        ev <- build_event_table(input, bolus_amt_override = scenarios$bolus_amt[i])
      } else {
        # sweep overrides
        pi[scenarios$p1[i]] <- scenarios$v1[i]
        if (!is.na(scenarios$p2[i])) {
          pi[scenarios$p2[i]] <- scenarios$v2[i]
        }
        ev <- build_event_table(input)
      }
      
      res <- tryCatch(as.data.frame(rxode2::rxSolve(m, pi, ev)), error = function(e) e)
      if (inherits(res, "error")) {
        return(list(err = paste("Simulation error in scenario:", scenarios$scenario[i], ":", res$message), df = NULL))
      }
      
      # Attach scenario metadata
      res$scenario <- scenarios$scenario[i]
      # attach the key knobs for export
      if (input$scenario_mode == "dose") {
        res$bolus_amt <- scenarios$bolus_amt[i]
      } else {
        res[[scenarios$p1[i]]] <- scenarios$v1[i]
        if (!is.na(scenarios$p2[i])) res[[scenarios$p2[i]]] <- scenarios$v2[i]
      }
      
      out[[i]] <- res
    }
    
    df_all <- bind_rows(out)
    list(err = NULL, df = df_all)
  })
  
  scen_summary <- reactive({
    sd <- scen_data()
    if (!is.null(sd$err) || is.null(sd$df)) return(NULL)
    summarize_scenarios(sd$df, tau = tau_ss())
  })
  
  
  output$overlay_head <- renderTable({
    sd <- scen_data()
    if (!is.null(sd$err)) return(data.frame(error = sd$err))
    head(sd$df, 12)
  })
  
  output$overlay_pk_plot <- renderPlot({
    sd <- scen_data()
    if (!is.null(sd$err)) { plot.new(); text(0.5, 0.5, sd$err); return() }
    df <- sd$df
    if (!("C1" %in% names(df))) { plot.new(); text(0.5, 0.5, "C1 not found. Define C1=central/v in ODE."); return() }
    
    ggplot(df, aes(x = time, y = C1, group = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "C1", title = "Overlay: C1 across scenarios") +
      theme_minimal()
  })
  
  output$overlay_pd_plot <- renderPlot({
    sd <- scen_data()
    if (!is.null(sd$err)) { plot.new(); text(0.5, 0.5, sd$err); return() }
    df <- sd$df
    
    # Only show if E exists
    if (!("E" %in% names(df))) {
      plot.new(); text(0.5, 0.5, "Effect (E) not available. Choose an effect-compartment model to overlay PD.")
      return()
    }
    
    ggplot(df, aes(x = time, y = E, group = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "E", title = "Overlay: Effect (E) across scenarios") +
      theme_minimal()
  })
  
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("rxode2_scenarios_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      sd <- scen_data()
      if (!is.null(sd$err) || is.null(sd$df)) {
        write.csv(data.frame(error = ifelse(is.null(sd$err), "No data", sd$err)), file, row.names = FALSE)
      } else {
        write.csv(sd$df, file, row.names = FALSE)
      }
    }
  )
}

shinyApp(ui, server)
