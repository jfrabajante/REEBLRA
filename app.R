# app.R
# Bayesian Logistic R0 Shiny App
# Goal: estimate an R0-like threshold equation from ODE simulations using
# Bayesian logistic regression, then compare the learned equation to analytic or derived R0.
#
# Recommended packages:
# install.packages(c(
#   "shiny", "tidyverse", "deSolve", "lhs", "DT", "plotly", "pROC",
#   "caret", "posterior", "rstanarm", "nnet"
# ))
#
# Notes:
# 1) Binary Bayesian logistic regression estimates:
#      Pr(R0 > cutoff or cases increase) = logit^{-1}(a + b1*x1 + ... + bp*xp)
#    This is an empirical threshold equation.
# 2) Successive binary logistic regression estimates ordered brackets:
#      <1, 1-(1+step), ..., >10, with user-selected step size 0.25, 0.5, 1, 2, or 5
# 3) For complex models, the app reports both analytic R0/Re where available and
#    a numerical growth-rate R0estimate. REEBLRA estimates r from early major-host
#    infection growth, reports R-squared and a 95% CI for r, uses r > 0 for
#    positive R0 conversion, adaptively relaxes the susceptible threshold when too few
#    early-window points are available, sets growth-based R0estimate < 1 to 0
#    for non-growing/subthreshold runs, and then learns a surrogate/alternative
#    symbolic closed-form threshold formula for R0 using logistic regression.
# 4) The learned equation is not a replacement for mechanistic derivation; it is
#    a simulation-calibrated approximation useful for screening, sensitivity, and
#    comparison with the next generation matrix result.

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(deSolve)
  library(lhs)
  library(DT)
  library(plotly)
  library(pROC)
  library(caret)
  library(posterior)
})

options(stringsAsFactors = FALSE)
options(mc.cores = max(1, parallel::detectCores() - 1))

safe_require <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

safe_deparse_formula <- function(form) {
  paste(deparse(form), collapse = " ")
}

collapse_long_equation <- function(x) {
  paste(as.character(x), collapse = " ")
}

top_binary_terms_text <- function(bf, top_n = 10) {
  tab <- posterior_coef_table(bf$fit, z_predictors = bf$z_predictors)
  tab <- tab[tab$Term %in% bf$z_predictors, , drop = FALSE]
  if (nrow(tab) == 0) return("No predictor coefficients available.")
  tab$Parameter <- sub("^z_", "", tab$Term)
  tab$AbsMean <- abs(tab$Mean)
  tab <- tab[order(tab$AbsMean, decreasing = TRUE), , drop = FALSE]
  tab <- utils::head(tab, top_n)
  paste0(
    apply(tab, 1, function(row) {
      paste0("\t", row[["Parameter"]], " (standardized beta = ", round(as.numeric(row[["Mean"]]), 5), ")")
    }),
    collapse = "\n"
  )
}

# =========================================================
# MODEL DEFINITIONS
# =========================================================

model_parameter_table <- function(model_name) {
  # Default ranges are intentionally simple for teaching/demo runs.
  # Most ODE rates/proportions are set to min = 0.01 and max = 1.
  # Exception: m is a vector-to-human ratio, not a probability/rate, so values > 1 are appropriate.
  switch(
    model_name,
    "Closed SIR" = tibble::tribble(
      ~parameter, ~min, ~max,
      "beta", 0.01, 1.00,
      "gamma", 0.01, 1.00
    ),
    "Open SIR" = tibble::tribble(
      ~parameter, ~min, ~max,
      "beta", 0.01, 1.00,
      "gamma", 0.01, 1.00,
      "mu", 0.01, 1.00
    ),
    "Open SEIR" = tibble::tribble(
      ~parameter, ~min, ~max,
      "beta", 0.01, 1.00,
      "sigma", 0.01, 1.00,
      "gamma", 0.01, 1.00,
      "mu", 0.01, 1.00
    ),
    "SEIR with vaccination" = tibble::tribble(
      ~parameter, ~min, ~max,
      "beta", 0.01, 1.00,
      "sigma", 0.01, 1.00,
      "gamma", 0.01, 1.00,
      "mu", 0.01, 1.00,
      "v", 0.01, 1.00
    ),
    "Age-structured SIR" = tibble::tribble(
      ~parameter, ~min, ~max,
      "beta", 0.01, 1.00,
      "gamma1", 0.01, 1.00,
      "gamma2", 0.01, 1.00,
      "c11", 0.01, 1.00,
      "c12", 0.01, 1.00,
      "c21", 0.01, 1.00,
      "c22", 0.01, 1.00
    ),
    "Vector-borne" = tibble::tribble(
      ~parameter, ~min, ~max,
      "beta_hv", 0.01, 1.00,
      "beta_vh", 0.01, 1.00,
      "gamma_h", 0.01, 1.00,
      "mu_v", 0.01, 1.00,
      "m", 0.01, 10.00
    )
  )
}

parameter_definition_table <- function(model_name) {
  common <- list(
    beta = "Effective transmission/contact rate.",
    gamma = "Removal rate from the infectious class; this may include recovery, disease-induced death, isolation, natural death if modeled separately in the equation, or loss of infectiousness. 1/gamma is the mean infectious duration when gamma is the only removal process.",
    mu = "Per-capita demographic turnover or natural mortality rate.",
    sigma = "Progression rate from exposed to infectious; 1/sigma is the mean latent period.",
    v = "Vaccination rate or vaccination pressure applied to susceptible individuals.",
    gamma1 = "Recovery rate of infectious individuals in age group 1.",
    gamma2 = "Recovery rate of infectious individuals in age group 2.",
    c11 = "Within-group contact coefficient: contacts from group 1 infectives to group 1 susceptibles.",
    c12 = "Cross-group contact coefficient: contacts from group 2 infectives to group 1 susceptibles.",
    c21 = "Cross-group contact coefficient: contacts from group 1 infectives to group 2 susceptibles.",
    c22 = "Within-group contact coefficient: contacts from group 2 infectives to group 2 susceptibles.",
    beta_hv = "Human-to-vector transmission rate; infectious humans generate newly infected vectors.",
    beta_vh = "Vector-to-human transmission rate; infectious vectors generate newly infected humans.",
    gamma_h = "Human infectious-class removal rate in the vector-borne model.",
    mu_v = "Vector death/removal rate; 1/mu_v is the mean infectious vector lifetime.",
    m = "Vector-to-human ratio. This is not a probability/rate, so values above 1 can be appropriate."
  )

  tibble::tibble(
    parameter = model_parameter_table(model_name)$parameter,
    definition = unname(unlist(common[model_parameter_table(model_name)$parameter]))
  )
}


default_initial_state <- function(model_name, I0 = 0.01) {
  switch(
    model_name,
    "Closed SIR" = c(S = 1 - I0, I = I0, R = 0, C = I0, inc = 0),
    "Open SIR" = c(S = 1 - I0, I = I0, R = 0, C = I0, inc = 0),
    "Open SEIR" = c(S = 1 - I0, E = 0, I = I0, R = 0, C = I0, inc = 0),
    "SEIR with vaccination" = c(S = 1 - I0, V = 0, E = 0, I = I0, R = 0, C = I0, inc = 0),
    "Age-structured SIR" = c(S1 = 0.50 - I0/2, I1 = I0/2, R1 = 0,
                             S2 = 0.50 - I0/2, I2 = I0/2, R2 = 0,
                             C = I0, inc = 0),
    "Vector-borne" = c(Sh = 1 - I0, Ih = I0, Rh = 0,
                       Sv = 1, Iv = 0,
                       Ch = I0, Cv = 0, inc_h = 0, inc_v = 0)
  )
}

ode_closed_sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- max(S + I + R, 1e-12)
    incidence <- beta * S * I / N
    dS <- -incidence
    dI <- incidence - gamma * I
    dR <- gamma * I
    dC <- incidence
    list(c(dS, dI, dR, dC, incidence))
  })
}

ode_open_sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- max(S + I + R, 1e-12)
    incidence <- beta * S * I / N
    births <- mu * N
    dS <- births - incidence - mu * S
    dI <- incidence - gamma * I - mu * I
    dR <- gamma * I - mu * R
    dC <- incidence
    list(c(dS, dI, dR, dC, incidence))
  })
}

ode_open_seir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- max(S + E + I + R, 1e-12)
    incidence <- beta * S * I / N
    births <- mu * N
    dS <- births - incidence - mu * S
    dE <- incidence - sigma * E - mu * E
    dI <- sigma * E - gamma * I - mu * I
    dR <- gamma * I - mu * R
    dC <- incidence
    list(c(dS, dE, dI, dR, dC, incidence))
  })
}

ode_seir_vaccination <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- max(S + V + E + I + R, 1e-12)
    incidence <- beta * S * I / N
    births <- mu * N
    dS <- births - incidence - v * S - mu * S
    dV <- v * S - mu * V
    dE <- incidence - sigma * E - mu * E
    dI <- sigma * E - gamma * I - mu * I
    dR <- gamma * I - mu * R
    dC <- incidence
    list(c(dS, dV, dE, dI, dR, dC, incidence))
  })
}

ode_age_sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- max(S1 + I1 + R1, 1e-12)
    N2 <- max(S2 + I2 + R2, 1e-12)
    lambda1 <- beta * (c11 * I1 / N1 + c12 * I2 / N2)
    lambda2 <- beta * (c21 * I1 / N1 + c22 * I2 / N2)
    inc1 <- lambda1 * S1
    inc2 <- lambda2 * S2
    dS1 <- -inc1
    dI1 <- inc1 - gamma1 * I1
    dR1 <- gamma1 * I1
    dS2 <- -inc2
    dI2 <- inc2 - gamma2 * I2
    dR2 <- gamma2 * I2
    dC <- inc1 + inc2
    list(c(dS1, dI1, dR1, dS2, dI2, dR2, dC, inc1 + inc2))
  })
}

ode_vector_borne <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Human compartments are proportions. Vector compartments are proportions
    # scaled by vector-to-human ratio m in the force on humans.
    inc_h <- beta_vh * m * Sh * Iv
    inc_v <- beta_hv * Sv * Ih
    dSh <- -inc_h
    dIh <- inc_h - gamma_h * Ih
    dRh <- gamma_h * Ih
    dSv <- mu_v * (1 - Sv) - inc_v
    dIv <- inc_v - mu_v * Iv
    dCh <- inc_h
    dCv <- inc_v
    list(c(dSh, dIh, dRh, dSv, dIv, dCh, dCv, inc_h, inc_v))
  })
}

get_ode_function <- function(model_name) {
  switch(
    model_name,
    "Closed SIR" = ode_closed_sir,
    "Open SIR" = ode_open_sir,
    "Open SEIR" = ode_open_seir,
    "SEIR with vaccination" = ode_seir_vaccination,
    "Age-structured SIR" = ode_age_sir,
    "Vector-borne" = ode_vector_borne
  )
}

# =========================================================
# ANALYTIC AND NUMERIC R0
# =========================================================

analytic_r0 <- function(model_name, par) {
  par <- as.list(par)
  switch(
    model_name,
    "Closed SIR" = with(par, beta / gamma),
    "Open SIR" = with(par, beta / (gamma + mu)),
    "Open SEIR" = with(par, beta * sigma / ((sigma + mu) * (gamma + mu))),
    "SEIR with vaccination" = with(par, {
      # Effective/control reproduction number under continuous vaccination rate v.
      # For dS/dt = mu*N - infection - v*S - mu*S, the disease-free susceptible
      # fraction is S*/N = mu/(mu + v), not 1 - v. The basic R0 without
      # vaccination excludes this susceptible-fraction multiplier.
      (mu / (mu + v)) * beta * sigma / ((sigma + mu) * (gamma + mu))
    }),
    "Age-structured SIR" = with(par, {
      K <- matrix(c(
        beta * c11 / gamma1, beta * c12 / gamma2,
        beta * c21 / gamma1, beta * c22 / gamma2
      ), nrow = 2, byrow = TRUE)
      max(Re(eigen(K)$values))
    }),
    "Vector-borne" = with(par, {
      sqrt((beta_hv * beta_vh * m) / (gamma_h * mu_v))
    })
  )
}

r0_formula_text <- function(model_name) {
  switch(
    model_name,
    "Closed SIR" = "Analytic R0 = beta / gamma",
    "Open SIR" = "Analytic R0 = beta / (gamma + mu)",
    "Open SEIR" = "Analytic R0 = beta*sigma / ((sigma + mu)*(gamma + mu))",
    "SEIR with vaccination" = "Effective R under continuous vaccination = (mu/(mu + v))*beta*sigma / ((sigma + mu)*(gamma + mu)); basic R0 without vaccination excludes the susceptible-fraction multiplier.",
    "Age-structured SIR" = "Analytic R0 = spectral radius of K, where K_ij = beta*c_ij/gamma_j, assuming the two groups are equally sized and fully susceptible at the disease-free state.",
    "Vector-borne" = "Analytic R0 = sqrt(beta_hv*beta_vh*m/(gamma_h*mu_v)); geometric mean of human-to-vector and vector-to-human transmission."
  )
}


ode_equation_text <- function(model_name) {
  switch(
    model_name,
    "Closed SIR" = paste(
      "Closed SIR model:",
      "dS/dt = - beta*S*I/N",
      "dI/dt = beta*S*I/N - gamma*I",
      "dR/dt = gamma*I",
      "dC/dt = beta*S*I/N",
      "N = S + I + R",
      "Analytic R0 = beta/gamma",
      sep = "\n"
    ),
    "Open SIR" = paste(
      "Open SIR model with demographic turnover:",
      "dS/dt = mu*N - beta*S*I/N - mu*S",
      "dI/dt = beta*S*I/N - gamma*I - mu*I",
      "dR/dt = gamma*I - mu*R",
      "dC/dt = beta*S*I/N",
      "N = S + I + R",
      "Analytic R0 = beta/(gamma + mu)",
      sep = "\n"
    ),
    "Open SEIR" = paste(
      "Open SEIR model with exposed class and demographic turnover:",
      "dS/dt = mu*N - beta*S*I/N - mu*S",
      "dE/dt = beta*S*I/N - sigma*E - mu*E",
      "dI/dt = sigma*E - gamma*I - mu*I",
      "dR/dt = gamma*I - mu*R",
      "dC/dt = beta*S*I/N",
      "N = S + E + I + R",
      "Analytic R0 = beta*sigma/((sigma + mu)*(gamma + mu))",
      sep = "\n"
    ),
    "SEIR with vaccination" = paste(
      "SEIR model with vaccination and demographic turnover:",
      "dS/dt = mu*N - beta*S*I/N - v*S - mu*S",
      "dV/dt = v*S - mu*V",
      "dE/dt = beta*S*I/N - sigma*E - mu*E",
      "dI/dt = sigma*E - gamma*I - mu*I",
      "dR/dt = gamma*I - mu*R",
      "dC/dt = beta*S*I/N",
      "N = S + V + E + I + R",
      "Effective R under continuous vaccination = (mu/(mu + v))*beta*sigma/((sigma + mu)*(gamma + mu))",
      sep = "\n"
    ),
    "Age-structured SIR" = paste(
      "Two-age-group SIR model:",
      "lambda1 = beta*(c11*I1/N1 + c12*I2/N2)",
      "lambda2 = beta*(c21*I1/N1 + c22*I2/N2)",
      "dS1/dt = -lambda1*S1",
      "dI1/dt = lambda1*S1 - gamma1*I1",
      "dR1/dt = gamma1*I1",
      "dS2/dt = -lambda2*S2",
      "dI2/dt = lambda2*S2 - gamma2*I2",
      "dR2/dt = gamma2*I2",
      "dC/dt = lambda1*S1 + lambda2*S2",
      "Analytic R0 = spectral radius of K, where K_ij = beta*c_ij/gamma_j, assuming equally sized fully susceptible groups",
      sep = "\n"
    ),
    "Vector-borne" = paste(
      "Vector-borne model:",
      "dSh/dt = - beta_vh*m*Sh*Iv",
      "dIh/dt = beta_vh*m*Sh*Iv - gamma_h*Ih",
      "dRh/dt = gamma_h*Ih",
      "dSv/dt = mu_v*(1 - Sv) - beta_hv*Sv*Ih",
      "dIv/dt = beta_hv*Sv*Ih - mu_v*Iv",
      "dCh/dt = beta_vh*m*Sh*Iv",
      "dCv/dt = beta_hv*Sv*Ih",
      "Analytic R0 = sqrt(beta_hv*beta_vh*m/(gamma_h*mu_v))",
      sep = "\n"
    )
  )
}

infectious_column <- function(model_name) {
  switch(
    model_name,
    "Closed SIR" = "I",
    "Open SIR" = "I",
    "Open SEIR" = "I",
    "SEIR with vaccination" = "I",
    "Age-structured SIR" = "Itotal",
    "Vector-borne" = "Ih"
  )
}


susceptible_column <- function(model_name) {
  switch(
    model_name,
    "Closed SIR" = "S",
    "Open SIR" = "S",
    "Open SEIR" = "S",
    "SEIR with vaccination" = "S",
    "Age-structured SIR" = "Stotal",
    "Vector-borne" = "Sh"
  )
}

prepare_sim_output <- function(model_name, sim) {
  sim <- as.data.frame(sim)
  if (model_name == "Age-structured SIR") {
    sim$Itotal <- sim$I1 + sim$I2
    sim$Stotal <- sim$S1 + sim$S2
  }
  sim
}

average_removal_rate <- function(model_name, par) {
  par <- as.list(par)
  switch(
    model_name,
    "Closed SIR" = par$gamma,
    "Open SIR" = par$gamma + par$mu,
    "Open SEIR" = par$gamma + par$mu,
    "SEIR with vaccination" = par$gamma + par$mu,
    "Age-structured SIR" = mean(c(par$gamma1, par$gamma2), na.rm = TRUE),
    "Vector-borne" = par$gamma_h
  )
}

average_latency_rate <- function(model_name, par) {
  par <- as.list(par)
  switch(
    model_name,
    "Closed SIR" = NA_real_,
    "Open SIR" = NA_real_,
    "Open SEIR" = par$sigma + par$mu,
    "SEIR with vaccination" = par$sigma + par$mu,
    "Age-structured SIR" = NA_real_,
    "Vector-borne" = NA_real_
  )
}

growth_rate_r0 <- function(model_name, sim, par, I0, susceptible_threshold = 0.95, min_points = 5) {
  # Growth-rate numerical R0 estimate for REEBLRA.
  #
  # REEBLRA first estimates the early exponential growth rate r from the
  # major-host infectious trajectory:
  #   log(I_major(t)) = a + r*t.
  #
  # The fitting window is restricted to times when the susceptible major-host
  # population is still approximately undepleted:
  #   S_major(t) >= susceptible_threshold * (S_major(0)).
  #
  # Then REEBLRA converts r into a numeric R0estimate using a simple
  # model-family approximation:
  #   SIR-like:  R0estimate = 1 + r / gamma_removal
  #   SEIR-like: R0estimate = ((r + sigma_latency)*(r + gamma_removal)) /
  #                            (sigma_latency*gamma_removal)
  #
  # Here gamma_removal is removal from the infectious state, including recovery,
  # disease-induced death, isolation, natural death if modeled, or loss of
  # infectiousness. For vector-borne systems, the major host is fitted and vector
  # transmission is treated as hidden in the observed early growth of I_major.
  #
  # The function reports r, 95% CI for r, R^2 of the log-linear fit, the number
  # of points used, threshold actually used, and warning flags. If too few points
  # are available at the chosen threshold, REEBLRA progressively relaxes the
  # susceptible threshold and, if needed, uses the earliest positive infectious
  # points within the simulation window. If the resulting growth-based R0estimate
  # is below 1, it is set to 0 to flag a non-growing/subthreshold numerical run.
  par <- as.list(par)
  sim <- as.data.frame(sim)
  Icol <- infectious_column(model_name)
  Scol <- susceptible_column(model_name)

  if (!Icol %in% names(sim) || !Scol %in% names(sim)) {
    return(tibble::tibble(
      numeric_secondary_R0 = NA_real_, early_growth_r = NA_real_,
      early_growth_r_lwr95 = NA_real_, early_growth_r_upr95 = NA_real_,
      early_growth_R2 = NA_real_, early_growth_n = 0L,
      early_growth_threshold_used = NA_real_,
      early_growth_window_end = NA_real_, early_growth_valid = FALSE,
      early_growth_warning = "Missing infectious or susceptible column."
    ))
  }

  I <- as.numeric(sim[[Icol]])
  S <- as.numeric(sim[[Scol]])
  time <- as.numeric(sim$time)
  gamma_removal <- average_removal_rate(model_name, par)
  sigma_latency <- average_latency_rate(model_name, par)

  S0_adj <- max(S[1], 1e-12)

  # Adaptive fitting-window rule:
  # start with the user-selected susceptible threshold, then relax stepwise if
  # fewer than min_points are available. This preserves early-growth fitting
  # whenever possible while avoiding unnecessary failed runs when the time step
  # is coarse or the epidemic depletes susceptibles quickly.
  candidate_thresholds <- sort(unique(c(
    susceptible_threshold,
    0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50
  )), decreasing = TRUE)
  candidate_thresholds <- candidate_thresholds[candidate_thresholds <= susceptible_threshold + 1e-12]
  if (length(candidate_thresholds) == 0) candidate_thresholds <- susceptible_threshold

  d <- data.frame()
  threshold_used <- NA_real_
  used_fallback_points <- FALSE

  for (thr in candidate_thresholds) {
    keep <- is.finite(time) & is.finite(I) & is.finite(S) &
      I > 0 & S >= thr * S0_adj
    d_try <- data.frame(time = time[keep], I = I[keep])
    d_try <- d_try[order(d_try$time), , drop = FALSE]
    if (nrow(d_try) >= min_points) {
      d <- d_try
      threshold_used <- thr
      break
    }
    if (nrow(d_try) > nrow(d)) {
      d <- d_try
      threshold_used <- thr
    }
  }

  if (nrow(d) < min_points) {
    keep_any <- is.finite(time) & is.finite(I) & I > 0
    d_any <- data.frame(time = time[keep_any], I = I[keep_any])
    d_any <- d_any[order(d_any$time), , drop = FALSE]
    if (nrow(d_any) >= min_points) {
      d <- utils::head(d_any, min_points)
      threshold_used <- NA_real_
      used_fallback_points <- TRUE
    }
  }

  if (nrow(d) < min_points) {
    return(tibble::tibble(
      numeric_secondary_R0 = NA_real_, early_growth_r = NA_real_,
      early_growth_r_lwr95 = NA_real_, early_growth_r_upr95 = NA_real_,
      early_growth_R2 = NA_real_, early_growth_n = nrow(d),
      early_growth_threshold_used = threshold_used,
      early_growth_window_end = ifelse(nrow(d) > 0, max(d$time), NA_real_),
      early_growth_valid = FALSE,
      early_growth_warning = paste0("Too few positive infectious points for r fitting within the simulation window; need at least ", min_points, ".")
    ))
  }

  fit <- tryCatch(stats::lm(log(I) ~ time, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(
      numeric_secondary_R0 = NA_real_, early_growth_r = NA_real_,
      early_growth_r_lwr95 = NA_real_, early_growth_r_upr95 = NA_real_,
      early_growth_R2 = NA_real_, early_growth_n = nrow(d),
      early_growth_threshold_used = threshold_used,
      early_growth_window_end = max(d$time),
      early_growth_valid = FALSE,
      early_growth_warning = "Log-linear r fit failed."
    ))
  }

  sm <- summary(fit)
  r_hat <- unname(stats::coef(fit)[["time"]])
  r2 <- unname(sm$r.squared)
  ci <- tryCatch(stats::confint(fit, "time", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
  r_lwr <- as.numeric(ci[1, 1])
  r_upr <- as.numeric(ci[1, 2])

  warning_msg <- character(0)
  if (isTRUE(used_fallback_points)) {
    warning_msg <- c(warning_msg, "Used earliest positive infectious points because no susceptible-threshold window had enough points.")
  } else if (is.finite(threshold_used) && threshold_used < susceptible_threshold) {
    warning_msg <- c(warning_msg, paste0("Susceptible threshold was relaxed from ", susceptible_threshold, " to ", threshold_used, " to obtain enough fitting points."))
  }
  valid <- TRUE
  if (!is.finite(r_hat) || r_hat <= 0) {
    valid <- FALSE
    warning_msg <- c(warning_msg, "Estimated r is not positive; growth-based R0estimate is set to 0 by convention.")
  }
  if (!is.finite(gamma_removal) || gamma_removal <= 0) {
    valid <- FALSE
    warning_msg <- c(warning_msg, "Removal rate is missing or non-positive.")
  }
  if (is.finite(r_lwr) && is.finite(r_upr) && r_lwr <= 0 && r_upr >= 0) {
    warning_msg <- c(warning_msg, "The 95% CI for r crosses 0.")
  }
  if (is.finite(r2) && r2 < 0.80) {
    warning_msg <- c(warning_msg, "Low R-squared for early log-linear fit; inspect fitting window.")
  }

  r0_est <- NA_real_
  if (is.finite(r_hat) && r_hat <= 0) {
    r0_est <- 0
  } else if (isTRUE(valid)) {
    if (model_name %in% c("Open SEIR", "SEIR with vaccination") &&
        is.finite(sigma_latency) && sigma_latency > 0) {
      r0_est <- ((r_hat + sigma_latency) * (r_hat + gamma_removal)) /
        (sigma_latency * gamma_removal)
    } else {
      r0_est <- 1 + r_hat / gamma_removal
    }
    if (is.finite(r0_est) && r0_est < 1) {
      r0_est <- 0
      warning_msg <- c(warning_msg, "Growth-based R0estimate was below 1 and was set to 0 by convention.")
    }
  }

  tibble::tibble(
    numeric_secondary_R0 = r0_est,
    early_growth_r = r_hat,
    early_growth_r_lwr95 = r_lwr,
    early_growth_r_upr95 = r_upr,
    early_growth_R2 = r2,
    early_growth_n = nrow(d),
    early_growth_threshold_used = threshold_used,
    early_growth_window_end = max(d$time),
    early_growth_valid = isTRUE(valid),
    early_growth_warning = ifelse(length(warning_msg) == 0, "OK", paste(unique(warning_msg), collapse = " "))
  )
}

format_bracket_number <- function(x) {
  # Clean labels: 1 instead of 1.00; 1.25 preserved; keep 10 as 10.
  x_chr <- format(round(as.numeric(x), 4), trim = TRUE, scientific = FALSE)
  x_chr <- sub("\\.0+$", "", x_chr)
  x_chr <- sub("(\\.[0-9]*?)0+$", "\\1", x_chr)
  x_chr
}

make_r0_brackets <- function(x, step = 0.25) {
  # User-selectable ordered R0 brackets:
  #   step = 0.25: <1, 1-1.25, 1.25-1.5, ..., >10
  #   step = 0.50: <1, 1-1.5, 1.5-2, ..., >10
  #   step = 1.00: <1, 1-2, 2-3, ..., >10
  step <- as.numeric(step)
  if (is.na(step) || !step %in% c(0.25, 0.5, 1, 2, 5)) step <- 0.25

  thresholds <- c(1, seq(1 + step, 10, by = step))
  thresholds <- thresholds[thresholds <= 10 + 1e-9]
  if (tail(thresholds, 1) < 10) thresholds <- c(thresholds, 10)

  breaks <- c(-Inf, thresholds, Inf)
  interval_labels <- paste0(
    vapply(head(thresholds, -1), format_bracket_number, character(1)),
    "-",
    vapply(tail(thresholds, -1), format_bracket_number, character(1))
  )
  labels <- c("<1", interval_labels, ">10")

  cut(x, breaks = breaks, labels = labels, right = TRUE, include.lowest = TRUE)
}

latin_hypercube_samples <- function(param_table, n) {
  samples <- lhs::randomLHS(n, nrow(param_table))
  colnames(samples) <- param_table$parameter
  samples <- as.data.frame(samples)
  for (j in seq_len(nrow(param_table))) {
    pname <- param_table$parameter[j]
    samples[[pname]] <- param_table$min[j] + samples[[pname]] * (param_table$max[j] - param_table$min[j])
  }
  samples
}

simulate_one <- function(model_name, par_row, times, I0, growth_epsilon, susceptible_threshold = 0.95) {
  state <- default_initial_state(model_name, I0 = I0)
  ode_fun <- get_ode_function(model_name)
  sim <- ode(
    y = state,
    times = times,
    func = ode_fun,
    parms = as.list(par_row),
    method = "lsoda"
  )
  sim <- prepare_sim_output(model_name, sim)
  Icol <- infectious_column(model_name)
  max_I <- max(sim[[Icol]], na.rm = TRUE)
  initial_I <- sim[[Icol]][1]
  cases_increase <- as.integer(max_I > initial_I + growth_epsilon)
  true_R0 <- analytic_r0(model_name, par_row)
  growth_r0 <- growth_rate_r0(
    model_name = model_name,
    sim = sim,
    par = par_row,
    I0 = I0,
    susceptible_threshold = susceptible_threshold
  )

  tibble::as_tibble(par_row) %>%
    mutate(
      analytic_R0 = true_R0,
      cases_increase = cases_increase,
      max_I = max_I,
      final_size = dplyr::last(sim[[grep("^C$|^Ch$", names(sim), value = TRUE)[1]]])
    ) %>%
    dplyr::bind_cols(growth_r0)
}

simulate_many <- function(model_name, param_table, n, tmax, dt, I0, growth_epsilon, susceptible_threshold = 0.95, show_progress = FALSE) {
  samples <- latin_hypercube_samples(param_table, n)
  times <- seq(0, tmax, by = dt)
  out <- vector("list", n)

  for (i in seq_len(n)) {
    if (isTRUE(show_progress)) {
      shiny::incProgress(1 / n, detail = paste("Simulation", i, "of", n))
    }
    out[[i]] <- tryCatch(
      simulate_one(model_name, samples[i, , drop = FALSE], times, I0, growth_epsilon, susceptible_threshold = susceptible_threshold),
      error = function(e) {
        tibble::as_tibble(samples[i, , drop = FALSE]) %>%
          mutate(
            analytic_R0 = NA_real_,
            numeric_secondary_R0 = NA_real_,
            early_growth_r = NA_real_,
            early_growth_r_lwr95 = NA_real_,
            early_growth_r_upr95 = NA_real_,
            early_growth_R2 = NA_real_,
            early_growth_n = NA_integer_,
            early_growth_threshold_used = NA_real_,
            early_growth_window_end = NA_real_,
            early_growth_valid = FALSE,
            early_growth_warning = "Simulation or growth-rate fitting failed.",
            cases_increase = NA_integer_,
            max_I = NA_real_,
            final_size = NA_real_
          )
      }
    )
  }

  # Keep all attempted simulation rows. Do not use drop_na() here because
  # invalid/uncertain growth-rate fits intentionally produce NA R0estimate or
  # NA confidence limits. Dropping all NA rows can make the displayed number of
  # simulations much smaller than the user-requested number. Downstream fitting
  # functions already use complete cases for the selected outcome/predictors.
  bind_rows(out) %>%
    mutate(
      R0_bracket_analytic = make_r0_brackets(analytic_R0),
      R0_bracket_numeric = make_r0_brackets(numeric_secondary_R0)
    )
}

# =========================================================
# MODEL FITTING HELPERS
# =========================================================

standardize_predictors <- function(df, predictors) {
  centers <- sapply(df[predictors], mean, na.rm = TRUE)
  scales <- sapply(df[predictors], sd, na.rm = TRUE)
  scales[is.na(scales) | scales == 0] <- 1
  out <- df
  for (p in predictors) {
    out[[paste0("z_", p)]] <- (out[[p]] - centers[[p]]) / scales[[p]]
  }
  list(data = out, centers = centers, scales = scales, z_predictors = paste0("z_", predictors))
}

binary_metrics <- function(y_true, prob, cutoff = 0.5, validation_label = "Apparent / in-sample") {
  y_true <- as.integer(y_true)
  pred <- ifelse(prob >= cutoff, 1, 0)
  out <- tibble::tibble(
    Validation = validation_label,
    Metric = c("Accuracy", "AUC", "Sensitivity", "Specificity", "Precision", "MSE"),
    Value = NA_real_
  )
  out$Value[out$Metric == "Accuracy"] <- mean(pred == y_true, na.rm = TRUE)
  out$Value[out$Metric == "MSE"] <- mean((prob - y_true)^2, na.rm = TRUE)
  out$Value[out$Metric == "AUC"] <- tryCatch(as.numeric(pROC::auc(pROC::roc(y_true, prob, quiet = TRUE))), error = function(e) NA_real_)

  cm <- tryCatch(caret::confusionMatrix(factor(pred, levels = c(0, 1)),
                                        factor(y_true, levels = c(0, 1)),
                                        positive = "1"), error = function(e) NULL)
  if (!is.null(cm)) {
    out$Value[out$Metric == "Sensitivity"] <- unname(cm$byClass["Sensitivity"])
    out$Value[out$Metric == "Specificity"] <- unname(cm$byClass["Specificity"])
    out$Value[out$Metric == "Precision"] <- unname(cm$byClass["Pos Pred Value"])
  }
  out
}

make_stratified_folds <- function(y, k = 5, seed = 123) {
  y <- as.factor(y)
  n <- length(y)
  k <- as.integer(k)
  if (is.na(k) || k < 2) return(rep(1L, n))
  k <- min(k, n)
  set.seed(seed)
  folds <- integer(n)
  for (cls in levels(y)) {
    idx <- which(y == cls)
    if (length(idx) == 0) next
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

fit_stan_glm_binary <- function(form, data, chains, iter, seed, algorithm = "meanfield") {
  if (algorithm == "sampling") {
    rstanarm::stan_glm(
      form,
      data = data,
      family = stats::binomial(link = "logit"),
      prior = rstanarm::normal(location = 0, scale = 2, autoscale = TRUE),
      prior_intercept = rstanarm::normal(location = 0, scale = 2, autoscale = TRUE),
      chains = chains,
      iter = iter,
      seed = seed,
      refresh = 0,
      algorithm = "sampling"
    )
  } else {
    rstanarm::stan_glm(
      form,
      data = data,
      family = stats::binomial(link = "logit"),
      prior = rstanarm::normal(location = 0, scale = 2, autoscale = TRUE),
      prior_intercept = rstanarm::normal(location = 0, scale = 2, autoscale = TRUE),
      iter = iter,
      seed = seed,
      refresh = 0,
      algorithm = algorithm
    )
  }
}


posterior_epred_mean <- function(fit, newdata) {
  # rstanarm::posterior_epred usually returns a posterior-draws x observations matrix.
  # For some algorithms/versions it may return a vector. This wrapper handles both.
  pp <- rstanarm::posterior_epred(fit, newdata = newdata)
  if (is.null(dim(pp))) {
    return(as.numeric(pp))
  }
  if (length(dim(pp)) == 2) {
    return(as.numeric(colMeans(pp)))
  }
  # Defensive fallback for unexpected arrays.
  as.numeric(apply(pp, 2, mean))
}

extract_bayes_coef_means <- function(fit, z_predictors = NULL) {
  # Robust coefficient extraction for rstanarm::stan_glm.
  # Avoids summary(fit)$coefficients because some rstanarm/algorithm combinations
  # return summary objects that are not list-like, causing "$ operator is invalid
  # for atomic vectors".
  est <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if (!is.null(est)) {
    est <- unlist(est)
    est <- stats::setNames(as.numeric(est), names(est))
    if (length(est) > 0 && !any(is.na(names(est))) && all(nzchar(names(est)))) {
      return(est)
    }
  }

  draws <- tryCatch(as.matrix(fit), error = function(e) NULL)
  if (!is.null(draws)) {
    wanted <- unique(c("(Intercept)", z_predictors))
    terms <- intersect(wanted, colnames(draws))
    if (length(terms) == 0) {
      # Last fallback: use all non-auxiliary columns that look like coefficients.
      terms <- grep("^\\(Intercept\\)$|^z_", colnames(draws), value = TRUE)
    }
    if (length(terms) > 0) {
      return(colMeans(draws[, terms, drop = FALSE]))
    }
  }

  stop("Could not extract posterior coefficient means from the Bayesian logistic model.")
}

posterior_coef_table <- function(fit, z_predictors = NULL) {
  draws <- tryCatch(as.matrix(fit), error = function(e) NULL)
  wanted <- unique(c("(Intercept)", z_predictors))

  if (!is.null(draws)) {
    terms <- intersect(wanted, colnames(draws))
    if (length(terms) == 0) {
      terms <- grep("^\\(Intercept\\)$|^z_", colnames(draws), value = TRUE)
    }
    if (length(terms) > 0) {
      out <- data.frame(
        Term = terms,
        Mean = colMeans(draws[, terms, drop = FALSE]),
        SD = apply(draws[, terms, drop = FALSE], 2, stats::sd),
        Q2.5 = apply(draws[, terms, drop = FALSE], 2, stats::quantile, probs = 0.025),
        Median = apply(draws[, terms, drop = FALSE], 2, stats::median),
        Q97.5 = apply(draws[, terms, drop = FALSE], 2, stats::quantile, probs = 0.975),
        row.names = NULL,
        check.names = FALSE
      )
      return(out)
    }
  }

  est <- extract_bayes_coef_means(fit, z_predictors)
  data.frame(
    Term = names(est),
    Mean = as.numeric(est),
    SD = NA_real_,
    Q2.5 = NA_real_,
    Median = NA_real_,
    Q97.5 = NA_real_,
    row.names = NULL,
    check.names = FALSE
  )
}

posterior_equation <- function(fit, centers = NULL, scales = NULL, z_predictors = NULL) {
  # Extract posterior mean coefficients robustly from a stan_glm object.
  est <- extract_bayes_coef_means(fit, z_predictors = z_predictors)

  intercept_name <- if ("(Intercept)" %in% names(est)) "(Intercept)" else names(est)[1]
  intercept_z <- unname(est[intercept_name])
  betas_z <- est[setdiff(names(est), intercept_name)]

  # Keep only the standardized predictors in the equation. This prevents auxiliary
  # Stan columns from entering the printed formula.
  if (!is.null(z_predictors)) {
    betas_z <- betas_z[names(betas_z) %in% z_predictors]
  }

  build_rhs <- function(intercept, betas) {
    rhs <- ""
    if (length(betas) > 0) {
      rhs <- paste(
        paste0(
          ifelse(betas >= 0, " + ", " - "),
          abs(round(as.numeric(betas), 5)),
          "*", names(betas)
        ),
        collapse = ""
      )
    }
    paste0(round(as.numeric(intercept), 5), rhs)
  }

  rhs_z <- build_rhs(intercept_z, betas_z)

  # Convert from standardized predictors back to original ODE parameters:
  # z_x = (x - mean_x)/sd_x. Therefore,
  # eta = a_z + sum b_z*z_x = a_original + sum b_original*x.
  rhs_original <- NA_character_
  original_available <- !is.null(centers) && !is.null(scales) && length(betas_z) > 0

  if (original_available) {
    param_names <- sub("^z_", "", names(betas_z))
    valid <- param_names %in% names(centers) & param_names %in% names(scales)
    if (all(valid)) {
      b_original <- as.numeric(betas_z) / as.numeric(scales[param_names])
      names(b_original) <- param_names
      a_original <- intercept_z - sum(as.numeric(betas_z) * as.numeric(centers[param_names]) / as.numeric(scales[param_names]))
      rhs_original <- build_rhs(a_original, b_original)
    }
  }

  list(
    logit_standardized = paste0("logit(p) = ", rhs_z),
    probability_standardized = paste0("p = 1 / (1 + exp(-(", rhs_z, ")))"),
    threshold_standardized = paste0("Estimated threshold score: eta = ", rhs_z, "; classify as threshold TRUE when eta > 0 or p > 0.5."),
    logit_original = ifelse(is.na(rhs_original), "Original-scale equation unavailable; check coefficient extraction.", paste0("logit(p) = ", rhs_original)),
    probability_original = ifelse(is.na(rhs_original), "Original-scale probability equation unavailable; check coefficient extraction.", paste0("p = 1 / (1 + exp(-(", rhs_original, ")))")),
    threshold_original = ifelse(is.na(rhs_original), "Original-scale threshold equation unavailable; check coefficient extraction.", paste0("Estimated R0-threshold score on original parameter scale: eta = ", rhs_original, "; classify as threshold TRUE when eta > 0 or p > 0.5."))
  )
}


ordinary_logistic_equation <- function(fit, centers = NULL, scales = NULL, z_predictors = NULL) {
  co <- stats::coef(fit)
  intercept_name <- if ("(Intercept)" %in% names(co)) "(Intercept)" else names(co)[1]
  intercept_z <- unname(co[intercept_name])
  betas_z <- co[setdiff(names(co), intercept_name)]

  if (!is.null(z_predictors)) {
    betas_z <- betas_z[names(betas_z) %in% z_predictors]
  }

  build_rhs <- function(intercept, betas) {
    rhs <- ""
    if (length(betas) > 0) {
      rhs <- paste(
        paste0(
          ifelse(betas >= 0, " + ", " - "),
          abs(round(as.numeric(betas), 5)),
          "*", names(betas)
        ),
        collapse = ""
      )
    }
    paste0(round(as.numeric(intercept), 5), rhs)
  }

  rhs_z <- build_rhs(intercept_z, betas_z)
  rhs_original <- NA_character_

  if (!is.null(centers) && !is.null(scales) && length(betas_z) > 0) {
    param_names <- sub("^z_", "", names(betas_z))
    valid <- param_names %in% names(centers) & param_names %in% names(scales)
    if (all(valid)) {
      b_original <- as.numeric(betas_z) / as.numeric(scales[param_names])
      names(b_original) <- param_names
      a_original <- intercept_z -
        sum(as.numeric(betas_z) * as.numeric(centers[param_names]) / as.numeric(scales[param_names]))
      rhs_original <- build_rhs(a_original, b_original)
    }
  }

  list(
    eta_standardized = paste0("eta = ", rhs_z),
    probability_standardized = paste0("p = 1 / (1 + exp(-(", rhs_z, ")))") ,
    eta_original = ifelse(
      is.na(rhs_original),
      "Original-scale eta equation unavailable; check coefficient extraction.",
      paste0("eta = ", rhs_original)
    ),
    probability_original = ifelse(
      is.na(rhs_original),
      "Original-scale probability equation unavailable; check coefficient extraction.",
      paste0("p = 1 / (1 + exp(-(", rhs_original, ")))")
    )
  )
}

fit_binary_bayes <- function(df, predictors, outcome, chains, iter, seed, algorithm = "meanfield", k_folds = 2) {
  if (!safe_require("rstanarm")) {
    stop("Package 'rstanarm' is required for binary Bayesian logistic regression. Install it first with install.packages('rstanarm').")
  }

  keep_model <- c(predictors, outcome)
  d_full <- df[stats::complete.cases(df[, keep_model, drop = FALSE]), , drop = FALSE]
  d_full[[outcome]] <- as.integer(d_full[[outcome]])

  if (nrow(d_full) < 20) {
    stop("Too few complete rows for Bayesian logistic regression. Increase sample size or check parameter ranges.")
  }

  class_counts <- table(factor(d_full[[outcome]], levels = c(0, 1)))
  if (length(unique(d_full[[outcome]])) < 2) {
    stop(paste0(
      "Binary outcome has only one class after simulation (0 = ", class_counts[["0"]],
      ", 1 = ", class_counts[["1"]], "). A binary logistic model cannot be fitted. ",
      "Use the Simulation results tab to check Threshold positives/negatives, then widen parameter ranges, ",
      "increase n_samples, change the R0 cutoff, or switch the threshold source."
    ))
  }

  k_folds <- as.integer(k_folds)
  if (is.na(k_folds) || k_folds < 1) k_folds <- 1
  k_folds <- min(k_folds, nrow(d_full))

  # Final model for the reported equation: fit on all available complete data.
  scaled <- standardize_predictors(d_full, predictors)
  d2 <- scaled$data
  zpred <- scaled$z_predictors
  form <- as.formula(paste(outcome, "~", paste(zpred, collapse = " + ")))

  fit <- fit_stan_glm_binary(
    form = form,
    data = d2,
    chains = chains,
    iter = iter,
    seed = seed,
    algorithm = algorithm
  )

  ordinary_fit <- tryCatch(
    stats::glm(form, data = d2, family = stats::binomial(link = "logit")),
    error = function(e) NULL
  )
  ordinary_equation <- if (is.null(ordinary_fit)) {
    NULL
  } else {
    ordinary_logistic_equation(
      ordinary_fit,
      centers = scaled$centers,
      scales = scaled$scales,
      z_predictors = zpred
    )
  }

  prob <- posterior_epred_mean(fit, d2)
  metrics_apparent <- binary_metrics(d2[[outcome]], prob, validation_label = "Apparent / in-sample")
  cv_predictions <- NULL
  metrics <- metrics_apparent

  # K-fold validation: fit each fold on training data and predict held-out rows.
  # This gives a less optimistic validation metric than the apparent/in-sample values.
  if (k_folds >= 2) {
    folds <- make_stratified_folds(d_full[[outcome]], k = k_folds, seed = seed)
    cv_rows <- vector("list", k_folds)

    for (fold in seq_len(k_folds)) {
      train_idx <- which(folds != fold)
      test_idx <- which(folds == fold)
      if (length(test_idx) == 0 || length(unique(d_full[[outcome]][train_idx])) < 2) next

      train_raw <- d_full[train_idx, , drop = FALSE]
      test_raw <- d_full[test_idx, , drop = FALSE]
      train_scaled <- standardize_predictors(train_raw, predictors)
      train_d <- train_scaled$data
      test_d <- standardize_newdata(test_raw, train_scaled$centers, train_scaled$scales, predictors)
      fold_form <- as.formula(paste(outcome, "~", paste(train_scaled$z_predictors, collapse = " + ")))

      fold_fit <- fit_stan_glm_binary(
        form = fold_form,
        data = train_d,
        chains = chains,
        iter = iter,
        seed = seed + fold,
        algorithm = algorithm
      )

      fold_prob <- posterior_epred_mean(fold_fit, test_d)
      cv_rows[[fold]] <- test_raw %>%
        dplyr::mutate(
          .row_id = test_idx,
          .fold = fold,
          .observed = .data[[outcome]],
          .cv_predicted_prob = as.numeric(fold_prob),
          .cv_predicted_class = as.integer(.cv_predicted_prob >= 0.5)
        )
    }

    cv_predictions <- dplyr::bind_rows(cv_rows)
    if (!is.null(cv_predictions) && nrow(cv_predictions) > 0) {
      metrics_cv <- binary_metrics(
        cv_predictions$.observed,
        cv_predictions$.cv_predicted_prob,
        validation_label = paste0(k_folds, "-fold cross-validation")
      )
      metrics <- dplyr::bind_rows(metrics_cv, metrics_apparent)
    }
  }

  eq <- posterior_equation(fit, centers = scaled$centers, scales = scaled$scales, z_predictors = zpred)

  list(
    fit = fit,
    ordinary_fit = ordinary_fit,
    ordinary_equation = ordinary_equation,
    data = d2,
    predictors = predictors,
    z_predictors = zpred,
    centers = scaled$centers,
    scales = scaled$scales,
    prob = prob,
    metrics = metrics,
    cv_predictions = cv_predictions,
    equation = eq,
    outcome = outcome,
    formula = form,
    algorithm = algorithm,
    class_counts = class_counts,
    k_folds = k_folds
  )
}



fit_multiclass_model <- function(df, predictors, outcome_class, use_bayes = FALSE, chains = NULL, iter = NULL, seed = NULL, bracket_step = 0.25, k_folds = 1) {
  # Successive binary logistic R0 bracket estimation.
  #
  # Instead of fitting one multinomial model in a single step, this approach
  # fits a sequence of binary logistic equations:
  #   R0 > 1?
  #   R0 > 1 + step?
  #   R0 > 1 + 2*step?
  #   ...
  #   R0 > 10?
  #
  # The cumulative probabilities are converted into bracket probabilities, and
  # the predicted class is the bracket with the largest probability. This is
  # often more stable and easier to interpret for ordered R0 brackets because
  # the R0 categories are naturally ordinal.
  if (isTRUE(use_bayes)) {
    warning("use_bayes is currently not implemented for multiclass brackets; using frequentist successive logistic regression.")
  }

  keep <- c(predictors, outcome_class)

  # Robust numeric R0 source for successive threshold fitting.
  # For generated/uploaded data, the primary numeric column is R0estimate.
  # For built-in ODE simulations, the source is analytic_R0 or numeric_secondary_R0.
  if (outcome_class == "R0_bracket_numeric") {
    if ("R0estimate" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$R0estimate))
    } else if ("numeric_secondary_R0" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$numeric_secondary_R0))
    } else if ("analytic_R0" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$analytic_R0))
    } else {
      r0_values <- rep(NA_real_, nrow(df))
    }
  } else if (outcome_class == "R0_bracket_analytic") {
    if ("analytic_R0" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$analytic_R0))
    } else {
      r0_values <- rep(NA_real_, nrow(df))
    }
  } else {
    # Defensive fallback: use the best available numeric R0-like column.
    if ("R0estimate" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$R0estimate))
    } else if ("numeric_secondary_R0" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$numeric_secondary_R0))
    } else if ("analytic_R0" %in% names(df)) {
      r0_values <- suppressWarnings(as.numeric(df$analytic_R0))
    } else {
      r0_values <- rep(NA_real_, nrow(df))
    }
  }

  d0 <- df[, keep, drop = FALSE]
  d0$R0_for_successive <- r0_values
  d <- d0[stats::complete.cases(d0[, c(predictors, outcome_class, "R0_for_successive"), drop = FALSE]), , drop = FALSE]

  if (nrow(d) < 20) {
    stop(paste0(
      "Too few complete rows for successive binary R0 bracket modeling. ",
      "Rows available after cleaning = ", nrow(d), ". ",
      "Check that R0estimate/numeric_secondary_R0 is finite and that selected predictor columns are numeric and non-missing."
    ))
  }

  k_folds <- as.integer(k_folds)
  if (is.na(k_folds) || k_folds < 1) k_folds <- 1
  k_folds <- min(k_folds, nrow(d))

  bracket_step <- as.numeric(bracket_step)
  if (is.na(bracket_step) || !bracket_step %in% c(0.25, 0.5, 1, 2, 5)) bracket_step <- 0.25
  thresholds <- c(1, seq(1 + bracket_step, 10, by = bracket_step))
  thresholds <- thresholds[thresholds <= 10 + 1e-9]
  if (tail(thresholds, 1) < 10) thresholds <- c(thresholds, 10)
  labels <- c(
    "<1",
    paste0(
      vapply(head(thresholds, -1), format_bracket_number, character(1)),
      "-",
      vapply(tail(thresholds, -1), format_bracket_number, character(1))
    ),
    ">10"
  )

  fit_successive <- function(train_raw) {
    scaled <- standardize_predictors(train_raw, predictors)
    d2 <- scaled$data
    zpred <- scaled$z_predictors
    fits <- list()
    threshold_defaults <- setNames(rep(NA_real_, length(thresholds)), as.character(thresholds))
    usable_thresholds <- c()

    for (thr in thresholds) {
      y <- as.integer(d2$R0_for_successive > thr)
      threshold_defaults[[as.character(thr)]] <- if (all(y == 1, na.rm = TRUE)) 1 else if (all(y == 0, na.rm = TRUE)) 0 else NA_real_
      if (length(unique(y)) < 2) {
        next
      }
      response_name <- paste0("gt_", gsub("\\.", "_", as.character(thr)))
      d2[[response_name]] <- y
      form <- as.formula(paste(response_name, "~", paste(zpred, collapse = " + ")))
      fits[[as.character(thr)]] <- stats::glm(form, data = d2, family = stats::binomial(link = "logit"))
      usable_thresholds <- c(usable_thresholds, thr)
    }

    if (length(fits) < 1) {
      stop("No usable binary threshold models were fitted. Try wider R0 range or more samples.")
    }

    list(
      fit = fits,
      data = d2,
      predictors = predictors,
      z_predictors = zpred,
      centers = scaled$centers,
      scales = scaled$scales,
      thresholds = thresholds,
      threshold_defaults = threshold_defaults,
      usable_thresholds = usable_thresholds,
      labels = labels
    )
  }

  final <- fit_successive(d)
  prob_mat <- successive_probabilities_from_fits(final$fit, final$data, thresholds, labels, threshold_defaults = final$threshold_defaults)
  pred_class <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
  observed_bracket <- make_r0_brackets(final$data$R0_for_successive, step = bracket_step)
  accuracy <- mean(pred_class == as.character(observed_bracket), na.rm = TRUE)

  cv_predictions <- NULL
  if (k_folds >= 2) {
    fold_seed <- if (is.null(seed) || is.na(seed)) 123 else seed
    folds <- make_stratified_folds(make_r0_brackets(d$R0_for_successive, step = bracket_step), k = k_folds, seed = fold_seed)
    cv_rows <- vector("list", k_folds)
    for (fold in seq_len(k_folds)) {
      train_idx <- which(folds != fold)
      test_idx <- which(folds == fold)
      if (length(test_idx) == 0) next

      train_raw <- d[train_idx, , drop = FALSE]
      test_raw <- d[test_idx, , drop = FALSE]
      fold_fit <- tryCatch(fit_successive(train_raw), error = function(e) NULL)
      if (is.null(fold_fit)) next
      test_z <- standardize_newdata(test_raw, fold_fit$centers, fold_fit$scales, predictors)
      fold_prob <- successive_probabilities_from_fits(fold_fit$fit, test_z, thresholds, labels, threshold_defaults = fold_fit$threshold_defaults)
      fold_pred <- colnames(fold_prob)[max.col(fold_prob, ties.method = "first")]
      fold_obs <- make_r0_brackets(test_raw$R0_for_successive, step = bracket_step)
      prob_df <- as.data.frame(fold_prob, check.names = FALSE)
      names(prob_df) <- paste0("prob_", make.names(names(prob_df)))
      cv_rows[[fold]] <- cbind(
        test_raw,
        .row_id = test_idx,
        .fold = fold,
        observed_bracket = as.character(fold_obs),
        predicted_class = fold_pred,
        prob_df
      )
    }
    cv_predictions <- dplyr::bind_rows(cv_rows)
  }

  method <- "Successive binary logistic regression for ordered R0 brackets"

  list(
    fit = final$fit,
    data = final$data,
    predictors = predictors,
    z_predictors = final$z_predictors,
    centers = final$centers,
    scales = final$scales,
    probs = as.data.frame(prob_mat, check.names = FALSE),
    pred_class = pred_class,
    accuracy = accuracy,
    method = method,
    outcome_class = outcome_class,
    class_counts = table(observed_bracket),
    thresholds = thresholds,
    threshold_defaults = final$threshold_defaults,
    usable_thresholds = final$usable_thresholds,
    labels = labels,
    bracket_step = bracket_step,
    observed_bracket = observed_bracket,
    cv_predictions = cv_predictions,
    k_folds = k_folds
  )
}

successive_probabilities_from_fits <- function(fits, newdata, thresholds, labels, threshold_defaults = NULL) {
  n <- nrow(newdata)
  p_gt <- matrix(NA_real_, nrow = n, ncol = length(thresholds))
  colnames(p_gt) <- as.character(thresholds)

  if (is.null(threshold_defaults)) {
    threshold_defaults <- setNames(rep(NA_real_, length(thresholds)), as.character(thresholds))
  }

  for (j in seq_along(thresholds)) {
    thr <- thresholds[j]
    thr_name <- as.character(thr)
    fit <- fits[[thr_name]]
    if (!is.null(fit)) {
      p_gt[, j] <- as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
    } else if (thr_name %in% names(threshold_defaults) && !is.na(threshold_defaults[[thr_name]])) {
      # If a threshold was not fit because all training values were on one side,
      # use the correct constant probability: 0 if no training values exceeded
      # the threshold; 1 if all training values exceeded it.
      p_gt[, j] <- threshold_defaults[[thr_name]]
    }
  }

  # Remaining missing values are filled conservatively from the previous threshold.
  for (j in seq_along(thresholds)) {
    if (all(is.na(p_gt[, j]))) {
      if (j == 1) {
        p_gt[, j] <- 0.5
      } else {
        p_gt[, j] <- p_gt[, j - 1]
      }
    }
  }

  # Enforce non-increasing P(R0 > threshold) as thresholds rise.
  for (j in 2:ncol(p_gt)) {
    p_gt[, j] <- pmin(p_gt[, j - 1], p_gt[, j])
  }

  prob <- matrix(0, nrow = n, ncol = length(labels))
  colnames(prob) <- labels

  prob[, 1] <- 1 - p_gt[, 1]
  for (j in 2:length(thresholds)) {
    prob[, j] <- p_gt[, j - 1] - p_gt[, j]
  }
  prob[, length(labels)] <- p_gt[, length(thresholds)]

  prob[prob < 0] <- 0
  rs <- rowSums(prob)
  rs[rs == 0 | is.na(rs)] <- 1
  prob <- prob / rs
  prob
}

successive_predict_class <- function(mf, znew) {
  prob_mat <- successive_probabilities_from_fits(mf$fit, znew, mf$thresholds, mf$labels, threshold_defaults = mf$threshold_defaults)
  pred_class <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
  list(predicted_class = pred_class, probabilities = as.data.frame(prob_mat, check.names = FALSE))
}

multiclass_equation_text <- function(mf) {
  # Successive binary logistic equations for ordered R0 brackets.
  # Show every requested threshold through R0 > 10. A finite logistic
  # equation is displayed only when observations occur on both sides.
  build_rhs <- function(intercept, betas) {
    rhs <- ""
    if (length(betas) > 0) {
      rhs <- paste(
        paste0(
          ifelse(betas >= 0, " + ", " - "),
          abs(round(as.numeric(betas), 5)),
          "*", names(betas)
        ),
        collapse = ""
      )
    }
    paste0(round(as.numeric(intercept), 5), rhs)
  }

  threshold_names <- as.character(mf$thresholds)
  r0_values <- suppressWarnings(as.numeric(mf$data$R0_for_successive))
  finite_r0 <- r0_values[is.finite(r0_values)]
  range_text <- if (length(finite_r0) > 0) {
    paste0(
      "Observed R0estimate range used for fitting: ",
      round(min(finite_r0), 6), " to ", round(max(finite_r0), 6), "."
    )
  } else {
    "Observed R0estimate range used for fitting: unavailable."
  }

  unavailable_line <- function(thr, scale_label = NULL) {
    n_gt <- sum(r0_values > as.numeric(thr), na.rm = TRUE)
    n_le <- sum(r0_values <= as.numeric(thr), na.rm = TRUE)
    constant_p <- if (n_gt == 0) 0 else if (n_le == 0) 1 else NA_real_
    prefix <- paste0("eta_R0_gt_", format_bracket_number(thr), " = ")
    if (is.finite(constant_p)) {
      paste0(
        prefix,
        "not estimable as a finite logistic equation",
        if (!is.null(scale_label)) paste0(" on the ", scale_label, " scale") else "",
        " (N[R0 <= ", format_bracket_number(thr), "] = ", n_le,
        "; N[R0 > ", format_bracket_number(thr), "] = ", n_gt,
        "; empirical P(R0 > ", format_bracket_number(thr), ") = ", constant_p, ")."
      )
    } else {
      paste0(
        prefix,
        "not estimable",
        if (!is.null(scale_label)) paste0(" on the ", scale_label, " scale") else "",
        " (N[R0 <= ", format_bracket_number(thr), "] = ", n_le,
        "; N[R0 > ", format_bracket_number(thr), "] = ", n_gt, ")."
      )
    }
  }

  lines <- c(
    "SUCCESSIVE BINARY LOGISTIC R0 BRACKET EQUATIONS",
    "",
    range_text,
    "",
    "The app does not use one-time multinomial classification here.",
    "Instead, it fits ordered threshold equations:",
    "  eta_c = log(P(R0 > c) / P(R0 <= c))",
    "  P(R0 > c) = 1 / (1 + exp(-eta_c))",
    "",
    "Prediction rule:",
    "  The fitted cumulative probabilities P(R0 > c) are converted into bracket probabilities.",
    "  The predicted R0 bracket is the bracket with the largest predicted probability.",
    "  This is not a one-time nominal multinomial model; it uses the ordered R0 thresholds.",
    "",
    "Equations on standardized predictor scale:"
  )

  for (thr in threshold_names) {
    fit_thr <- mf$fit[[thr]]
    if (is.null(fit_thr)) {
      lines <- c(lines, unavailable_line(thr, "standardized predictor"))
      next
    }
    co <- stats::coef(fit_thr)
    intercept <- if ("(Intercept)" %in% names(co)) co["(Intercept)"] else co[1]
    betas <- co[setdiff(names(co), "(Intercept)")]
    betas <- betas[names(betas) %in% mf$z_predictors]
    lines <- c(
      lines,
      paste0("eta_R0_gt_", format_bracket_number(thr), " = ", build_rhs(intercept, betas))
    )
  }

  lines <- c(lines, "", "Equations on original ODE/uploaded parameter scale:")

  for (thr in threshold_names) {
    fit_thr <- mf$fit[[thr]]
    if (is.null(fit_thr)) {
      lines <- c(lines, unavailable_line(thr, "original parameter"))
      next
    }
    co <- stats::coef(fit_thr)
    intercept_z <- if ("(Intercept)" %in% names(co)) co["(Intercept)"] else co[1]
    betas_z <- co[setdiff(names(co), "(Intercept)")]
    betas_z <- betas_z[names(betas_z) %in% mf$z_predictors]
    param_names <- sub("^z_", "", names(betas_z))
    valid <- param_names %in% names(mf$centers) & param_names %in% names(mf$scales)
    if (length(betas_z) > 0 && all(valid)) {
      b_original <- as.numeric(betas_z) / as.numeric(mf$scales[param_names])
      names(b_original) <- param_names
      a_original <- as.numeric(intercept_z) -
        sum(as.numeric(betas_z) * as.numeric(mf$centers[param_names]) / as.numeric(mf$scales[param_names]))
      lines <- c(
        lines,
        paste0("eta_R0_gt_", format_bracket_number(thr), " = ", build_rhs(a_original, b_original))
      )
    } else {
      lines <- c(
        lines,
        paste0(
          "eta_R0_gt_", format_bracket_number(thr),
          " = original-scale equation unavailable; check coefficient extraction."
        )
      )
    }
  }

  paste(lines, collapse = "\n")
}

# =========================================================
# LATEX EQUATION EXPORT HELPERS
# =========================================================

latex_escape <- function(x) {
  # Escape text for LaTeX text mode without corrupting LaTeX commands.
  x <- as.character(x)
  x <- gsub("\\\\", "<<BACKSLASH>>", x, perl = TRUE)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("<", "\\\\textless{}", x, fixed = TRUE)
  x <- gsub(">", "$>$", x, fixed = TRUE)
  x <- gsub("<<BACKSLASH>>", "\\\\textbackslash{}", x, fixed = TRUE)
  x
}

latex_math_name <- function(x) {
  # Convert R parameter names into compact valid LaTeX math names.
  # Examples: beta -> \beta, gamma -> \gamma, beta_hv -> \beta_{\mathrm{hv}}, c11 -> \mathrm{c}_{11}.
  x <- as.character(x)
  x <- sub("^z_", "", x)
  greek <- c(
    beta = "\\beta", gamma = "\\gamma", mu = "\\mu", sigma = "\\sigma",
    lambda = "\\lambda", alpha = "\\alpha", delta = "\\delta", rho = "\\rho",
    theta = "\\theta", phi = "\\phi", omega = "\\omega", eta = "\\eta"
  )
  one_name <- function(nm) {
    parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    base <- parts[1]
    suffix <- if (length(parts) > 1) paste(parts[-1], collapse = "_") else ""

    if (base %in% names(greek)) {
      base_tex <- greek[[base]]
    } else if (grepl("^[A-Za-z]+[0-9]+$", base)) {
      base_letters <- sub("^([A-Za-z]+)[0-9]+$", "\\1", base)
      base_digits <- sub("^[A-Za-z]+([0-9]+)$", "\\1", base)
      base_tex <- paste0("\\mathrm{", base_letters, "}_{", base_digits, "}")
    } else if (grepl("^[A-Za-z]+$", base)) {
      base_tex <- paste0("\\mathrm{", base, "}")
    } else {
      base_tex <- paste0("\\mathrm{", gsub("_", "\\\\_", base, fixed = TRUE), "}")
    }

    if (nzchar(suffix)) {
      suffix_tex <- gsub("_", "\\,", suffix, fixed = TRUE)
      paste0(base_tex, "_{\\mathrm{", suffix_tex, "}}")
    } else {
      base_tex
    }
  }
  vapply(x, one_name, character(1), USE.NAMES = FALSE)
}

latex_predictor_math_name <- function(x, standardized = FALSE) {
  # For standardized-scale fitted equations, keep the z-prefix explicit:
  # z_beta is printed as z_{\beta}, not as \beta.
  base_tex <- latex_math_name(x)
  if (isTRUE(standardized)) {
    paste0("z_{", base_tex, "}")
  } else {
    base_tex
  }
}


latex_safe_subscript <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9]+", "", x)
  x
}

format_latex_number <- function(x, digits = 5, math_mode = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) {
    if (isTRUE(math_mode)) return("0")
    return("\\textemdash{}")
  }
  if (abs(x) < 10^(-digits)) x <- 0
  formatC(x, digits = digits, format = "fg", drop0trailing = TRUE)
}

latex_rhs_terms <- function(intercept, betas, digits = 5, standardized = FALSE) {
  terms <- data.frame(
    Sign = "",
    Coefficient = format_latex_number(intercept, digits = digits, math_mode = TRUE),
    Variable = "",
    stringsAsFactors = FALSE
  )
  if (length(betas) > 0) {
    for (nm in names(betas)) {
      b <- as.numeric(betas[[nm]])
      terms <- rbind(
        terms,
        data.frame(
          Sign = ifelse(b >= 0, "+", "-"),
          Coefficient = format_latex_number(abs(b), digits = digits, math_mode = TRUE),
          Variable = latex_predictor_math_name(nm, standardized = standardized),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  terms
}


latex_rhs_inline <- function(intercept, betas, digits = 5, standardized = FALSE) {
  terms <- latex_rhs_terms(intercept, betas, digits = digits, standardized = standardized)
  out <- terms$Coefficient[1]
  if (nrow(terms) > 1) {
    for (i in 2:nrow(terms)) {
      out <- paste0(out, " ", terms$Sign[i], " ", terms$Coefficient[i], "\\,", terms$Variable[i])
    }
  }
  out
}


latex_aligned_equation_lines <- function(lhs, intercept, betas, digits = 5, standardized = FALSE) {
  # This intentionally does not use alignment tabs (&). It is meant to be placed
  # inside a gathered/equation environment so continuation terms do not get pushed
  # far to the right by long probability lines.
  terms <- latex_rhs_terms(intercept, betas, digits = digits, standardized = standardized)
  lines <- paste0(lhs, " = ", terms$Coefficient[1])
  if (nrow(terms) > 1) {
    for (i in 2:nrow(terms)) {
      lines <- c(lines, paste0("{} ", terms$Sign[i], " ", terms$Coefficient[i], "\\,", terms$Variable[i]))
    }
  }
  paste(paste0(lines, c(rep(",\\\\", length(lines) - 1), ".")), collapse = "\n")
}


latex_metric_lines <- function(metrics) {
  if (is.null(metrics) || nrow(metrics) == 0) return("No metrics available.\\\\")
  apply(metrics, 1, function(row) {
    paste0(latex_escape(row[["Validation"]]), " & ",
           latex_escape(row[["Metric"]]), " & ",
           format_latex_number(as.numeric(row[["Value"]]), digits = 4), "\\\\")
  })
}

binary_latex_document <- function(bf, model_name = "REEBLRA model", data_mode = "simulate", cutoff = 1) {
  coef_tab <- posterior_coef_table(bf$fit, z_predictors = bf$z_predictors)
  est <- extract_bayes_coef_means(bf$fit, z_predictors = bf$z_predictors)
  intercept_name <- if ("(Intercept)" %in% names(est)) "(Intercept)" else names(est)[1]
  intercept_z <- unname(est[intercept_name])
  betas_z <- est[setdiff(names(est), intercept_name)]
  betas_z <- betas_z[names(betas_z) %in% bf$z_predictors]

  param_names <- sub("^z_", "", names(betas_z))
  original_available <- length(betas_z) > 0 &&
    all(param_names %in% names(bf$centers)) &&
    all(param_names %in% names(bf$scales))

  if (isTRUE(original_available)) {
    b_original <- as.numeric(betas_z) / as.numeric(bf$scales[param_names])
    names(b_original) <- param_names
    a_original <- as.numeric(intercept_z) -
      sum(as.numeric(betas_z) * as.numeric(bf$centers[param_names]) / as.numeric(bf$scales[param_names]))
  }

  coef_lines <- character(0)
  if (!is.null(coef_tab) && nrow(coef_tab) > 0) {
    coef_lines <- apply(coef_tab, 1, function(row) {
      paste0(latex_escape(row[["Term"]]), " & ",
             format_latex_number(as.numeric(row[["Mean"]]), digits = 5), " & ",
             format_latex_number(as.numeric(row[["SD"]]), digits = 5), " & ",
             format_latex_number(as.numeric(row[["Q2.5"]]), digits = 5), " & ",
             format_latex_number(as.numeric(row[["Q97.5"]]), digits = 5), "\\\\")
    })
  } else {
    coef_lines <- "No coefficient table available & \\textemdash{} & \\textemdash{} & \\textemdash{} & \\textemdash{}\\\\"
  }

  metric_lines <- latex_metric_lines(bf$metrics)
  model_label <- ifelse(identical(data_mode, "upload"), "uploaded R0 table", model_name)
  outcome_text <- ifelse(identical(bf$outcome, "threshold_binary"),
                         paste0("the selected threshold outcome, usually R0 > ", cutoff),
                         bf$outcome)

  body <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{amsmath,amssymb,longtable,array}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage{lmodern}",
    "\\usepackage{microtype}",
    "\\setlength{\\parskip}{0.5em}",
    "\\setlength{\\parindent}{0pt}",
    "\\sloppy",
    "\\emergencystretch=3em",
    "\\title{Binary REEBLRA Logistic Threshold Equation}",
    paste0("\\author{Generated from the REEBLRA Shiny app: ", latex_escape(model_label), "}"),
    "\\date{\\today}",
    "\\begin{document}",
    "\\maketitle",
    "\\section*{Purpose}",
    paste0("This file reports the fitted binary REEBLRA equation for ", latex_escape(outcome_text), ". ",
           "The equation is a simulation-calibrated surrogate threshold model. It approximates threshold behavior in the simulated or uploaded parameter space and should not be read as a replacement for a mechanistic next-generation-matrix derivation."),
    "\\section*{Binary logistic model}",
    "Let \\(Y=1\\) denote that the selected threshold event is true and let \\(p=\\Pr(Y=1\\mid\\boldsymbol{x})\\). REEBLRA fits",
    "\\begin{align}",
    "p &= \\frac{1}{1+\\exp(-\\eta)},\\\\",
    "\\operatorname{logit}(p) &= \\eta.",
    "\\end{align}",
    "\\section*{Equation on standardized predictor scale}",
    "For standardized predictors \\(z_j=(x_j-\\bar{x}_j)/s_j\\), the posterior-mean linear predictor is",
    "\\begin{equation}",
    "\\begin{gathered}",
    latex_aligned_equation_lines("\\eta_z", intercept_z, betas_z, standardized = TRUE),
    "\\end{gathered}",
    "\\end{equation}",
    "and therefore",
    "\\begin{equation}",
    "p_z = \\frac{1}{1+\\exp(-\\eta_z)}.",
    "\\end{equation}",
    "\\section*{Equation on original parameter scale}"
  )

  if (isTRUE(original_available)) {
    body <- c(body,
      "After substituting \\(z_j=(x_j-\\bar{x}_j)/s_j\\), the equation on the original parameter scale is",
      "\\begin{equation}",
      "\\begin{gathered}",
      latex_aligned_equation_lines("\\eta_x", a_original, b_original, standardized = FALSE),
      "\\end{gathered}",
      "\\end{equation}",
      "with",
      "\\begin{equation}",
      "p_x = \\frac{1}{1+\\exp(-\\eta_x)}.",
      "\\end{equation}"
    )
  } else {
    body <- c(body, "The original-scale equation was unavailable because the fitted model did not contain the required scaling information.")
  }

  body <- c(body,
    "\\section*{Classification rule}",
    "Using the default probability cutoff of 0.5,",
    "\\begin{equation}",
    "\\widehat{Y}=\\begin{cases}1, & \\eta>0\\quad (p>0.5),\\\\ 0, & \\eta\\le 0\\quad (p\\le 0.5).\\end{cases}",
    "\\end{equation}",
    "A positive coefficient increases the fitted probability of crossing the threshold, holding the other fitted predictors fixed. A negative coefficient decreases that probability. Interpret coefficient signs carefully when predictors are correlated.",
    "\\section*{Posterior coefficient summary}",
    "\\begin{longtable}{@{}p{0.32\\linewidth}rrrr@{}}",
    "\\hline",
    "Term & Posterior mean & SD & 2.5\\% & 97.5\\%\\\\",
    "\\hline",
    coef_lines,
    "\\hline",
    "\\end{longtable}",
    "\\section*{Model performance summary}",
    "\\begin{longtable}{@{}p{0.45\\linewidth}p{0.25\\linewidth}r@{}}",
    "\\hline",
    "Validation & Metric & Value\\\\",
    "\\hline",
    metric_lines,
    "\\hline",
    "\\end{longtable}",
    "\\section*{Interpretation note}",
    "The binary equation is most useful for screening, sensitivity analysis, and locating the sampled parameter region where the system moves from subthreshold to superthreshold behavior. Strong apparent performance should be compared with k-fold validation because in-sample metrics can be optimistic.",
    "\\end{document}"
  )

  paste(body, collapse = "\n")
}

multiclass_latex_document <- function(mf, model_name = "REEBLRA model", data_mode = "simulate") {
  coef_tab <- multiclass_coef_table(mf)

  threshold_label <- function(thr) {
    paste0("\\(R_0>", latex_escape(format_bracket_number(thr)), "\\)")
  }

  threshold_rows <- split(vapply(mf$thresholds, format_bracket_number, character(1)),
                          ceiling(seq_along(mf$thresholds) / 6))
  threshold_lines <- vapply(threshold_rows, function(x) {
    paste0("\\(", paste(x, collapse = ", "), "\\)\\\\")
  }, character(1))

  make_equation_lines <- function(thr, intercept, betas, suffix = "", standardized = FALSE) {
    sub <- latex_safe_subscript(thr)
    if (nzchar(suffix)) sub <- paste0(sub, suffix)

    term_strings <- character(0)
    if (length(betas) > 0) {
      term_strings <- vapply(seq_along(betas), function(i) {
        b <- as.numeric(betas[[i]])
        paste0(ifelse(b >= 0, "+ ", "- "),
               format_latex_number(abs(b), digits = 5), "\\,",
               latex_predictor_math_name(names(betas)[i], standardized = standardized))
      }, character(1))
    }

    # Break equations across lines without alignment tabs. This prevents the
    # continuation terms from being pushed far to the right by long probability
    # definitions and keeps standardized predictors visibly marked as z_{...}.
    eta_out <- c(
      "\\begin{equation*}",
      "\\begin{gathered}",
      paste0("\\eta_{", sub, "} = ", format_latex_number(intercept, digits = 5))
    )

    if (length(term_strings) > 0) {
      for (trm in term_strings) {
        eta_out <- c(eta_out, paste0("{} ", trm))
      }
    }

    eta_out[length(eta_out)] <- paste0(eta_out[length(eta_out)], ".")
    eta_out <- c(eta_out, "\\end{gathered}", "\\end{equation*}")

    prob_out <- c(
      "\\begin{equation*}",
      paste0("q_{", sub, "} = \\Pr(R_0>", latex_escape(format_bracket_number(thr)),
             "\\,\\mid\\,\\boldsymbol{x}) = \\frac{1}{1+\\exp(-\\eta_{", sub, "})}."),
      "\\end{equation*}"
    )

    c(eta_out, prob_out)
  }

  standardized_eq_lines <- character(0)
  original_eq_lines <- character(0)

  for (thr in names(mf$fit)) {
    co <- stats::coef(mf$fit[[thr]])
    intercept_z <- if ("(Intercept)" %in% names(co)) co["(Intercept)"] else co[1]
    betas_z <- co[setdiff(names(co), "(Intercept)")]
    betas_z <- betas_z[names(betas_z) %in% mf$z_predictors]

    standardized_eq_lines <- c(
      standardized_eq_lines,
      paste0("\\paragraph{Threshold ", threshold_label(thr), ".}"),
      make_equation_lines(thr, intercept_z, betas_z, standardized = TRUE),
      ""
    )

    param_names <- sub("^z_", "", names(betas_z))
    valid <- length(betas_z) > 0 &&
      all(param_names %in% names(mf$centers)) &&
      all(param_names %in% names(mf$scales))

    if (isTRUE(valid)) {
      b_original <- as.numeric(betas_z) / as.numeric(mf$scales[param_names])
      names(b_original) <- param_names
      a_original <- as.numeric(intercept_z) -
        sum(as.numeric(betas_z) * as.numeric(mf$centers[param_names]) / as.numeric(mf$scales[param_names]))
      original_eq_lines <- c(
        original_eq_lines,
        paste0("\\paragraph{Threshold ", threshold_label(thr), ".}"),
        make_equation_lines(thr, a_original, b_original, suffix = "x", standardized = FALSE),
        ""
      )
    }
  }

  coef_lines <- character(0)
  if (!is.null(coef_tab) && nrow(coef_tab) > 0) {
    coef_lines <- apply(coef_tab, 1, function(row) {
      term <- as.character(row[["Term"]])
      term_math <- if (identical(term, "(Intercept)")) "Intercept" else paste0("\\(", latex_predictor_math_name(term, standardized = startsWith(term, "z_")), "\\)")
      thr <- sub("^R0 > ", "", as.character(row[["Threshold"]]))
      paste0(threshold_label(thr), " & ",
             term_math, " & ",
             format_latex_number(as.numeric(row[["Estimate"]]), digits = 5), "\\\\")
    })
  } else {
    coef_lines <- "No coefficient table available & \\textemdash{} & \\textemdash{}\\\\"
  }

  original_rows <- lapply(names(mf$fit), function(thr) {
    co <- stats::coef(mf$fit[[thr]])
    intercept_z <- if ("(Intercept)" %in% names(co)) co["(Intercept)"] else co[1]
    betas_z <- co[setdiff(names(co), "(Intercept)")]
    betas_z <- betas_z[names(betas_z) %in% mf$z_predictors]
    param_names <- sub("^z_", "", names(betas_z))
    valid <- length(betas_z) > 0 &&
      all(param_names %in% names(mf$centers)) &&
      all(param_names %in% names(mf$scales))
    if (isTRUE(valid)) {
      b_original <- as.numeric(betas_z) / as.numeric(mf$scales[param_names])
      names(b_original) <- param_names
      a_original <- as.numeric(intercept_z) -
        sum(as.numeric(betas_z) * as.numeric(mf$centers[param_names]) / as.numeric(mf$scales[param_names]))
      data.frame(
        Threshold = thr,
        Term = c("Intercept", names(b_original)),
        Estimate = c(a_original, as.numeric(b_original)),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(Threshold = thr, Term = "Unavailable", Estimate = NA_real_, stringsAsFactors = FALSE)
    }
  })
  original_tab <- dplyr::bind_rows(original_rows)

  original_lines <- apply(original_tab, 1, function(row) {
    term <- as.character(row[["Term"]])
    term_math <- if (identical(term, "Intercept") || identical(term, "Unavailable")) latex_escape(term) else paste0("\\(", latex_math_name(term), "\\)")
    paste0(threshold_label(row[["Threshold"]]), " & ",
           term_math, " & ",
           format_latex_number(as.numeric(row[["Estimate"]]), digits = 5), "\\\\")
  })

  metric_tab <- multiclass_metrics(mf)
  metric_summary <- metric_tab
  if (!is.null(metric_summary) && nrow(metric_summary) > 0) {
    metric_summary <- metric_summary %>%
      dplyr::filter(Class %in% c("Overall", "Weighted average") |
                      Metric %in% c("Observed support", "One-vs-rest AUC")) %>%
      dplyr::filter(!(Metric == "Observed support" & Support == 0))
  }
  metric_lines <- character(0)
  if (!is.null(metric_summary) && nrow(metric_summary) > 0) {
    metric_lines <- apply(metric_summary, 1, function(row) {
      paste0(latex_escape(row[["Validation"]]), " & ",
             latex_escape(row[["Class"]]), " & ",
             latex_escape(row[["Metric"]]), " & ",
             format_latex_number(as.numeric(row[["Value"]]), digits = 4), " & ",
             format_latex_number(as.numeric(row[["Support"]]), digits = 4), "\\\\")
    })
  } else {
    metric_lines <- "No metrics available & \\textemdash{} & \\textemdash{} & \\textemdash{} & \\textemdash{}\\\\"
  }

  model_label <- ifelse(identical(data_mode, "upload"), "uploaded R0 table", model_name)

  body <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=0.85in]{geometry}",
    "\\usepackage{amsmath,amssymb,longtable,array}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage{lmodern}",
    "\\usepackage{microtype}",
    "\\setlength{\\parskip}{0.5em}",
    "\\setlength{\\parindent}{0pt}",
    "\\sloppy",
    "\\emergencystretch=4em",
    "\\allowdisplaybreaks",
    "\\title{Multiclass REEBLRA Ordered R0-Bracket Equations}",
    paste0("\\author{Generated from the REEBLRA Shiny app: ", latex_escape(model_label), "}"),
    "\\date{\\today}",
    "\\begin{document}",
    "\\maketitle",
    "\\section*{Purpose}",
    "This file reports the multiclass REEBLRA equations for ordered \\(R_0\\)-bracket estimation. The app fits a sequence of binary threshold models rather than a one-time nominal multinomial model.",
    "\\section*{Ordered-threshold model}",
    "The fitted thresholds are:",
    "\\begin{center}",
    "\\begin{tabular}{l}",
    threshold_lines,
    "\\end{tabular}",
    "\\end{center}",
    "For each threshold \\(c_k\\), the model estimates",
    "\\begin{align}",
    "q_k &= \\Pr(R_0>c_k\\mid \\boldsymbol{x}) = \\frac{1}{1+\\exp(-\\eta_k)},\\\\",
    "\\eta_k &= a_k + \\sum_j b_{kj} z_j.",
    "\\end{align}",
    "The fitted cumulative probabilities are forced to be non-increasing as the threshold increases. They are then converted into bracket probabilities.",
    "\\section*{Conversion from cumulative probabilities to bracket probabilities}",
    "For thresholds \\(c_1<c_2<\\cdots<c_K\\),",
    "\\begin{align}",
    "\\Pr(R_0\\le c_1) &= 1-q_1,\\\\",
    "\\Pr(c_{j-1}<R_0\\le c_j) &= q_{j-1}-q_j,\\qquad j=2,\\ldots,K,\\\\",
    "\\Pr(R_0>c_K) &= q_K.",
    "\\end{align}",
    "The predicted bracket is the bracket with the largest predicted bracket probability.",
    "\\section*{Standardized-scale fitted equations}",
    "These equations use standardized predictors \\(z_j=(x_j-\\bar{x}_j)/s_j\\). Each threshold has its own logistic equation. Long equations are line-broken to keep them within the page margins.",
    "\\small",
    standardized_eq_lines,
    "\\normalsize",
    "\\section*{Original-scale fitted equations}",
    "These equations substitute the scaling constants back into the model so that the original ODE or uploaded parameters can be entered directly.",
    "\\small",
    original_eq_lines,
    "\\normalsize",
    "\\section*{Standardized-scale coefficient table}",
    "\\small",
    "\\begin{longtable}{@{}p{0.28\\linewidth}p{0.34\\linewidth}r@{}}",
    "\\hline",
    "Threshold equation & Term & Estimate\\\\",
    "\\hline",
    coef_lines,
    "\\hline",
    "\\end{longtable}",
    "\\normalsize",
    "\\section*{Original-scale coefficient table}",
    "\\small",
    "\\begin{longtable}{@{}p{0.28\\linewidth}p{0.34\\linewidth}r@{}}",
    "\\hline",
    "Threshold equation & Term & Estimate\\\\",
    "\\hline",
    original_lines,
    "\\hline",
    "\\end{longtable}",
    "\\normalsize",
    "\\section*{Model performance summary}",
    "The table reports overall and weighted summaries, plus observed support and one-vs-rest AUC where available. Full per-class sensitivity, specificity, and precision remain available in the CSV downloads.",
    "\\small",
    "\\begin{longtable}{@{}p{0.30\\linewidth}p{0.25\\linewidth}p{0.22\\linewidth}rr@{}}",
    "\\hline",
    "Validation & Class & Metric & Value & Support\\\\",
    "\\hline",
    metric_lines,
    "\\hline",
    "\\end{longtable}",
    "\\normalsize",
    "\\section*{Interpretation note}",
    "The multiclass REEBLRA output should be read as an ordered threshold approximation. Sparse or imbalanced brackets may have low sensitivity or precision, so class distribution and one-vs-rest AUC should be checked before making strong bracket-level claims. Very large coefficients may indicate quasi-separation or sparse brackets; in that case, increase the sample size, use a larger bracket step, or combine sparse brackets.",
    "\\end{document}"
  )

  paste(body, collapse = "\n")
}

multiclass_coef_table <- function(mf) {
  threshold_names <- as.character(mf$thresholds)
  r0_values <- suppressWarnings(as.numeric(mf$data$R0_for_successive))
  expected_terms <- c("(Intercept)", mf$z_predictors)

  rows <- lapply(threshold_names, function(thr) {
    fit_thr <- mf$fit[[thr]]
    n_gt <- sum(r0_values > as.numeric(thr), na.rm = TRUE)
    n_le <- sum(r0_values <= as.numeric(thr), na.rm = TRUE)

    if (is.null(fit_thr)) {
      status <- if (n_gt == 0) {
        "Not estimable: no observations above this threshold"
      } else if (n_le == 0) {
        "Not estimable: no observations at or below this threshold"
      } else {
        "Not estimable"
      }
      return(data.frame(
        Threshold = paste0("R0 > ", format_bracket_number(thr)),
        Term = expected_terms,
        Estimate = NA_real_,
        N_at_or_below = n_le,
        N_above = n_gt,
        Status = status,
        row.names = NULL,
        check.names = FALSE
      ))
    }

    co <- stats::coef(fit_thr)
    estimates <- stats::setNames(rep(NA_real_, length(expected_terms)), expected_terms)
    common_terms <- intersect(expected_terms, names(co))
    estimates[common_terms] <- as.numeric(co[common_terms])

    data.frame(
      Threshold = paste0("R0 > ", format_bracket_number(thr)),
      Term = expected_terms,
      Estimate = as.numeric(estimates),
      N_at_or_below = n_le,
      N_above = n_gt,
      Status = "Fitted",
      row.names = NULL,
      check.names = FALSE
    )
  })

  dplyr::bind_rows(rows)
}


multiclass_coefficient_heatmap_data <- function(mf) {
  threshold_order <- as.character(mf$thresholds)
  if (length(threshold_order) == 0 || length(mf$predictors) == 0) {
    return(data.frame())
  }

  r0_values <- suppressWarnings(as.numeric(mf$data$R0_for_successive))

  rows <- lapply(threshold_order, function(thr) {
    fit_thr <- mf$fit[[thr]]
    n_gt <- sum(r0_values > as.numeric(thr), na.rm = TRUE)
    n_le <- sum(r0_values <= as.numeric(thr), na.rm = TRUE)

    if (is.null(fit_thr)) {
      estimates <- rep(NA_real_, length(mf$predictors))
      status <- if (n_gt == 0) {
        "Not estimable: no observations above this threshold"
      } else if (n_le == 0) {
        "Not estimable: no observations at or below this threshold"
      } else {
        "Not estimable"
      }
    } else {
      co <- stats::coef(fit_thr)
      z_terms <- paste0("z_", mf$predictors)
      estimates <- as.numeric(co[z_terms])
      status <- "Fitted"
    }

    data.frame(
      R0_level = paste0("R0 > ", format_bracket_number(thr)),
      Threshold = suppressWarnings(as.numeric(thr)),
      Parameter = mf$predictors,
      Estimate = estimates,
      N_at_or_below = n_le,
      N_above = n_gt,
      Status = status,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  out <- dplyr::bind_rows(rows) %>%
    dplyr::group_by(R0_level, Threshold) %>%
    dplyr::mutate(
      Row_max_abs = {
        finite_abs <- abs(Estimate[is.finite(Estimate)])
        if (length(finite_abs) == 0) NA_real_ else max(finite_abs)
      },
      Row_normalized = dplyr::if_else(
        is.finite(Estimate) & is.finite(Row_max_abs) & Row_max_abs > 0,
        Estimate / Row_max_abs,
        NA_real_
      )
    ) %>%
    dplyr::ungroup()

  level_order <- paste0(
    "R0 > ",
    vapply(threshold_order, format_bracket_number, character(1))
  )
  out$R0_level <- factor(out$R0_level, levels = rev(level_order))
  out
}


# =========================================================
# DIAGNOSTICS, UPLOAD DATA, MULTICLASS METRICS, AND TRIAL HELPERS
# =========================================================

prepare_uploaded_r0_data <- function(file_path, r0_col = "R0estimate", max_params = 50, bracket_step = 0.25) {
  df <- readr::read_csv(file_path, show_col_types = FALSE)
  if (!r0_col %in% names(df)) {
    stop("Uploaded CSV must contain a column named R0estimate, or select the correct R0 column.")
  }
  df <- as.data.frame(df)
  df[[r0_col]] <- as.numeric(df[[r0_col]])
  candidate_predictors <- setdiff(names(df), r0_col)
  numeric_predictors <- candidate_predictors[sapply(df[candidate_predictors], is.numeric)]
  if (length(numeric_predictors) < 1) {
    stop("Uploaded table must contain at least one numeric parameter column in addition to R0estimate.")
  }
  if (length(numeric_predictors) > max_params) {
    stop(paste0("Uploaded table has ", length(numeric_predictors), " numeric parameter columns. The suggested maximum is ", max_params, ". Please reduce parameters first."))
  }
  out <- df[, c(r0_col, numeric_predictors), drop = FALSE]
  names(out)[names(out) == r0_col] <- "R0estimate"
  out <- out[stats::complete.cases(out), , drop = FALSE]
  if (nrow(out) < 20) stop("Uploaded table has fewer than 20 complete rows after cleaning.")
  out$analytic_R0 <- NA_real_
  out$numeric_secondary_R0 <- out$R0estimate
  out$cases_increase <- as.integer(out$R0estimate > 1)
  out$max_I <- NA_real_
  out$final_size <- NA_real_
  out$R0_bracket_analytic <- NA
  out$R0_bracket_numeric <- make_r0_brackets(out$R0estimate, step = bracket_step)
  out
}

multiclass_metrics_core <- function(observed, predicted, probs, classes, validation_label = "Apparent / in-sample") {
  observed <- factor(as.character(observed), levels = classes)
  predicted <- factor(as.character(predicted), levels = classes)
  probs <- as.data.frame(probs, check.names = FALSE)
  accuracy <- mean(predicted == observed, na.rm = TRUE)

  support_tbl <- table(observed)
  support_df <- tibble::tibble(
    Class = names(support_tbl),
    Support = as.numeric(support_tbl)
  )

  rows <- lapply(classes, function(cls) {
    y_bin <- as.integer(observed == cls)
    pred_bin <- as.integer(predicted == cls)
    TP <- sum(y_bin == 1 & pred_bin == 1, na.rm = TRUE)
    TN <- sum(y_bin == 0 & pred_bin == 0, na.rm = TRUE)
    FP <- sum(y_bin == 0 & pred_bin == 1, na.rm = TRUE)
    FN <- sum(y_bin == 1 & pred_bin == 0, na.rm = TRUE)
    sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
    spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
    prec <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
    auc <- NA_real_
    if (cls %in% names(probs) && length(unique(y_bin)) == 2) {
      auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(y_bin, probs[[cls]], quiet = TRUE))), error = function(e) NA_real_)
    }
    tibble::tibble(
      Validation = validation_label,
      Class = cls,
      Metric = c("Sensitivity", "Specificity", "Precision", "One-vs-rest AUC"),
      Value = c(sens, spec, prec, auc),
      Support = as.numeric(support_tbl[[cls]])
    )
  })

  per_class <- dplyr::bind_rows(rows)

  # Macro average = simple average across classes.
  # This treats rare and common R0 brackets equally, so it is useful for checking
  # whether the model performs reasonably across all brackets, including sparse ones.
  macro <- per_class %>%
    dplyr::group_by(Validation, Metric) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Class = "Macro average", Support = sum(support_df$Support, na.rm = TRUE)) %>%
    dplyr::select(Validation, Class, Metric, Value, Support)

  # Weighted average = average across classes weighted by observed class frequency.
  # This is often more stable when R0 brackets are imbalanced, but it can be dominated
  # by common brackets. Report it together with macro average and per-class metrics.
  weighted <- per_class %>%
    dplyr::group_by(Validation, Metric) %>%
    dplyr::summarise(
      Value = {
        ok <- !is.na(Value) & !is.na(Support) & Support > 0
        if (sum(ok) == 0) NA_real_ else stats::weighted.mean(Value[ok], w = Support[ok])
      },
      Support = sum(Support, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Class = "Weighted average") %>%
    dplyr::select(Validation, Class, Metric, Value, Support)

  support_rows <- support_df %>%
    dplyr::mutate(Validation = validation_label, Metric = "Observed support", Value = Support) %>%
    dplyr::select(Validation, Class, Metric, Value, Support)

  dplyr::bind_rows(
    tibble::tibble(Validation = validation_label, Class = "Overall", Metric = "Accuracy", Value = accuracy, Support = sum(support_df$Support, na.rm = TRUE)),
    macro,
    weighted,
    per_class,
    support_rows
  ) %>%
    dplyr::mutate(Value = round(Value, 4))
}

multiclass_metrics <- function(mf) {
  observed <- if (!is.null(mf$observed_bracket)) {
    factor(as.character(mf$observed_bracket), levels = mf$labels)
  } else {
    factor(mf$data[[mf$outcome_class]])
  }
  probs <- as.data.frame(mf$probs, check.names = FALSE)
  apparent <- multiclass_metrics_core(
    observed = observed,
    predicted = mf$pred_class,
    probs = probs,
    classes = mf$labels,
    validation_label = "Apparent / in-sample"
  )

  if (!is.null(mf$cv_predictions) && nrow(mf$cv_predictions) > 0) {
    cv_prob_names <- paste0("prob_", make.names(mf$labels))
    cv_probs <- mf$cv_predictions[, cv_prob_names, drop = FALSE]
    names(cv_probs) <- mf$labels
    cv <- multiclass_metrics_core(
      observed = mf$cv_predictions$observed_bracket,
      predicted = mf$cv_predictions$predicted_class,
      probs = cv_probs,
      classes = mf$labels,
      validation_label = paste0(mf$k_folds, "-fold cross-validation")
    )
    return(dplyr::bind_rows(cv, apparent))
  }

  apparent
}

mcmc_diagnostics_table <- function(bf) {
  if (is.null(bf) || is.null(bf$fit)) return(data.frame(Note = "No Bayesian fit available."))
  if (!identical(bf$algorithm, "sampling")) {
    return(data.frame(Note = "MCMC diagnostics are not available because the selected algorithm was variational Bayes / meanfield, not MCMC sampling."))
  }

  # Best option: use rstan summary from the underlying Stan fit to retain chain structure.
  diag <- NULL
  if (safe_require("rstan")) {
    diag <- tryCatch({
      sm <- rstan::summary(bf$fit$stanfit)$summary
      out <- as.data.frame(sm)
      out$variable <- rownames(out)
      keep <- grepl("^\\(Intercept\\)$|^z_", out$variable)
      if (any(keep)) out <- out[keep, , drop = FALSE]
      # rstan names are n_eff and Rhat.
      out <- out[, intersect(c("variable", "mean", "sd", "n_eff", "Rhat"), names(out)), drop = FALSE]
      if ("Rhat" %in% names(out)) out$Convergence_flag <- ifelse(is.na(out$Rhat), "Rhat unavailable", ifelse(out$Rhat <= 1.01, "OK", "Check Rhat > 1.01"))
      if ("n_eff" %in% names(out)) out$ESS_flag <- ifelse(is.na(out$n_eff), "ESS unavailable", ifelse(out$n_eff >= 100, "OK", "Low effective sample size"))
      out
    }, error = function(e) NULL)
  }
  if (!is.null(diag)) return(diag)

  # Fallback: posterior package. Depending on the object, Rhat may be unavailable.
  draws <- tryCatch(posterior::as_draws_array(bf$fit), error = function(e) NULL)
  if (is.null(draws)) draws <- tryCatch(posterior::as_draws_matrix(as.matrix(bf$fit)), error = function(e) NULL)
  if (is.null(draws)) return(data.frame(Note = "Could not extract MCMC draws for diagnostics."))
  diag <- tryCatch(posterior::summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail"), error = function(e) NULL)
  if (is.null(diag)) return(data.frame(Note = "Could not compute Rhat/ESS diagnostics."))
  diag <- as.data.frame(diag)
  keep <- grepl("^\\(Intercept\\)$|^z_", diag$variable)
  if (any(keep)) diag <- diag[keep, , drop = FALSE]
  diag$Convergence_flag <- ifelse(is.na(diag$rhat), "Rhat unavailable", ifelse(diag$rhat <= 1.01, "OK", "Check Rhat > 1.01"))
  diag$ESS_flag <- ifelse(is.na(diag$ess_bulk), "ESS unavailable", ifelse(diag$ess_bulk >= 100, "OK", "Low ESS"))
  diag
}

make_trial_inputs <- function(id_prefix, predictors, defaults = NULL) {
  if (length(predictors) == 0) return(helpText("Fit a model first to show parameter inputs."))
  tagList(lapply(predictors, function(p) {
    val <- if (!is.null(defaults) && p %in% names(defaults)) as.numeric(defaults[[p]]) else 0.5
    numericInput(paste0(id_prefix, p), p, value = round(val, 5), step = 0.01)
  }))
}

trial_newdata <- function(input, id_prefix, predictors) {
  vals <- lapply(predictors, function(p) {
    v <- input[[paste0(id_prefix, p)]]
    if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
  })
  names(vals) <- predictors
  as.data.frame(vals, check.names = FALSE)
}

standardize_newdata <- function(newd, centers, scales, predictors) {
  out <- newd
  for (p in predictors) {
    out[[paste0("z_", p)]] <- (as.numeric(out[[p]]) - as.numeric(centers[[p]])) / as.numeric(scales[[p]])
  }
  out
}

predict_binary_trial <- function(bf, newd) {
  znew <- standardize_newdata(newd, bf$centers, bf$scales, bf$predictors)
  prob <- posterior_epred_mean(bf$fit, znew)
  tibble::tibble(
    Predicted_probability_threshold_TRUE = as.numeric(prob),
    Predicted_class = ifelse(prob >= 0.5, "TRUE / above cutoff", "FALSE / below cutoff")
  )
}

predict_multiclass_trial <- function(mf, newd) {
  znew <- standardize_newdata(newd, mf$centers, mf$scales, mf$predictors)
  successive_predict_class(mf, znew)
}

safe_analytic_for_trial <- function(model_name, newd, data_mode = "ode") {
  if (!identical(data_mode, "ode")) return(NA_real_)
  tryCatch(analytic_r0(model_name, newd[1, , drop = FALSE]), error = function(e) NA_real_)
}

format_trial_comparison <- function(pred_text, analytic_value = NA_real_) {
  if (is.na(analytic_value)) {
    paste0(pred_text, "\nAnalytic R0/Re: not available for uploaded datasets or this input.")
  } else {
    paste0(pred_text, "\nAnalytic R0/Re for these parameter values: ", round(analytic_value, 6),
           "\nAnalytic threshold interpretation: ", ifelse(analytic_value > 1, "R0/Re > 1", "R0/Re <= 1"))
  }
}


sensitivity_from_binary_fit <- function(bf) {
  if (is.null(bf)) return(data.frame(Note = "Fit the binary Bayesian logistic model first."))

  draws <- tryCatch(as.matrix(bf$fit), error = function(e) NULL)
  terms <- character(0)
  if (!is.null(draws)) {
    terms <- intersect(bf$z_predictors, colnames(draws))
  }

  if (!is.null(draws) && length(terms) > 0) {
    out <- lapply(terms, function(term) {
      coef_draws <- as.numeric(draws[, term])
      sens_draws <- tanh(coef_draws)
      data.frame(
        Parameter = sub("^z_", "", term),
        Normalized_Bayesian_sensitivity_coefficient_mean = mean(sens_draws, na.rm = TRUE),
        Normalized_Bayesian_sensitivity_coefficient_median = stats::median(sens_draws, na.rm = TRUE),
        Normalized_Bayesian_sensitivity_coefficient_lwr95 = as.numeric(stats::quantile(sens_draws, probs = 0.025, na.rm = TRUE)),
        Normalized_Bayesian_sensitivity_coefficient_upr95 = as.numeric(stats::quantile(sens_draws, probs = 0.975, na.rm = TRUE)),
        Abs_normalized_Bayesian_sensitivity_coefficient_median = abs(stats::median(sens_draws, na.rm = TRUE)),
        Standardized_Bayesian_coefficient_mean = mean(coef_draws, na.rm = TRUE),
        Standardized_Bayesian_coefficient_median = stats::median(coef_draws, na.rm = TRUE),
        Standardized_Bayesian_coefficient_lwr95 = as.numeric(stats::quantile(coef_draws, probs = 0.025, na.rm = TRUE)),
        Standardized_Bayesian_coefficient_upr95 = as.numeric(stats::quantile(coef_draws, probs = 0.975, na.rm = TRUE)),
        Transformation = "tanh(standardized Bayesian coefficient)",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    tab <- dplyr::bind_rows(out)
  } else {
    tab <- posterior_coef_table(bf$fit, z_predictors = bf$z_predictors)
    tab <- tab[tab$Term %in% bf$z_predictors, , drop = FALSE]
    if (nrow(tab) == 0) return(data.frame(Note = "No coefficient terms available."))
    tab$Parameter <- sub("^z_", "", tab$Term)
    tab$Normalized_Bayesian_sensitivity_coefficient_mean <- tanh(tab$Mean)
    tab$Normalized_Bayesian_sensitivity_coefficient_median <- tanh(ifelse(is.na(tab$Median), tab$Mean, tab$Median))
    tab$Normalized_Bayesian_sensitivity_coefficient_lwr95 <- tanh(tab$Q2.5)
    tab$Normalized_Bayesian_sensitivity_coefficient_upr95 <- tanh(tab$Q97.5)
    tab$Abs_normalized_Bayesian_sensitivity_coefficient_median <- abs(tab$Normalized_Bayesian_sensitivity_coefficient_median)
    tab$Standardized_Bayesian_coefficient_mean <- tab$Mean
    tab$Standardized_Bayesian_coefficient_median <- tab$Median
    tab$Standardized_Bayesian_coefficient_lwr95 <- tab$Q2.5
    tab$Standardized_Bayesian_coefficient_upr95 <- tab$Q97.5
    tab$Transformation <- "tanh(standardized Bayesian coefficient)"
    tab <- tab[, c(
      "Parameter",
      "Normalized_Bayesian_sensitivity_coefficient_mean",
      "Normalized_Bayesian_sensitivity_coefficient_median",
      "Normalized_Bayesian_sensitivity_coefficient_lwr95",
      "Normalized_Bayesian_sensitivity_coefficient_upr95",
      "Abs_normalized_Bayesian_sensitivity_coefficient_median",
      "Standardized_Bayesian_coefficient_mean",
      "Standardized_Bayesian_coefficient_median",
      "Standardized_Bayesian_coefficient_lwr95",
      "Standardized_Bayesian_coefficient_upr95",
      "Transformation"
    ), drop = FALSE]
  }

  tab <- tab[order(tab$Abs_normalized_Bayesian_sensitivity_coefficient_median, decreasing = TRUE), , drop = FALSE]
  row.names(tab) <- NULL
  tab
}

collinearity_input_matrix <- function(df, predictors) {
  predictors <- predictors[predictors %in% names(df)]
  predictors <- predictors[sapply(df[predictors], is.numeric)]
  if (length(predictors) == 0) return(NULL)

  x <- df[, predictors, drop = FALSE]
  x <- x[stats::complete.cases(x), , drop = FALSE]
  if (nrow(x) < 3) return(NULL)

  nzv <- sapply(x, function(v) stats::sd(v, na.rm = TRUE) > 0)
  x <- x[, nzv, drop = FALSE]
  if (ncol(x) == 0) return(NULL)
  x
}

vif_table_from_data <- function(df, predictors) {
  x <- collinearity_input_matrix(df, predictors)
  if (is.null(x)) {
    return(data.frame(Note = "No usable numeric parameter columns for VIF."))
  }

  if (ncol(x) < 2) {
    return(data.frame(
      Parameter = names(x),
      VIF = 1,
      Tolerance = 1,
      R2_against_other_parameters = 0,
      Interpretation = "Only one usable parameter; VIF is 1 by definition.",
      row.names = NULL,
      check.names = FALSE
    ))
  }

  out <- lapply(names(x), function(p) {
    others <- setdiff(names(x), p)
    form <- stats::as.formula(paste(p, "~", paste(others, collapse = " + ")))
    fit <- tryCatch(stats::lm(form, data = x), error = function(e) NULL)
    r2 <- if (is.null(fit)) NA_real_ else summary(fit)$r.squared
    vif <- ifelse(is.na(r2), NA_real_, ifelse(r2 >= 1, Inf, 1 / (1 - r2)))
    tol <- ifelse(is.finite(vif) && vif > 0, 1 / vif, ifelse(is.infinite(vif), 0, NA_real_))
    interpretation <- dplyr::case_when(
      is.na(vif) ~ "VIF unavailable.",
      vif < 2.5 ~ "Low collinearity.",
      vif < 5 ~ "Moderate collinearity; usually acceptable, but interpret coefficients with caution.",
      vif < 10 ~ "High collinearity; coefficients may be unstable.",
      TRUE ~ "Very high collinearity; consider removing, combining, or reparameterizing predictors."
    )
    data.frame(
      Parameter = p,
      VIF = vif,
      Tolerance = tol,
      R2_against_other_parameters = r2,
      Interpretation = interpretation,
      row.names = NULL,
      check.names = FALSE
    )
  })

  dplyr::bind_rows(out) %>%
    dplyr::arrange(dplyr::desc(VIF))
}

correlation_matrix_from_data <- function(df, predictors) {
  x <- collinearity_input_matrix(df, predictors)
  if (is.null(x)) return(data.frame(Note = "No usable numeric parameter columns for correlation analysis."))
  stats::cor(x, use = "pairwise.complete.obs")
}

correlation_long_from_data <- function(df, predictors) {
  cor_mat <- correlation_matrix_from_data(df, predictors)
  if (is.data.frame(cor_mat) && "Note" %in% names(cor_mat)) return(cor_mat)
  as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
    dplyr::rename(Parameter_1 = Var1, Parameter_2 = Var2, Correlation = Freq)
}

plot_collinearity_corr_gg <- function(df, predictors) {
  cor_long <- correlation_long_from_data(df, predictors)
  if ("Note" %in% names(cor_long)) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0, label = cor_long$Note[1]) +
             ggplot2::theme_void())
  }
  ggplot2::ggplot(cor_long, ggplot2::aes(x = Parameter_1, y = Parameter_2, fill = Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1), midpoint = 0) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Parameter correlation matrix",
      x = NULL,
      y = NULL,
      fill = "Correlation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_collinearity_vif_gg <- function(df, predictors) {
  vif_tab <- vif_table_from_data(df, predictors)
  if ("Note" %in% names(vif_tab)) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0, label = vif_tab$Note[1]) +
             ggplot2::theme_void())
  }
  vif_tab <- vif_tab %>%
    dplyr::mutate(Parameter = factor(Parameter, levels = Parameter[order(VIF)]))

  ggplot2::ggplot(vif_tab, ggplot2::aes(x = Parameter, y = VIF)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = c(2.5, 5, 10), linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Collinearity check using variance inflation factor (VIF)",
      x = NULL,
      y = "VIF"
    ) +
    ggplot2::theme_minimal()
}



# =========================================================
# DOWNLOAD TABLE AND PLOT HELPERS
# =========================================================

binary_predictions_df <- function(bf) {
  d <- bf$data
  d$predicted_prob <- bf$prob
  d$predicted_class <- as.integer(d$predicted_prob >= 0.5)
  d
}

multiclass_predictions_df <- function(mf) {
  d <- mf$data
  d$observed_bracket <- if (!is.null(mf$observed_bracket)) as.character(mf$observed_bracket) else as.character(mf$data[[mf$outcome_class]])
  d$predicted_class <- mf$pred_class
  probs <- as.data.frame(mf$probs, check.names = FALSE)
  names(probs) <- paste0("prob_", make.names(names(probs)))
  cbind(d, probs)
}

multiclass_confusion_df <- function(mf) {
  observed_raw <- if (!is.null(mf$observed_bracket)) {
    as.character(mf$observed_bracket)
  } else {
    as.character(mf$data[[mf$outcome_class]])
  }
  observed <- factor(observed_raw, levels = mf$labels)
  predicted <- factor(as.character(mf$pred_class), levels = mf$labels)

  as.data.frame.matrix(table(Observed = observed, Predicted = predicted)) %>%
    tibble::rownames_to_column("Observed")
}

multiclass_distribution_df <- function(mf) {
  observed <- if (!is.null(mf$observed_bracket)) {
    factor(as.character(mf$observed_bracket), levels = mf$labels)
  } else {
    factor(as.character(mf$data[[mf$outcome_class]]))
  }
  predicted <- factor(as.character(mf$pred_class), levels = levels(observed))

  dplyr::bind_rows(
    as.data.frame(table(Bracket = observed), stringsAsFactors = FALSE) %>%
      dplyr::mutate(Distribution = "Observed"),
    as.data.frame(table(Bracket = predicted), stringsAsFactors = FALSE) %>%
      dplyr::mutate(Distribution = "Predicted")
  ) %>%
    dplyr::rename(Count = Freq) %>%
    dplyr::mutate(Bracket = as.character(Bracket))
}

multiclass_auc_df <- function(mf, validation_label = "Apparent / in-sample") {
  multiclass_metrics(mf) %>%
    dplyr::filter(
      Metric == "One-vs-rest AUC",
      Validation == validation_label,
      !Class %in% c("Overall", "Macro average", "Weighted average")
    ) %>%
    dplyr::mutate(Class = as.character(Class), Value = as.numeric(Value))
}

multiclass_all_auc_df <- function(mf) {
  multiclass_metrics(mf) %>%
    dplyr::filter(
      Metric == "One-vs-rest AUC",
      !Class %in% c("Overall", "Macro average", "Weighted average")
    ) %>%
    dplyr::mutate(Class = as.character(Class), Value = as.numeric(Value))
}

binary_roc_df <- function(bf) {
  rows <- list()

  roc_app <- tryCatch(
    pROC::roc(bf$data[[bf$outcome]], bf$prob, quiet = TRUE),
    error = function(e) NULL
  )
  if (!is.null(roc_app)) {
    rows[["apparent"]] <- data.frame(
      Validation = "Apparent / in-sample",
      AUC = as.numeric(pROC::auc(roc_app)),
      Specificity = roc_app$specificities,
      Sensitivity = roc_app$sensitivities,
      False_positive_rate = 1 - roc_app$specificities,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(bf$cv_predictions) && nrow(bf$cv_predictions) > 0 &&
      ".observed" %in% names(bf$cv_predictions) &&
      ".cv_predicted_prob" %in% names(bf$cv_predictions) &&
      length(unique(bf$cv_predictions$.observed)) == 2) {
    roc_cv <- tryCatch(
      pROC::roc(bf$cv_predictions$.observed, bf$cv_predictions$.cv_predicted_prob, quiet = TRUE),
      error = function(e) NULL
    )
    if (!is.null(roc_cv)) {
      rows[["cv"]] <- data.frame(
        Validation = paste0(bf$k_folds, "-fold cross-validation"),
        AUC = as.numeric(pROC::auc(roc_cv)),
        Specificity = roc_cv$specificities,
        Sensitivity = roc_cv$sensitivities,
        False_positive_rate = 1 - roc_cv$specificities,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) return(data.frame(Note = "ROC unavailable."))
  dplyr::bind_rows(rows)
}

safe_png <- function(file, plot_expr, width = 1200, height = 900, res = 150) {
  grDevices::png(file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(plot_expr)
}

plot_binary_probability_gg <- function(bf) {
  d <- binary_predictions_df(bf)
  predictors <- bf$z_predictors
  if (length(predictors) >= 2) {
    ggplot2::ggplot(d, ggplot2::aes(x = .data[[predictors[1]]], y = .data[[predictors[2]]], color = predicted_prob)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::labs(
        title = "2D predicted threshold probability scatter",
        x = predictors[1], y = predictors[2], color = "Predicted\nprobability"
      ) +
      ggplot2::theme_minimal()
  } else {
    ggplot2::ggplot(d, ggplot2::aes(x = .data[[predictors[1]]], y = predicted_prob, color = predicted_prob)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::labs(title = "Predicted threshold probability", x = predictors[1], y = "Predicted probability") +
      ggplot2::theme_minimal()
  }
}

plot_binary_roc_gg <- function(bf) {
  roc_df <- binary_roc_df(bf)
  if (!all(c("False_positive_rate", "Sensitivity", "Validation") %in% names(roc_df))) {
    stop("ROC curve unavailable. Check whether both binary classes are present.")
  }
  lab_df <- roc_df %>%
    dplyr::group_by(Validation) %>%
    dplyr::summarise(AUC = dplyr::first(AUC), .groups = "drop") %>%
    dplyr::mutate(Label = paste0(Validation, " (AUC = ", round(AUC, 4), ")"))
  roc_df <- roc_df %>%
    dplyr::left_join(lab_df, by = c("Validation", "AUC"))

  ggplot2::ggplot(roc_df, ggplot2::aes(x = False_positive_rate, y = Sensitivity, linetype = Label)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Binary ROC curve",
      x = "False positive rate",
      y = "True positive rate",
      linetype = "Validation"
    ) +
    ggplot2::theme_minimal()
}

plot_multiclass_prediction_gg <- function(mf) {
  d <- mf$data
  d$predicted_class <- mf$pred_class
  predictors <- mf$z_predictors
  if (length(predictors) >= 2) {
    ggplot2::ggplot(d, ggplot2::aes(x = .data[[predictors[1]]], y = .data[[predictors[2]]], color = predicted_class)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::labs(title = "Predicted R0 bracket", x = predictors[1], y = predictors[2], color = "Predicted\nbracket") +
      ggplot2::theme_minimal()
  } else {
    ggplot2::ggplot(d, ggplot2::aes(x = .data[[predictors[1]]], fill = predicted_class)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.85) +
      ggplot2::labs(title = "Predicted R0 bracket", x = predictors[1], y = "Count", fill = "Predicted\nbracket") +
      ggplot2::theme_minimal()
  }
}

plot_multiclass_distribution_gg <- function(mf) {
  dist_df <- multiclass_distribution_df(mf)
  ggplot2::ggplot(dist_df, ggplot2::aes(x = Bracket, y = Count, fill = Distribution)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(title = "Observed and predicted R0 bracket distribution", x = "R0 bracket", y = "Number of simulations") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_multiclass_auc_gg <- function(mf, validation_label = "Apparent / in-sample") {
  auc_df <- multiclass_auc_df(mf, validation_label = validation_label) %>%
    dplyr::filter(!is.na(Value), is.finite(Value))
  if (nrow(auc_df) == 0) stop("One-vs-rest AUC unavailable for this validation setting.")
  ggplot2::ggplot(auc_df, ggplot2::aes(x = Class, y = Value)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = paste0(validation_label, " one-vs-rest AUC per R0 bracket"), x = "R0 bracket", y = "One-vs-rest AUC") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_sensitivity_gg <- function(bf) {
  s <- sensitivity_from_binary_fit(bf)
  score_col <- "Normalized_Bayesian_sensitivity_coefficient_median"
  lwr_col <- "Normalized_Bayesian_sensitivity_coefficient_lwr95"
  upr_col <- "Normalized_Bayesian_sensitivity_coefficient_upr95"
  if (!score_col %in% names(s)) stop("No sensitivity table available.")
  s$Parameter <- factor(s$Parameter, levels = rev(s$Parameter))
  ggplot2::ggplot(s, ggplot2::aes(x = .data[[score_col]], y = Parameter)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lwr_col]], xmax = .data[[upr_col]]), height = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
    ggplot2::labs(
      title = "Normalized Bayesian sensitivity coefficient",
      x = "Median normalized Bayesian sensitivity coefficient with 95% credible interval",
      y = ""
    ) +
    ggplot2::theme_minimal()
}




# =========================================================
# CSV GENERATOR HELPERS FOR UPLOADED TIME-SERIES DATA
# =========================================================

csvgen_read_uploaded_table <- function(fileinfo) {
  req(fileinfo)
  ext <- tolower(tools::file_ext(fileinfo$name))
  if (ext == "csv") {
    readr::read_csv(fileinfo$datapath, show_col_types = FALSE)
  } else if (ext %in% c("txt", "tsv")) {
    readr::read_tsv(fileinfo$datapath, show_col_types = FALSE)
  } else {
    stop("Unsupported file type. Please upload .csv, .txt, or .tsv.")
  }
}

csvgen_safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

csvgen_fit_growth_one_series <- function(
    time,
    y,
    min_positive_value = 1e-12,
    window_mode = c("automatic", "first_n", "time_range", "below_fraction_of_peak", "all_positive"),
    first_n = 8,
    time_min = NULL,
    time_max = NULL,
    peak_fraction = 0.30,
    min_points = 5
) {
  window_mode <- match.arg(window_mode)

  d <- tibble::tibble(
    t = csvgen_safe_numeric(time),
    y = csvgen_safe_numeric(y)
  ) %>%
    dplyr::filter(is.finite(t), is.finite(y), y > min_positive_value) %>%
    dplyr::arrange(t)

  min_points <- as.integer(min_points)
  first_n <- as.integer(first_n)
  if (!is.finite(min_points) || min_points < 3) min_points <- 5
  if (nrow(d) < min_points) {
    return(tibble::tibble(
      r = NA_real_, r_lwr95 = NA_real_, r_upr95 = NA_real_, growth_R2 = NA_real_,
      n_points_used = nrow(d), window_start = NA_real_, window_end = NA_real_,
      valid_growth_fit = FALSE,
      warning = paste0("Too few positive points. Need at least ", min_points, ".")
    ))
  }

  if (!is.finite(first_n) || first_n < min_points) first_n <- min_points
  first_n <- min(first_n, nrow(d))

  if (window_mode == "automatic") {
    peak_idx <- which.max(d$y)
    d_auto_base <- d[seq_len(max(peak_idx, min_points)), , drop = FALSE]
    if (nrow(d_auto_base) < min_points) d_auto_base <- d
    max_len <- min(first_n, nrow(d_auto_base))
    candidate_lengths <- seq(min_points, max_len)

    score_one <- function(k) {
      dd <- d_auto_base[seq_len(k), , drop = FALSE]
      ff <- tryCatch(stats::lm(log(y) ~ t, data = dd), error = function(e) NULL)
      if (is.null(ff)) return(tibble::tibble(k = k, r = NA_real_, r2 = NA_real_, score = -Inf))
      ss <- summary(ff)
      rr <- unname(stats::coef(ff)[["t"]])
      r2v <- unname(ss$r.squared)
      early_penalty <- 0.001 * max(dd$t, na.rm = TRUE)
      positive_bonus <- ifelse(is.finite(rr) && rr > 0, 1, 0)
      tibble::tibble(k = k, r = rr, r2 = r2v, score = positive_bonus + r2v - early_penalty)
    }

    scores <- dplyr::bind_rows(lapply(candidate_lengths, score_one))
    if (nrow(scores) > 0 && any(is.finite(scores$score))) {
      best_k <- scores$k[which.max(scores$score)]
      d_fit <- d_auto_base[seq_len(best_k), , drop = FALSE]
    } else {
      d_fit <- dplyr::slice_head(d, n = first_n)
    }
  } else if (window_mode == "first_n") {
    d_fit <- dplyr::slice_head(d, n = first_n)
  } else if (window_mode == "time_range") {
    d_fit <- d %>% dplyr::filter(t >= time_min, t <= time_max)
  } else if (window_mode == "below_fraction_of_peak") {
    ymax <- max(d$y, na.rm = TRUE)
    d_fit <- d %>% dplyr::filter(y <= peak_fraction * ymax)
    if (nrow(d_fit) > first_n) d_fit <- dplyr::slice_head(d_fit, n = first_n)
  } else {
    d_fit <- d
  }

  if (nrow(d_fit) < min_points) d_fit <- dplyr::slice_head(d, n = first_n)

  if (nrow(d_fit) < min_points) {
    return(tibble::tibble(
      r = NA_real_, r_lwr95 = NA_real_, r_upr95 = NA_real_, growth_R2 = NA_real_,
      n_points_used = nrow(d_fit),
      window_start = ifelse(nrow(d_fit) > 0, min(d_fit$t), NA_real_),
      window_end = ifelse(nrow(d_fit) > 0, max(d_fit$t), NA_real_),
      valid_growth_fit = FALSE,
      warning = paste0("Too few points in selected early-growth window. Need at least ", min_points, ".")
    ))
  }

  fit <- tryCatch(stats::lm(log(y) ~ t, data = d_fit), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(
      r = NA_real_, r_lwr95 = NA_real_, r_upr95 = NA_real_, growth_R2 = NA_real_,
      n_points_used = nrow(d_fit), window_start = min(d_fit$t), window_end = max(d_fit$t),
      valid_growth_fit = FALSE, warning = "Log-linear growth-rate fit failed."
    ))
  }

  sm <- summary(fit)
  r_hat <- unname(stats::coef(fit)[["t"]])
  r2 <- unname(sm$r.squared)
  ci <- tryCatch(stats::confint(fit, "t", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))

  warn <- character(0)
  if (!is.finite(r_hat)) warn <- c(warn, "Estimated r is not finite.")
  if (is.finite(r_hat) && r_hat <= 0) warn <- c(warn, "Estimated r is non-positive; epidemic growth was not detected in the chosen window.")
  if (is.finite(ci[1, 1]) && is.finite(ci[1, 2]) && ci[1, 1] <= 0 && ci[1, 2] >= 0) warn <- c(warn, "The 95% CI for r crosses 0.")
  if (is.finite(r2) && r2 < 0.80) warn <- c(warn, "Low R-squared for log-linear early-growth fit.")

  tibble::tibble(
    r = r_hat,
    r_lwr95 = as.numeric(ci[1, 1]),
    r_upr95 = as.numeric(ci[1, 2]),
    growth_R2 = r2,
    n_points_used = nrow(d_fit),
    window_start = min(d_fit$t),
    window_end = max(d_fit$t),
    valid_growth_fit = is.finite(r_hat) && r_hat > 0,
    warning = ifelse(length(warn) == 0, "OK", paste(unique(warn), collapse = " "))
  )
}

csvgen_convert_r_to_R0 <- function(
    r,
    conversion_method,
    gamma_removal = NA_real_,
    sigma_latency = NA_real_,
    generation_time = NA_real_,
    nonpositive_rule = c("zero", "one", "keep_formula", "missing")
) {
  nonpositive_rule <- match.arg(nonpositive_rule)
  if (!is.finite(r)) return(NA_real_)
  if (r <= 0) {
    if (nonpositive_rule == "zero") return(0)
    if (nonpositive_rule == "one") return(1)
    if (nonpositive_rule == "missing") return(NA_real_)
  }
  out <- NA_real_
  if (conversion_method == "SIR-like: R0 = 1 + r/gamma_removal") {
    if (is.finite(gamma_removal) && gamma_removal > 0) out <- 1 + r / gamma_removal
  } else if (conversion_method == "SEIR-like with latency: R0 = ((r+sigma)(r+gamma))/(sigma*gamma)") {
    if (is.finite(gamma_removal) && gamma_removal > 0 && is.finite(sigma_latency) && sigma_latency > 0) {
      out <- ((r + sigma_latency) * (r + gamma_removal)) / (sigma_latency * gamma_removal)
    }
  } else if (conversion_method == "Generation-time approximation: R0 = exp(r*Tgen)") {
    if (is.finite(generation_time) && generation_time > 0) out <- exp(r * generation_time)
  } else if (conversion_method == "Unit-time Euler multiplier: R0estimate = 1 + r (one time unit)") {
    # This preserves the epidemic threshold r > 0 <=> R0estimate > 1,
    # but should be interpreted as a one-time-unit Euler growth multiplier,
    # not as a mechanistic generation-interval conversion from r to R0.
    out <- 1 + r
  }
  if (is.finite(out) && out < 1 && r <= 0 && nonpositive_rule == "zero") return(0)
  out
}


# =========================================================
# CUSTOM USER-DEFINED ODE HELPERS
# =========================================================

custom_clean_lines <- function(txt) {
  lines <- unlist(strsplit(as.character(txt), "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!grepl("^#", lines)]
  lines
}

custom_parse_named_numeric_lines <- function(txt, value_name = "value") {
  lines <- custom_clean_lines(txt)
  if (length(lines) == 0) {
    stop(paste0("No ", value_name, " entries were provided."))
  }

  rows <- lapply(lines, function(line) {
    # Accept name = value, name: value, name, value, or name value.
    parts <- unlist(strsplit(line, "=|:|,|\\s+"))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) < 2) {
      stop(paste0("Could not parse line: ", line))
    }
    nm <- parts[1]
    val <- suppressWarnings(as.numeric(parts[2]))
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", nm)) {
      stop(paste0("Invalid variable/parameter name: ", nm, ". Use letters, numbers, and underscores; start with a letter."))
    }
    if (!is.finite(val)) {
      stop(paste0("Invalid numeric ", value_name, " for ", nm, "."))
    }
    data.frame(name = nm, value = val, stringsAsFactors = FALSE)
  })

  out <- dplyr::bind_rows(rows)
  if (anyDuplicated(out$name)) {
    stop(paste0("Duplicate name found in ", value_name, " entries."))
  }
  stats::setNames(out$value, out$name)
}

custom_parse_parameter_ranges <- function(txt, max_params = 50) {
  lines <- custom_clean_lines(txt)
  if (length(lines) == 0) stop("No parameter ranges were provided.")

  rows <- lapply(lines, function(line) {
    # Accept a = 0.01, 1; a: 0.01, 1; a 0.01 1
    parts <- unlist(strsplit(line, "=|:|,|\\s+"))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) < 3) {
      stop(paste0("Could not parse parameter range line: ", line, ". Use format such as a = 0.01, 1."))
    }
    nm <- parts[1]
    mn <- suppressWarnings(as.numeric(parts[2]))
    mx <- suppressWarnings(as.numeric(parts[3]))
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", nm)) {
      stop(paste0("Invalid parameter name: ", nm, ". Use letters, numbers, and underscores; start with a letter."))
    }
    if (!is.finite(mn) || !is.finite(mx)) {
      stop(paste0("Invalid numeric range for parameter ", nm, "."))
    }
    if (mn >= mx) {
      stop(paste0("Parameter ", nm, " must have min < max."))
    }
    data.frame(parameter = nm, min = mn, max = mx, stringsAsFactors = FALSE)
  })

  out <- dplyr::bind_rows(rows)
  if (anyDuplicated(out$parameter)) stop("Duplicate parameter name found.")
  if (nrow(out) > max_params) {
    stop(paste0("Maximum of ", max_params, " parameters allowed for custom ODE simulation."))
  }
  out
}

custom_parse_equations <- function(eq_text, max_equations = 30) {
  lines <- custom_clean_lines(eq_text)
  if (length(lines) == 0) stop("No custom ODE equations were provided.")
  if (length(lines) > max_equations) {
    stop(paste0("Maximum of ", max_equations, " state equations allowed."))
  }
  if (any(!grepl("=", lines, fixed = TRUE))) {
    stop("Each ODE line must contain '='. Example: dx/dt = a*x - b*x*y")
  }

  lhs <- trimws(sub("=.*$", "", lines))
  rhs <- trimws(sub("^[^=]*=", "", lines))
  if (any(!nzchar(rhs))) stop("Each ODE equation must have a non-empty right-hand side.")

  states <- lhs
  states <- gsub("\\s+", "", states)
  states <- sub("^d\\(([^)]+)\\)/dt$", "\\1", states)
  states <- sub("^d([A-Za-z][A-Za-z0-9_]*)/dt$", "\\1", states)
  states <- sub("^d_?([A-Za-z][A-Za-z0-9_]*)$", "\\1", states)

  if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", states))) {
    stop("Invalid state name detected. Use forms such as dx/dt, dI/dt, dx, dI, or I.")
  }
  if (anyDuplicated(states)) stop("Duplicate state equation detected.")

  parsed_rhs <- lapply(rhs, function(z) {
    tryCatch(parse(text = z)[[1]], error = function(e) {
      stop(paste0("Could not parse right-hand side: ", z, ". Error: ", e$message))
    })
  })

  list(lines = lines, states = states, rhs = rhs, parsed_rhs = parsed_rhs)
}

custom_parse_auxiliary_equations <- function(eq_text, max_equations = 30) {
  lines <- custom_clean_lines(eq_text)
  if (length(lines) == 0) {
    return(list(lines = character(0), variables = character(0), rhs = character(0), parsed_rhs = list()))
  }
  if (length(lines) > max_equations) {
    stop(paste0("Maximum of ", max_equations, " auxiliary equations allowed."))
  }
  if (any(!grepl("=", lines, fixed = TRUE))) {
    stop("Each non-ODE equation line must contain '='. Example: I_total = I1 + I2")
  }

  lhs <- trimws(sub("=.*$", "", lines))
  rhs <- trimws(sub("^[^=]*=", "", lines))
  lhs <- gsub("\\s+", "", lhs)
  if (any(!nzchar(rhs))) stop("Each non-ODE equation must have a non-empty right-hand side.")
  if (any(grepl("^d", lhs))) {
    stop("Non-ODE equations should define auxiliary/output variables, not derivatives. Example: I_total = I1 + I2")
  }
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", lhs))) {
    stop("Invalid auxiliary/output variable name detected. Use letters, numbers, and underscores; start with a letter.")
  }
  if (anyDuplicated(lhs)) stop("Duplicate auxiliary/output equation detected.")

  parsed_rhs <- lapply(rhs, function(z) {
    tryCatch(parse(text = z)[[1]], error = function(e) {
      stop(paste0("Could not parse non-ODE right-hand side: ", z, ". Error: ", e$message))
    })
  })

  list(lines = lines, variables = lhs, rhs = rhs, parsed_rhs = parsed_rhs)
}

custom_eval_auxiliary_values <- function(aux_info, state_values, parameters, time) {
  if (is.null(aux_info) || length(aux_info$variables) == 0) return(numeric(0))
  env <- new.env(parent = baseenv())
  for (nm in names(state_values)) assign(nm, as.numeric(state_values[[nm]]), envir = env)
  for (nm in names(parameters)) assign(nm, as.numeric(parameters[[nm]]), envir = env)
  assign("time", time, envir = env)
  assign("t", time, envir = env)
  if (!exists("N", envir = env, inherits = FALSE)) {
    assign("N", sum(as.numeric(state_values), na.rm = TRUE), envir = env)
  }

  vals <- numeric(length(aux_info$variables))
  names(vals) <- aux_info$variables
  for (i in seq_along(aux_info$variables)) {
    val <- tryCatch(eval(aux_info$parsed_rhs[[i]], envir = env), error = function(e) NA_real_)
    vals[i] <- as.numeric(val)[1]
    assign(aux_info$variables[i], vals[i], envir = env)
  }
  if (any(!is.finite(vals))) {
    stop("Custom non-ODE equation returned a non-finite value. Check auxiliary equations, parameter ranges, and denominators.")
  }
  vals
}

custom_add_auxiliary_outputs <- function(sim, aux_info, par_row) {
  sim <- as.data.frame(sim)
  if (is.null(aux_info) || length(aux_info$variables) == 0) return(sim)
  aux_mat <- matrix(NA_real_, nrow = nrow(sim), ncol = length(aux_info$variables))
  colnames(aux_mat) <- aux_info$variables
  state_cols <- setdiff(names(sim), "time")
  for (i in seq_len(nrow(sim))) {
    aux_mat[i, ] <- custom_eval_auxiliary_values(
      aux_info = aux_info,
      state_values = stats::setNames(as.numeric(sim[i, state_cols, drop = TRUE]), state_cols),
      parameters = as.list(par_row),
      time = sim$time[i]
    )
  }
  dplyr::bind_cols(sim, as.data.frame(aux_mat, check.names = FALSE))
}

custom_enforce_nonnegative_output <- function(sim, nonnegative_states = character(0)) {
  sim <- as.data.frame(sim)
  nonnegative_states <- intersect(as.character(nonnegative_states), names(sim))
  nonnegative_states <- setdiff(nonnegative_states, "time")
  if (length(nonnegative_states) == 0) return(sim)

  for (nm in nonnegative_states) {
    sim[[nm]] <- pmax(suppressWarnings(as.numeric(sim[[nm]])), 0)
  }
  sim
}

custom_validate_model_inputs <- function(eq_info, aux_info, init_state, param_table) {
  missing_init <- setdiff(eq_info$states, names(init_state))
  if (length(missing_init) > 0) {
    stop(paste0("Missing initial values for state(s): ", paste(missing_init, collapse = ", ")))
  }
  init_state <- init_state[eq_info$states]

  if (!is.null(aux_info) && length(aux_info$variables) > 0) {
    duplicated_model_names <- intersect(aux_info$variables, c(eq_info$states, param_table$parameter))
    if (length(duplicated_model_names) > 0) {
      stop(paste0("Auxiliary/output variable(s) cannot duplicate state or parameter names: ", paste(duplicated_model_names, collapse = ", ")))
    }
  }

  allowed_special <- c("time", "t", "N", "pi")
  aux_variables <- if (is.null(aux_info)) character(0) else aux_info$variables
  rhs_symbols <- unique(unlist(lapply(eq_info$parsed_rhs, all.vars)))
  unknown <- setdiff(rhs_symbols, c(eq_info$states, param_table$parameter, aux_variables, allowed_special))
  if (length(unknown) > 0) {
    stop(paste0(
      "Unknown symbol(s) in ODE equations: ", paste(unknown, collapse = ", "),
      ". Declare them as states, parameters, or non-ODE auxiliary equations."
    ))
  }

  available <- c(eq_info$states, param_table$parameter, allowed_special)
  if (!is.null(aux_info) && length(aux_info$variables) > 0) {
    for (i in seq_along(aux_info$variables)) {
      aux_symbols <- all.vars(aux_info$parsed_rhs[[i]])
      unknown_aux <- setdiff(aux_symbols, available)
      if (length(unknown_aux) > 0) {
        stop(paste0(
          "Unknown symbol(s) in non-ODE equation for ", aux_info$variables[i], ": ",
          paste(unknown_aux, collapse = ", "),
          ". Auxiliary equations are evaluated from top to bottom, so define dependencies earlier."
        ))
      }
      available <- c(available, aux_info$variables[i])
    }
  }

  list(eq_info = eq_info, aux_info = aux_info, init_state = init_state, param_table = param_table)
}

make_custom_ode_function <- function(eq_info, aux_info = NULL, nonnegative_states = character(0)) {
  states <- eq_info$states
  parsed_rhs <- eq_info$parsed_rhs
  nonnegative_states <- intersect(as.character(nonnegative_states), states)

  function(time, state, parameters) {
    env <- new.env(parent = baseenv())

    raw_state <- as.numeric(state)
    names(raw_state) <- names(state)
    eval_state <- raw_state
    if (length(nonnegative_states) > 0) {
      nn <- intersect(nonnegative_states, names(eval_state))
      eval_state[nn] <- pmax(eval_state[nn], 0)
    }

    for (nm in names(eval_state)) assign(nm, as.numeric(eval_state[[nm]]), envir = env)
    for (nm in names(parameters)) assign(nm, as.numeric(parameters[[nm]]), envir = env)
    assign("time", time, envir = env)
    assign("t", time, envir = env)
    if (!exists("N", envir = env, inherits = FALSE)) {
      assign("N", sum(as.numeric(eval_state), na.rm = TRUE), envir = env)
    }
    if (!is.null(aux_info) && length(aux_info$variables) > 0) {
      for (i in seq_along(aux_info$variables)) {
        val <- tryCatch(eval(aux_info$parsed_rhs[[i]], envir = env), error = function(e) NA_real_)
        val <- as.numeric(val)[1]
        if (!is.finite(val)) stop("Custom non-ODE equation returned a non-finite value during ODE simulation.")
        assign(aux_info$variables[i], val, envir = env)
      }
    }

    derivs <- numeric(length(states))
    names(derivs) <- states
    for (i in seq_along(states)) {
      val <- tryCatch(eval(parsed_rhs[[i]], envir = env), error = function(e) NA_real_)
      derivs[i] <- as.numeric(val)[1]
    }
    if (any(!is.finite(derivs))) {
      stop("Custom ODE returned a non-finite derivative. Check equations, parameter ranges, and denominators.")
    }

    if (length(nonnegative_states) > 0) {
      nn <- intersect(nonnegative_states, names(derivs))
      at_lower_bound <- is.finite(raw_state[nn]) & raw_state[nn] <= 0 & derivs[nn] < 0
      if (any(at_lower_bound, na.rm = TRUE)) {
        derivs[nn[at_lower_bound]] <- 0
      }
    }

    list(derivs)
  }
}

simulate_custom_one <- function(par_row, eq_info, aux_info, init_state, times, target_state,
                                nonnegative_states = character(0),
                                growth_epsilon = 1e-5,
                                conversion_method,
                                gamma_removal = NA_real_,
                                sigma_latency = NA_real_,
                                generation_time = NA_real_,
                                nonpositive_rule = "zero",
                                min_positive_value = 1e-12,
                                growth_window_first_n = 8,
                                min_points = 5) {
  nonnegative_states <- intersect(as.character(nonnegative_states), eq_info$states)
  ode_fun <- make_custom_ode_function(eq_info, aux_info = aux_info, nonnegative_states = nonnegative_states)
  sim <- deSolve::ode(
    y = init_state,
    times = times,
    func = ode_fun,
    parms = as.list(par_row),
    method = "lsoda"
  )
  sim <- custom_enforce_nonnegative_output(sim, nonnegative_states = nonnegative_states)
  sim <- custom_add_auxiliary_outputs(sim, aux_info, par_row)

  if (!target_state %in% names(sim)) {
    stop(paste0("Target infectious/output state '", target_state, "' is not a simulated state or auxiliary/output variable."))
  }

  y <- suppressWarnings(as.numeric(sim[[target_state]]))
  initial_y <- y[1]
  max_y <- suppressWarnings(max(y, na.rm = TRUE))
  cases_increase <- as.integer(is.finite(max_y) && is.finite(initial_y) && max_y > initial_y + growth_epsilon)

  fit_diag <- csvgen_fit_growth_one_series(
    time = sim$time,
    y = y,
    min_positive_value = min_positive_value,
    window_mode = "first_n",
    first_n = growth_window_first_n,
    min_points = min_points
  )

  r0_est <- csvgen_convert_r_to_R0(
    r = fit_diag$r,
    conversion_method = conversion_method,
    gamma_removal = gamma_removal,
    sigma_latency = sigma_latency,
    generation_time = generation_time,
    nonpositive_rule = nonpositive_rule
  )

  tibble::as_tibble(par_row) %>%
    dplyr::mutate(
      R0estimate = r0_est,
      analytic_R0 = NA_real_,
      numeric_secondary_R0 = R0estimate,
      cases_increase = cases_increase,
      max_I = max_y,
      final_size = dplyr::last(y),
      early_growth_r = fit_diag$r,
      early_growth_r_lwr95 = fit_diag$r_lwr95,
      early_growth_r_upr95 = fit_diag$r_upr95,
      early_growth_R2 = fit_diag$growth_R2,
      early_growth_n = fit_diag$n_points_used,
      early_growth_threshold_used = NA_real_,
      early_growth_window_end = fit_diag$window_end,
      early_growth_valid = fit_diag$valid_growth_fit,
      early_growth_warning = fit_diag$warning
    )
}

simulate_custom_many <- function(param_table, eq_info, aux_info, init_state, n, tmax, dt, target_state,
                                 nonnegative_states = character(0),
                                 growth_epsilon = 1e-5,
                                 conversion_method,
                                 gamma_removal = NA_real_,
                                 sigma_latency = NA_real_,
                                 generation_time = NA_real_,
                                 nonpositive_rule = "zero",
                                 min_positive_value = 1e-12,
                                 growth_window_first_n = 8,
                                 min_points = 5,
                                 show_progress = FALSE) {
  samples <- latin_hypercube_samples(param_table, n)
  times <- seq(0, tmax, by = dt)
  out <- vector("list", n)

  for (i in seq_len(n)) {
    if (isTRUE(show_progress)) {
      shiny::incProgress(1 / n, detail = paste("Custom ODE simulation", i, "of", n))
    }
    out[[i]] <- tryCatch(
      simulate_custom_one(
        par_row = samples[i, , drop = FALSE],
        eq_info = eq_info,
        aux_info = aux_info,
        init_state = init_state,
        times = times,
        target_state = target_state,
        nonnegative_states = nonnegative_states,
        growth_epsilon = growth_epsilon,
        conversion_method = conversion_method,
        gamma_removal = gamma_removal,
        sigma_latency = sigma_latency,
        generation_time = generation_time,
        nonpositive_rule = nonpositive_rule,
        min_positive_value = min_positive_value,
        growth_window_first_n = growth_window_first_n,
        min_points = min_points
      ),
      error = function(e) {
        tibble::as_tibble(samples[i, , drop = FALSE]) %>%
          dplyr::mutate(
            R0estimate = NA_real_, analytic_R0 = NA_real_, numeric_secondary_R0 = NA_real_,
            cases_increase = NA_integer_, max_I = NA_real_, final_size = NA_real_,
            early_growth_r = NA_real_, early_growth_r_lwr95 = NA_real_, early_growth_r_upr95 = NA_real_,
            early_growth_R2 = NA_real_, early_growth_n = NA_integer_, early_growth_threshold_used = NA_real_,
            early_growth_window_end = NA_real_, early_growth_valid = FALSE,
            early_growth_warning = paste0("Custom simulation or growth fitting failed: ", e$message)
          )
      }
    )
  }

  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      R0_bracket_analytic = NA,
      R0_bracket_numeric = make_r0_brackets(numeric_secondary_R0)
    )
}

custom_ode_equation_text <- function(eq_info = NULL, aux_info = NULL, target_state = NULL, conversion_method = NULL, nonnegative_states = character(0)) {
  if (is.null(eq_info)) {
    return("Custom ODE mode: enter equations, optional non-ODE auxiliary/output equations, initial values, parameter ranges, and choose the target infectious/output state. Analytic next-generation-matrix R0 is not derived automatically.")
  }
  aux_lines <- if (!is.null(aux_info) && length(aux_info$lines) > 0) paste(aux_info$lines, collapse = "\n") else "None"
  nonnegative_states <- intersect(as.character(nonnegative_states), eq_info$states)
  nonnegative_line <- if (length(nonnegative_states) > 0) paste(nonnegative_states, collapse = ", ") else "None"
  paste(
    "Custom user-defined ODE model:",
    paste(eq_info$lines, collapse = "\n"),
    "",
    "Custom non-ODE auxiliary/output equations:",
    aux_lines,
    "",
    paste0("Target infectious/output state for early growth: ", target_state),
    paste0("State variables forced to be nonnegative: ", nonnegative_line),
    paste0("Growth-rate-to-R0estimate conversion: ", conversion_method),
    "Note: if one species has multiple infectious classes, define an auxiliary total such as I_total = I1 + I2 and use that total as the target infectious/output state for growth fitting.",
    "Analytic R0/Re is not available in custom mode unless supplied externally.",
    sep = "\n"
  )
}

# =========================================================
# SHINY UI
# =========================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* REEBLRA_STYLE_VERSION_2026_06_11 */
      :root {
        --reeblra-maroon: #6f1726;
        --reeblra-maroon-dark: #4f0f1a;
        --reeblra-green: #0f6b3d;
        --reeblra-green-dark: #084b2b;
        --reeblra-gold: #d6a326;
        --reeblra-gold-soft: #fff3c7;
        --reeblra-blue: #1f5f99;
        --reeblra-blue-soft: #eaf4ff;
        --reeblra-cream: #fffaf0;
        --reeblra-panel: #ffffff;
        --reeblra-ink: #24302a;
        --reeblra-muted: #64706a;
        --reeblra-border: rgba(111, 23, 38, 0.16);
        --reeblra-shadow: 0 14px 34px rgba(79, 15, 26, 0.12);
      }

      body {
        background:
          radial-gradient(circle at top left, rgba(214, 163, 38, 0.20), transparent 28rem),
          radial-gradient(circle at top right, rgba(31, 95, 153, 0.16), transparent 30rem),
          linear-gradient(180deg, #fffaf0 0%, #f7fbf7 46%, #f6f8fb 100%);
        color: var(--reeblra-ink);
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
      }

      .container-fluid {
        max-width: 1720px;
        padding: 18px 24px 36px 24px;
      }

      .container-fluid > h2:first-child {
        margin: 0 0 14px 0;
        padding: 22px 26px;
        color: #ffffff;
        font-weight: 800;
        letter-spacing: 0.2px;
        line-height: 1.22;
        border-radius: 24px;
        background:
          linear-gradient(125deg, var(--reeblra-maroon-dark) 0%, var(--reeblra-maroon) 46%, var(--reeblra-green) 100%);
        box-shadow: var(--reeblra-shadow);
        border-bottom: 5px solid var(--reeblra-gold);
      }

      .app-header {
        position: relative;
        margin: 0 0 18px 0;
        padding: 22px 26px 20px 26px;
        border-radius: 24px;
        background:
          linear-gradient(135deg, rgba(255, 255, 255, 0.96), rgba(255, 250, 240, 0.93)),
          linear-gradient(135deg, rgba(111, 23, 38, 0.06), rgba(15, 107, 61, 0.07));
        border: 1px solid var(--reeblra-border);
        border-left: 8px solid var(--reeblra-maroon);
        box-shadow: var(--reeblra-shadow);
        overflow: hidden;
      }

      .app-header:after {
        content: '';
        position: absolute;
        right: -44px;
        top: -44px;
        width: 180px;
        height: 180px;
        border-radius: 999px;
        background: rgba(214, 163, 38, 0.18);
        box-shadow: -54px 54px 0 rgba(31, 95, 153, 0.10);
        pointer-events: none;
      }

      .app-header h5,
      .app-header h4 {
        position: relative;
        z-index: 1;
      }

      .app-header h5 {
        margin: 4px 0;
        color: var(--reeblra-maroon) !important;
        font-weight: 700;
      }

      .app-header h4:first-of-type {
        margin-top: 10px;
        color: var(--reeblra-green) !important;
        font-weight: 800 !important;
      }

      .app-header h4:not(:first-of-type) {
        color: var(--reeblra-ink);
        font-size: 16px;
        line-height: 1.55;
        font-weight: 500;
      }

      .app-header h5:last-child {
        margin-top: 14px;
        padding: 12px 14px;
        color: var(--reeblra-maroon-dark) !important;
        background: linear-gradient(90deg, var(--reeblra-gold-soft), rgba(234, 244, 255, 0.72));
        border-left: 5px solid var(--reeblra-gold);
        border-radius: 14px;
      }

      .well {
        background: rgba(255, 255, 255, 0.95);
        border: 1px solid var(--reeblra-border);
        border-radius: 22px;
        box-shadow: var(--reeblra-shadow);
      }

      .well h4,
      .well h5 {
        color: var(--reeblra-maroon);
        font-weight: 800;
      }

      .well h4 {
        padding-bottom: 8px;
        border-bottom: 3px solid rgba(214, 163, 38, 0.55);
      }

      .help-block {
        color: var(--reeblra-muted);
        line-height: 1.45;
      }

      label.control-label,
      .shiny-input-container > label {
        color: var(--reeblra-green-dark);
        font-weight: 700;
      }

      .form-control,
      .selectize-input,
      textarea.form-control {
        border-radius: 12px !important;
        border: 1px solid rgba(15, 107, 61, 0.24) !important;
        box-shadow: none !important;
        background-color: #ffffff;
      }

      .form-control:focus,
      .selectize-input.focus,
      textarea.form-control:focus {
        border-color: var(--reeblra-blue) !important;
        box-shadow: 0 0 0 3px rgba(31, 95, 153, 0.14) !important;
      }

      .radio label,
      .checkbox label {
        color: var(--reeblra-ink);
        font-weight: 500;
      }

      .mainPanel,
      .tab-content {
        background: rgba(255, 255, 255, 0.93);
        border: 1px solid var(--reeblra-border);
        border-radius: 22px;
        box-shadow: var(--reeblra-shadow);
      }

      .tab-content {
        padding: 22px;
      }

      .nav-tabs {
        border-bottom: 0;
        margin-bottom: 10px;
      }

      .nav-tabs > li > a {
        margin-right: 7px;
        margin-bottom: 8px;
        color: var(--reeblra-maroon);
        background: rgba(255, 255, 255, 0.82);
        border: 1px solid rgba(111, 23, 38, 0.14);
        border-radius: 999px;
        font-weight: 700;
        transition: all 0.16s ease-in-out;
      }

      .nav-tabs > li > a:hover,
      .nav-tabs > li > a:focus {
        color: var(--reeblra-green-dark);
        background: var(--reeblra-gold-soft);
        border-color: rgba(214, 163, 38, 0.60);
      }

      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover,
      .nav-tabs > li.active > a:focus {
        color: #ffffff;
        background: linear-gradient(135deg, var(--reeblra-maroon), var(--reeblra-green));
        border: 1px solid transparent;
        box-shadow: 0 8px 18px rgba(111, 23, 38, 0.22);
      }

      .tab-pane h3 {
        margin-top: 4px;
        color: var(--reeblra-maroon);
        font-weight: 850;
        padding-bottom: 10px;
        border-bottom: 4px solid rgba(214, 163, 38, 0.58);
      }

      .tab-pane h4 {
        color: var(--reeblra-green-dark);
        font-weight: 800;
        margin-top: 22px;
      }

      .tab-pane p,
      .tab-pane li {
        line-height: 1.55;
      }

      .alert-info {
        color: #143d63;
        background: linear-gradient(90deg, var(--reeblra-blue-soft), #ffffff);
        border: 1px solid rgba(31, 95, 153, 0.22);
        border-left: 6px solid var(--reeblra-blue);
        border-radius: 16px;
      }

      .btn,
      .btn-default {
        border-radius: 999px;
        font-weight: 800;
        border-width: 0;
        box-shadow: 0 8px 16px rgba(36, 48, 42, 0.12);
      }

      .btn-primary,
      .action-button.btn-primary {
        background: linear-gradient(135deg, var(--reeblra-maroon), var(--reeblra-maroon-dark));
        color: #ffffff;
      }

      .btn-success,
      .action-button.btn-success {
        background: linear-gradient(135deg, var(--reeblra-green), var(--reeblra-green-dark));
        color: #ffffff;
      }

      .btn-default,
      .btn-info,
      .download-button {
        color: #ffffff !important;
        background: linear-gradient(135deg, var(--reeblra-blue), var(--reeblra-green));
      }

      .btn:hover,
      .btn:focus {
        filter: brightness(1.04);
        transform: translateY(-1px);
      }

      table.dataTable,
      .dataTables_wrapper,
      .table {
        border-radius: 14px;
        overflow: hidden;
      }

      table.dataTable thead th,
      .table > thead > tr > th {
        color: #ffffff;
        background: var(--reeblra-maroon);
        border-bottom: 3px solid var(--reeblra-gold) !important;
      }

      .dataTables_wrapper .dataTables_filter input,
      .dataTables_wrapper .dataTables_length select {
        border-radius: 10px;
        border: 1px solid rgba(15, 107, 61, 0.24);
      }

      pre,
      .shiny-text-output {
        color: #173123;
        background: linear-gradient(180deg, #ffffff, #fffaf0);
        border: 1px solid rgba(214, 163, 38, 0.34);
        border-left: 5px solid var(--reeblra-green);
        border-radius: 15px;
        box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.7);
      }

      hr {
        border-top: 1px solid rgba(111, 23, 38, 0.16);
      }

      .plotly,
      .shiny-plot-output,
      .html-widget {
        border-radius: 18px;
      }

      @media (max-width: 992px) {
        .container-fluid {
          padding: 12px;
        }
        .container-fluid > h2:first-child,
        .app-header,
        .well,
        .tab-content {
          border-radius: 18px;
        }
      }
    "))
  ),

  titlePanel("REEBLRA: Interpretable Machine Learning-based Symbolic and Closed-form Polynomial Surrogate Equations for Outbreak Basic Reproduction Number R0 Estimation"),
  
  tags$div(
    class = "app-header",
    tags$h5(style = "color:maroon", "by Jomar F. Rabajante*, Ariel L. Babierra, Norvin P. Bansilan, Maria Czarina T. Lagura-Arreza, Destiny SM. Lutero, Allen L. Nazareno, and Emerson R. Rico"),
    tags$h5(style = "color:maroon", "Institute of Mathematical Sciences, UPLB"),
    tags$h6(style = "color:maroon", "*Lead and corresponding author: jfrabajante@up.edu.ph"),
    tags$h4(style = "color:green; font-weight:600;", "R0 Equation Estimation using Bayesian Logistic Regression Algorithm (REEBLRA, pronounced “ree-bra”)"),
    tags$h4("Deriving a closed-form symbolic basic reproduction number (R0) equation in infectious disease models often involves finding a nonlinear model-specific formula consisting of parameters with arbitrary/unknown values. This process can be complex and time-consuming, especially for high-dimensional or structurally complicated models. REEBLRA is an alternative and interpretable method for finding closed-form surrogate R0 symbolic equation that leverages numerical simulations and machine learning, specifically logistic regression. This approach not only simplifies the calculation of R0 equation (using logit-based polynomial) but also integrates parameter sensitivity analysis into a unified framework, enhancing the efficiency of model analysis."),
    tags$h5(
      style = "color:green;",
      tags$a(
        href = "https://tinyurl.com/REEBLRAManual",
        target = "_blank",
        rel = "noopener noreferrer",
        style = "color; text-decoration;",
        "User Manual"
      )
    ),
    tags$h5("Important note: REEBLRA does not replace the mechanistic, analytic R0 equation derived from the next-generation matrix or other established mathematical methods. It generates an interpretable, simulation-based symbolic polynomial surrogate for R0 estimation, screening, and exploratory model analysis.")
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Workflow"),
      helpText("Step 1: choose data source. Step 2: choose the R0 definition/outcome. Step 3: fit Bayesian logistic regression. Step 4: inspect the learned threshold equation and compare it with known R0 (if available)."),
      helpText("Tip: if Bayesian fitting fails because all runs are outbreak or non-outbreak, widen parameter ranges or change appropriate settings."),
      radioButtons(
        "data_mode", "1) Data source",
        choices = c("Simulate using built-in ODE model" = "ode",
                    "Custom user-defined ODE model" = "custom_ode",
                    "Upload simulated R0estimate + parameter table" = "upload"),
        selected = "ode"
      ),
      conditionalPanel(
        condition = "input.data_mode == 'ode'",
        selectInput(
          "model_name", "ODE model to simulate",
          choices = c("Closed SIR", "Open SIR", "Open SEIR", "SEIR with vaccination",
                      "Age-structured SIR", "Vector-borne")
        ),
        h5("Parameter ranges for LHS sampling"),
        helpText("You can directly control the minimum and maximum value of every ODE parameter before running simulations."),
        uiOutput("parameter_range_controls"),
        numericInput("n_samples", "Number of LHS simulations", value = 400, min = 50, step = 50),
        textOutput("lhs_recommendation"),
        numericInput("tmax", "Simulation time horizon", value = 365, min = 10, step = 10),
        numericInput("dt", "Time step", value = 1, min = 0.1, step = 0.1),
        numericInput("I0", "Initial infectious proportion", value = 0.01, min = 1e-6, max = 0.2, step = 0.001),
        numericInput("susceptible_threshold", "Early-growth susceptible threshold", value = 0.95, min = 0.50, max = 1.00, step = 0.01),
        helpText("REEBLRA first fits r while the major susceptible population remains at or above this fraction of its initial value after index cases. Default is 0.95. If fewer than 5 positive infectious points are available, the app automatically relaxes the threshold stepwise and warns the user; if there are still too few points in the simulation window, r is not fitted."),
        numericInput("growth_epsilon", "Case-increase tolerance", value = 1e-5, min = 0, step = 1e-5)
      ),
      conditionalPanel(
        condition = "input.data_mode == 'custom_ode'",
        h5("Custom ODE model"),
        helpText("Enter up to 30 state equations. Use R-style syntax: ^ for powers, * for multiplication. Accepted left sides include dx/dt, dI/dt, dx, dI, or I."),
        textAreaInput(
          "custom_equations", "State equations",
          value = "dS/dt = -beta*S*I/N
dI/dt = beta*S*I/N - gamma*I
dR/dt = gamma*I",
          rows = 7
        ),
        textAreaInput(
          "custom_auxiliary_equations", "Equation (not ODE) / auxiliary outputs, optional",
          value = "I_total = I",
          rows = 3
        ),
        helpText("Use non-ODE equations for derived quantities, totals, forcing terms, or intervention indicators. Example: I_total = I1 + I2 can be selected as the early-growth target when one species has multiple infectious classes."),
        helpText("N is automatically set to the sum of state variables unless N is declared in the equations through states/parameters/non-ODE equations. The symbols t and time are also available."),
        textAreaInput(
          "custom_initial_values", "Initial values, one per line",
          value = "S = 0.99
I = 0.01
R = 0",
          rows = 4
        ),
        uiOutput("custom_nonnegative_state_controls"),
        helpText("Checked state variables are forced to remain nonnegative during custom ODE simulation. Uncheck a state only if negative values are meaningful for that variable."),
        textAreaInput(
          "custom_parameter_ranges", "Parameter ranges, one per line: parameter = min, max",
          value = "beta = 0.01, 1
gamma = 0.01, 1",
          rows = 5
        ),
        textInput("custom_target_state", "Target infectious/output state for early growth", value = "I"),
        numericInput("custom_n_samples", "Number of LHS simulations", value = 400, min = 20, step = 50),
        textOutput("custom_lhs_recommendation"),
        numericInput("custom_tmax", "Simulation time horizon", value = 365, min = 1, step = 10),
        numericInput("custom_dt", "Time step", value = 1, min = 0.001, step = 0.1),
        numericInput("custom_growth_epsilon", "Case-increase tolerance", value = 1e-5, min = 0, step = 1e-5),
        numericInput("custom_growth_window_first_n", "Early-growth points to use", value = 8, min = 5, step = 1),
        numericInput("custom_min_positive_value", "Minimum positive target value for log-growth fit", value = 1e-12, min = 0, step = 1e-12),
        h5("Growth-rate-to-R0estimate conversion"),
        selectInput(
          "custom_conversion_method", "Conversion method",
          choices = c(
            "SIR-like: R0 = 1 + r/gamma_removal",
            "SEIR-like with latency: R0 = ((r+sigma)(r+gamma))/(sigma*gamma)",
            "Generation-time approximation: R0 = exp(r*Tgen)",
            "Unit-time Euler multiplier: R0estimate = 1 + r (one time unit)"
          ),
          selected = "SIR-like: R0 = 1 + r/gamma_removal"
        ),
        numericInput("custom_gamma_removal", "gamma_removal", value = 0.1, min = 0, step = 0.01),
        numericInput("custom_sigma_latency", "sigma_latency", value = 0.2, min = 0, step = 0.01),
        numericInput("custom_generation_time", "Generation time Tgen", value = 5, min = 0, step = 0.1),
        selectInput(
          "custom_nonpositive_rule", "Rule if estimated r <= 0",
          choices = c(
            "Set R0estimate to 0" = "zero",
            "Set R0estimate to 1" = "one",
            "Keep formula value" = "keep_formula",
            "Set R0estimate to missing" = "missing"
          ),
          selected = "zero"
        ),
        helpText("Custom ODE mode generates R0estimate from early growth of the selected target state or optional non-ODE auxiliary/output equation. If one species has multiple infectious classes, define their total in the non-ODE equation box and use that total for growth fitting. It does not automatically derive an analytic next-generation-matrix R0.")
      ),
      conditionalPanel(
        condition = "input.data_mode == 'upload'",
        fileInput("upload_r0_file", "Upload CSV table", accept = c(".csv")),
        helpText("CSV must contain one column named R0estimate and parameter columns such as param1, param2, ... . The user should estimate R0estimate numerically before upload (e.g., using early growth rate fitting, or numerical spectral radius NGM). Recommended: fit early growth rate r, check r > 0, report R-squared and 95% CI for r, then convert r to R0estimate using removal/latency assumptions. By REEBLRA convention, numerical growth-based R0estimate values below 1 may be set to 0 to flag non-growing/subthreshold runs. REEBLRA then learns a surrogate or alternative symbolic closed-form threshold formula for R0 using logistic regression. Analytic R0/Re is not available for uploaded datasets."),
        uiOutput("upload_column_controls"),
        helpText("Suggested maximum: about 50 parameter columns. More than 50 can work, but Bayesian fitting and interpretation may become slow or unstable; use screening/regularization before this app if many parameters are available.")
      ),

      hr(),
      conditionalPanel(
        condition = "input.data_mode == 'ode'",
        selectInput(
          "threshold_source", "2a) Binary threshold outcome",
          choices = c(
            "Cases increase" = "cases_increase",
            "Analytic R0 > cutoff" = "analytic",
            "Numerical growth-rate R0estimate > cutoff" = "numeric"
          )
        ),
        selectInput(
          "multiclass_source", "2b) Multiclass R0 bracket source",
          choices = c("Analytic R0/Re" = "R0_bracket_analytic",
                      "Numerical growth-rate R0estimate" = "R0_bracket_numeric")
        )
      ),
      conditionalPanel(
        condition = "input.data_mode == 'custom_ode'",
        helpText("For custom ODEs, the app uses the simulated numerical R0estimate from early growth. Binary uses R0estimate > cutoff and multiclass uses R0estimate brackets. Analytic R0/Re is not available automatically.")
      ),
      conditionalPanel(
        condition = "input.data_mode == 'upload'",
        helpText("For uploaded tables, the app uses R0estimate directly. Thus, there is no 2a/2b choice: binary uses R0estimate > cutoff and multiclass uses R0estimate brackets.")
      ),
      numericInput("r0_cutoff", "Binary R0 cutoff (Not Applicable to the “Cases Increase” Option)", value = 1, min = 0.25, step = 0.25),
      selectInput(
        "multiclass_step",
        "Multiclass R0 bracket step size",
        choices = c("0.25" = 0.25, "0.5" = 0.5, "1" = 1, "2" = 2, "5" = 5),
        selected = 0.25
      ),
      helpText("This controls the successive binary thresholds and bracket labels. Examples: 0.25 gives <1, 1-1.25, ...; 0.5 gives <1, 1-1.5, ...; 1 gives <1, 1-2, ... ."),

      hr(),
      numericInput("chains", "Bayesian chains", value = 2, min = 1, max = 4, step = 1),
      numericInput("iter", "Iterations per chain", value = 1000, min = 500, step = 500),
      helpText("Use MCMC sampling for final runs. Use variational Bayes only for fast exploratory checking because it is approximate."),
      numericInput("seed", "Random seed", value = 123, min = 1, step = 1),
      numericInput("k_folds", "K-fold validation: number of folds", value = 2, min = 1, max = 20, step = 1),
      helpText("Use k = 1 to report only apparent/in-sample metrics. Use k >= 2 for train/test validation by k-fold cross-validation. Binary Bayesian k-fold fitting can be slow because the model is refitted in each fold."),
      selectInput("stan_algorithm", "Binary Bayes algorithm",
                  choices = c("MCMC sampling" = "sampling",
                              "Variational Bayes / fast exploratory" = "meanfield"),
                  selected = "meanfield"),
      hr(),
      actionButton("run_sim", "Run simulations / load uploaded table", class = "btn-primary"),
      actionButton("fit_binary", "Fit binary Bayesian logistic R0", class = "btn-success"),
      actionButton("fit_multi", "Fit multiclass R0 bracket model", class = "btn-warning")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Overview",
          h3("Purpose and goal"),
          p("This app estimates an R0-like equation from infectious disease model (e.g., simulation results) or uploaded dataset using Bayesian logistic regression. The goal is to learn a simulation-based epidemic-threshold equation from simulated parameter combinations, then compare it with the known analytic R0/Re formula when available."),
          h4("Main workflow"),
          tags$ol(
            tags$li(strong("Choose a data source (e.g., an ODE model.)"), " The app currently supports closed SIR, open SIR, open SEIR, SEIR with vaccination, age-structured SIR, vector-borne transmission, custom user-defined ODEs, and uploaded data."),
            tags$li(strong("Simulate many parameter combinations."), " Latin Hypercube Sampling is used to generate parameter sets within the editable ranges."),
            tags$li(strong("Define the R0 target."), " You may classify case increase, analytic R0 > cutoff, or numerical growth-rate R0estimate > cutoff."),
            tags$li(strong("Fit Bayesian logistic regression."), " The fitted equation gives Pr(threshold is true) = logit^{-1}(a + b1*x1 + ... + bp*xp)."),
            tags$li(strong("Compare and interpret."), " The app reports posterior coefficients, the logistic equation, performance metrics, ROC/AUC, and analytic-versus-numeric R0 agreement.")
          ),
          h4("Two R0 concepts in this app"),
          tags$ul(
            tags$li(strong("Threshold R0:"), " asks whether cases increase or whether R0 exceeds a cutoff. This is modeled as a binary classification problem."),
            tags$li(strong("Growth-rate R0estimate:"), " estimates the early exponential growth rate r from log(I_major(t)) = a + r t while the major susceptible population is still approximately undepleted. The app reports r, 95% CI for r, R-squared, fitting-window size, and warnings. R0estimate is computed from r using the removal rate and, for SEIR-like built-in models, the latency/progression rate. If too few points are available at the selected susceptible threshold, the app relaxes the threshold stepwise to obtain at least 5 positive infectious points when possible. By convention, growth-based R0estimate values below 1 are set to 0 to flag non-growing/subthreshold numerical runs.")
          ),
          h4("Analytic formula for selected model"),
          verbatimTextOutput("formula_text"),
          h4("How to read the binary Bayesian logistic equation"),
          tags$ul(
            tags$li("If the selected cutoff is 1 and p > 0.5, the model classifies the parameter set as R0 > 1."),
            tags$li("Positive coefficients increase the posterior probability of crossing the threshold; negative coefficients decrease it."),
            tags$li("Predictors are standardized internally, so coefficients are interpreted on z-scored parameter scales."),
            tags$li("The equation is an in silico approximation, not a substitute for mechanistic derivation. Its usefulness should be judged by AUC, accuracy, sensitivity, specificity, and agreement with analytic R0 (if available).")
          ),
          h4("How to read the multiclass ordered-bracket equations"),
          tags$ul(
            tags$li("The multiclass output does not ask only whether R0 exceeds one cutoff. It estimates the likely R0 bracket, such as <1, 1-1.25, 1.25-1.5, and so on, depending on the selected step size."),
            tags$li("The app fits successive binary equations: R0 > 1, R0 > 1 + step, R0 > 1 + 2*step, and so on up to R0 > 10."),
            tags$li("These cumulative probabilities are converted into bracket probabilities. The predicted bracket is the bracket with the largest predicted probability."),
            tags$li("Multiclass performance should be judged using overall accuracy, weighted averages, class distribution, one-vs-rest AUC, and the confusion matrix. Sparse brackets can have unstable sensitivity and precision even when the overall accuracy is high.")
          )
        ),

        tabPanel(
          "ODE equations",
          h3("ODE equations used in the selected model"),
          p("This tab shows the exact deterministic ODE system used by the simulator. These equations define the mechanism that generates the training data for the Bayesian logistic R0 equation."),
          verbatimTextOutput("ode_text"),
          h4("Analytic R0/Re used for comparison"),
          verbatimTextOutput("formula_text2"),
          h4("Parameter definitions"),
          DTOutput("param_def_table")
        ),

        tabPanel(
          "Parameter ranges",
          p("Edit the ranges before simulation. The app samples uniformly within these bounds using Latin Hypercube Sampling."),
          DTOutput("param_table"),
          verbatimTextOutput("param_note")
        ),

        tabPanel(
          "Simulation results",
          fluidRow(
            column(6, h4("Simulation summary"), tableOutput("sim_summary")),
            column(6, h4("R0 source comparison"), tableOutput("r0_compare"))
          ),
          DTOutput("sim_table"),
          br(),
          downloadButton("download_sim", "Download simulation CSV")
        ),

        tabPanel(
          "Collinearity analysis",
          h3("Collinearity analysis"),
          p("This tab checks whether sampled or uploaded parameters are strongly correlated with each other. VIF is useful to report when the logistic equation is interpreted through individual coefficients, because collinear parameters can make coefficient signs and magnitudes unstable even when predictive performance is good."),
          p("Interpretation guide: VIF near 1 means little collinearity; VIF around 2.5 to 5 suggests moderate collinearity; VIF above 5 suggests high collinearity; VIF above 10 suggests serious collinearity. High VIF does not automatically invalidate prediction, but it weakens parameter-by-parameter interpretation."),
          h4("Variance inflation factor (VIF)"),
          DTOutput("vif_table"),
          h4("Parameter correlation matrix"),
          plotlyOutput("corr_matrix_plot", height = "620px"),
          h4("Correlation matrix table"),
          DTOutput("corr_matrix_table"),
          br(),
          downloadButton("download_vif_csv", "Download VIF CSV"),
          downloadButton("download_vif_png", "Download VIF PNG"),
          downloadButton("download_corr_matrix_csv", "Download correlation matrix CSV"),
          downloadButton("download_corr_matrix_png", "Download correlation matrix PNG")
        ),

        tabPanel(
          "Binary Bayesian R0 threshold (Main output)",
          h3("Binary Bayesian logistic regression"),
          p("This section estimates a threshold equation. The response is 1 if the selected condition is true, for example analytic R0 > cutoff or cases increase, and 0 otherwise."),
          p("The fitted Bayesian equation is useful as a simulation-based R0-threshold approximation. For cutoff = 1, it mimics the classical rule: outbreak if R0 > 1."),
          div(class = "alert alert-info", "After clicking Fit binary Bayesian logistic R0, please wait until the fitting notification disappears. Variational Bayes is faster; MCMC and k-fold validation may take longer because the model may be refitted several times."),
          h4("How to read the Bayesian logistic equation"),
          tags$ul(
            tags$li("If the selected cutoff is 1 and p > 0.5, the model classifies the parameter set as R0 > 1."),
            tags$li("Positive coefficients increase the posterior probability of crossing the threshold; negative coefficients decrease it."),
            tags$li("Predictors are standardized internally, so coefficients are interpreted on z-scored parameter scales."),
            tags$li("The equation is an in silico approximation, not a substitute for mechanistic derivation. Its usefulness should be judged by AUC, accuracy, sensitivity, specificity, and agreement with analytic R0 (if available).")
          ),
          h4("Estimated equation"),
          verbatimTextOutput("binary_equation"),
          h4("Performance metrics"),
          p("When k >= 2, the table includes k-fold cross-validation metrics from held-out folds, plus apparent/in-sample metrics from the final model fitted to all data."),
          tableOutput("binary_metrics"),
          h4("MCMC convergence diagnostics"),
          p("Shown for MCMC sampling. Variational Bayes is approximate and does not have MCMC convergence diagnostics."),
          DTOutput("mcmc_diagnostics"),
          h4("Interactive trial on ODE/original parameter scale"),
          p("Enter parameter values on the original ODE scale. The app evaluates the fitted Bayesian logistic R0-threshold equation and, for built-in ODE models, compares it with analytic R0/Re."),
          uiOutput("binary_trial_inputs"),
          verbatimTextOutput("binary_trial_output"),
          h4("Posterior coefficients"),
          DTOutput("binary_coef"),
          h4("ROC curve"),
          plotOutput("binary_roc", height = "360px"),
          h4("Predicted probability surface / scatter"),
          uiOutput("binary_axis_controls"),
          plotlyOutput("binary_plotly", height = "520px"),
          br(),
          downloadButton("download_binary_metrics", "Download binary metrics CSV"),
          downloadButton("download_binary_predictions", "Download binary predictions CSV"),
          hr(),
          h3("Ordinary (non-Bayesian) binary logistic regression"),
          p("For comparison, this model uses the same binary outcome, complete rows, standardized predictors, and full-data formula as the Bayesian binary model. The equations below show the linear predictor eta and the corresponding probability p."),
          verbatimTextOutput("binary_ordinary_equation")
        ),

        tabPanel(
          "Multiclass R0 bracket (Main output)",
          h3("Multiclass R0 bracket estimation"),
          p("This section classifies each parameter set into R0 brackets: <1, 1-(1+step), ..., >10. It helps estimate not only whether R0 exceeds 1, but also the likely R0 range."),
          p("This app uses successive binary logistic equations for ordered R0 brackets: R0 > 1, then the next selected step-size threshold, and so on until R0 > 10. The step size can be 0.25, 0.5, 1, 2, or 5. This is more appropriate for ordered R0 brackets than one-time nominal multiclass classification."),
          div(class = "alert alert-info", "After clicking Fit multiclass R0 bracket model, please wait until the fitting notification disappears. The result can appear with a short delay, especially when many thresholds, many parameters, or k-fold validation are used."),
          h4("How to read the multiclass ordered-bracket output"),
          tags$ul(
            tags$li("The model estimates ordered events of the form R0 > c for each selected threshold c."),
            tags$li("The estimated probabilities Pr(R0 > c) are converted into probabilities for each R0 bracket."),
            tags$li("The predicted R0 bracket is the bracket with the highest converted probability."),
            tags$li("Positive coefficients in a threshold equation increase the probability of being above that specific threshold; negative coefficients decrease it."),
            tags$li("Judge multiclass usefulness using overall accuracy, weighted averages, class distribution, one-vs-rest AUC, and the confusion matrix. Sparse brackets may have poor sensitivity or precision even when the overall accuracy is acceptable. In that case, increase the sample size or use a larger bracket step.")
          ),
          h4("Model used"),
          verbatimTextOutput("multi_method"),
          h4("Multiclass equation"),
          verbatimTextOutput("multi_equation"),
          h4("Coefficient matrix across R0 levels"),
          p("Columns are parameters and rows are successive R0 threshold levels. Cell values are standardized logistic coefficients. The color gradient is normalized separately within each row by dividing by that row's largest absolute coefficient, so it shows how relative parameter importance and direction change across R0 levels."),
          plotlyOutput("multi_coef_heatmap", height = "760px"),
          h4("Multiclass coefficients"),
          DTOutput("multi_coef"),
          h4("Interactive trial on ODE/original parameter scale"),
          p("Enter parameter values on the original scale. The app computes successive-threshold probabilities for all R0 brackets and, for built-in ODE models, compares the predicted bracket with analytic R0/Re."),
          uiOutput("multi_trial_inputs"),
          verbatimTextOutput("multi_trial_output"),
          h4("Multiclass accuracy and one-vs-rest metrics"),
          p("When k >= 2, the table includes k-fold cross-validation metrics from held-out folds, plus apparent/in-sample metrics from the final model fitted to all data."),
          tableOutput("multi_metrics"),
          h4("Class distribution plot"),
          p("Use this plot to check whether the R0 brackets are too sparse or imbalanced. Very small observed counts usually lead to unstable sensitivity and precision."),
          plotlyOutput("multi_class_distribution_plot", height = "420px"),
          h4("One-vs-rest AUC plot per class"),
          p("Each bar shows how well the model discriminates one R0 bracket against all other brackets. AUC near 0.5 is weak; higher values indicate better discrimination."),
          plotlyOutput("multi_auc_plot", height = "420px"),
          h4("K-fold one-vs-rest AUC plot per class"),
          p("This plot is shown when k >= 2 and held-out fold predictions are available."),
          plotlyOutput("multi_cv_auc_plot", height = "420px"),
          h4("Predicted vs observed brackets"),
          tableOutput("multi_confusion"),
          h4("Predicted bracket plot"),
          uiOutput("multi_axis_controls"),
          plotlyOutput("multi_plotly", height = "520px"),
          br(),
          downloadButton("download_multi_predictions", "Download multiclass predictions CSV")
        ),

        tabPanel(
          "Parameter sensitivity analysis",
          h3("Normalized Bayesian sensitivity coefficients"),
          p("Aim: this tab summarizes which parameters most strongly influence the fitted probability of crossing the selected R0 threshold in the binary Bayesian logistic model."),
          p("The app first fits the binary Bayesian logistic equation using standardized predictors. For each posterior draw of each standardized coefficient b, it computes a normalized Bayesian sensitivity coefficient s = tanh(b). This transformation maps any real-valued coefficient to the interval (-1, 1), making effects easier to compare across parameters."),
          p("Interpretation: positive values mean the parameter tends to increase the probability of crossing the selected threshold; negative values mean the parameter tends to decrease that probability. Values near 0 indicate weak model-based influence, while values near -1 or +1 indicate stronger negative or positive influence. The table and plot report the posterior median and 95% credible interval of the transformed coefficient."),
          plotlyOutput("sensitivity_plot", height = "520px"),
          h4("Sensitivity table"),
          DTOutput("sensitivity_table")
        ),



        tabPanel(
          "CSV generator",
          h3("CSV generator for uploaded time-series data"),
          p("This tab creates a REEBLRA-ready upload table from time-series trajectories and matching parameter combinations. The output format is: R0estimate, param1, param2, ..., paramp."),
          div(class = "alert alert-info",
              strong("Required format. "),
              "The first uploaded file must have time in column 1 and one epidemic trajectory per remaining column. The second uploaded file must have one row per trajectory column. Thus, sim1 corresponds to row 1, sim2 to row 2, and so on."),
          sidebarLayout(
            sidebarPanel(
              width = 4,
              h4("1. Upload input files"),
              fileInput("csvgen_ts_file", "Time-series file (.csv, .txt, .tsv)", accept = c(".csv", ".txt", ".tsv")),
              fileInput("csvgen_par_file", "Parameter-combination file (.csv, .txt, .tsv)", accept = c(".csv", ".txt", ".tsv")),
              hr(),
              h4("2. Early-growth window"),
              selectInput(
                "csvgen_window_mode",
                "How should the early-growth window be selected?",
                choices = c(
                  "Automatic: choose early log-linear window" = "automatic",
                  "First N positive points" = "first_n",
                  "Use time range" = "time_range",
                  "Use points below fraction of peak" = "below_fraction_of_peak",
                  "Use all positive points" = "all_positive"
                ),
                selected = "automatic"
              ),
              numericInput("csvgen_first_n", "Automatic maximum points / first-N fallback", value = 8, min = 2, step = 1),
              numericInput("csvgen_time_min", "Start time if using time range", value = 0, step = 1),
              numericInput("csvgen_time_max", "End time if using time range", value = 14, step = 1),
              numericInput("csvgen_peak_fraction", "Fraction of peak if using below-peak rule", value = 0.30, min = 0.01, max = 1, step = 0.05),
              numericInput("csvgen_min_points", "Minimum points required for growth fit", value = 5, min = 3, step = 1),
              numericInput("csvgen_min_positive_value", "Minimum positive value before log transform", value = 1e-12, min = 0, step = 1e-12),
              hr(),
              h4("3. Convert growth rate r to R0estimate"),
              selectInput(
                "csvgen_conversion_method",
                "Conversion method",
                choices = c(
                  "SIR-like: R0 = 1 + r/gamma_removal",
                  "SEIR-like with latency: R0 = ((r+sigma)(r+gamma))/(sigma*gamma)",
                  "Generation-time approximation: R0 = exp(r*Tgen)",
                  "Unit-time Euler multiplier: R0estimate = 1 + r (one time unit)"
                ),
                selected = "SIR-like: R0 = 1 + r/gamma_removal"
              ),
              selectInput(
                "csvgen_rate_source",
                "Use fixed conversion values or columns from parameter file?",
                choices = c("Use fixed values below" = "fixed", "Use parameter-file columns" = "columns"),
                selected = "fixed"
              ),
              conditionalPanel(
                condition = "input.csvgen_rate_source == 'fixed'",
                numericInput("csvgen_gamma_removal", "gamma_removal: recovery/removal rate", value = 0.20, min = 0, step = 0.01),
                numericInput("csvgen_sigma_latency", "sigma_latency: exposed-to-infectious progression rate", value = 0.25, min = 0, step = 0.01),
                numericInput("csvgen_generation_time", "Mean generation time Tgen", value = 5, min = 0, step = 0.5)
              ),
              conditionalPanel(
                condition = "input.csvgen_rate_source == 'columns'",
                uiOutput("csvgen_rate_column_controls"),
                helpText("Use this when gamma, sigma, or generation time differs by simulation row. For SIR-like conversion, choose a gamma/removal column. For SEIR-like conversion, choose gamma/removal and sigma/latency columns.")
              ),
              selectInput(
                "csvgen_nonpositive_rule",
                "What to do when estimated r <= 0?",
                choices = c(
                  "Set R0estimate to 0" = "zero",
                  "Set R0estimate to 1" = "one",
                  "Keep formula value" = "keep_formula",
                  "Set R0estimate to missing" = "missing"
                ),
                selected = "zero"
              ),
              hr(),
              actionButton("csvgen_run", "Generate REEBLRA CSV", class = "btn-primary"),
              br(), br(),
              actionButton("csvgen_use_in_app", "Use generated CSV in REEBLRA app", class = "btn-success"),
              br(), br(),
              downloadButton("csvgen_download_csv", "Download output CSV"),
              downloadButton("csvgen_download_diagnostics", "Download diagnostics CSV")
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                tabPanel(
                  "Instructions",
                  br(),
                  h4("Purpose"),
                  p("This tab estimates an early exponential growth rate for each uploaded trajectory and converts it to R0estimate. The generated table can then be used directly in the binary and multiclass REEBLRA tabs."),
                  h4("Input 1: time-series file"),
                  tags$ul(
                    tags$li("Column 1 must be time, such as t or time."),
                    tags$li("Columns 2 to n must be epidemic trajectories: sim1, sim2, ..., simN."),
                    tags$li("Values should be non-negative. Zero and negative values are excluded before fitting log(Y)."),
                    tags$li("Each trajectory should represent the same epidemic signal, preferably infectious prevalence, infected hosts, or incidence.")
                  ),
                  tags$pre("t,sim1,sim2,sim3\n0,1,1,1\n1,1.20,0.95,1.50\n2,1.44,0.90,2.25"),
                  h4("Input 2: parameter-combination file"),
                  tags$ul(
                    tags$li("Each row corresponds to one trajectory column in the time-series file."),
                    tags$li("Row 1 corresponds to sim1, row 2 to sim2, and so on."),
                    tags$li("Columns should be parameter names such as beta, gamma, sigma, mu, c11, c12, or other model parameters."),
                    tags$li("Do not include R0estimate. This tab generates it.")
                  ),
                  tags$pre("beta,gamma,mu\n0.30,0.20,0.01\n0.15,0.20,0.01\n0.55,0.20,0.01"),
                  h4("Method"),
                  p("For each trajectory, the app fits log(Y(t)) = a + r t over an early-growth window. The default Automatic option tests early windows up to the first peak and selects a short window with high log-linear R-squared and positive growth when possible."),
                  h4("Latency option"),
                  p("Latency is needed only when the model has an exposed/latent stage and the SEIR-like conversion is desired. Without latency, use the SIR-like conversion. If gamma and sigma vary by simulation, choose 'Use parameter-file columns' so conversion uses row-specific rates.")
                ),
                tabPanel("Input preview", br(), h4("Time-series file"), DTOutput("csvgen_ts_preview"), br(), h4("Parameter file"), DTOutput("csvgen_par_preview")),
                tabPanel("Generated output CSV", br(), DTOutput("csvgen_output_table")),
                tabPanel("Diagnostics", br(), DTOutput("csvgen_diagnostics_table")),
                tabPanel("Warnings", br(), verbatimTextOutput("csvgen_warnings_text"))
              )
            )
          )
        ),

        tabPanel(
          "Downloads",
          h3("Download CSV, PNG, and LaTeX outputs"),
          p("Download individual outputs below, or download a ZIP file containing all currently available CSV, PNG, and LaTeX outputs. Outputs that require a fitted model will be available after fitting that model."),
          h4("All outputs"),
          downloadButton("download_all_outputs_zip", "Download all available outputs as ZIP"),
          hr(),
          h4("CSV outputs"),
          fluidRow(
            column(4,
                   downloadButton("download_sim2", "Simulation / uploaded data CSV"),
                   br(), br(),
                   downloadButton("download_param_ranges_csv", "Parameter ranges / summary CSV"),
                   br(), br(),
                   downloadButton("download_param_definitions_csv", "Parameter definitions CSV"),
                   br(), br(),
                   downloadButton("download_sim_summary_csv", "Simulation summary CSV"),
                   br(), br(),
                   downloadButton("download_r0_compare_csv", "R0 comparison summary CSV"),
                   br(), br(),
                   downloadButton("download_vif_csv2", "Collinearity VIF CSV"),
                   br(), br(),
                   downloadButton("download_corr_matrix_csv2", "Parameter correlation matrix CSV")
            ),
            column(4,
                   downloadButton("download_binary_metrics2", "Binary metrics CSV"),
                   br(), br(),
                   downloadButton("download_binary_predictions2", "Binary predictions CSV"),
                   br(), br(),
                   downloadButton("download_binary_cv_predictions_csv", "Binary k-fold predictions CSV"),
                   br(), br(),
                   downloadButton("download_binary_roc_csv", "Binary ROC data CSV"),
                   br(), br(),
                   downloadButton("download_binary_coefficients_csv", "Binary posterior coefficients CSV"),
                   br(), br(),
                   downloadButton("download_mcmc_diagnostics_csv", "MCMC diagnostics CSV"),
                   br(), br(),
                   downloadButton("download_sensitivity_csv", "Parameter sensitivity CSV")
            ),
            column(4,
                   downloadButton("download_multi_predictions2", "Multiclass predictions CSV"),
                   br(), br(),
                   downloadButton("download_multi_cv_predictions_csv", "Multiclass k-fold predictions CSV"),
                   br(), br(),
                   downloadButton("download_multi_metrics_csv", "Multiclass metrics CSV"),
                   br(), br(),
                   downloadButton("download_multi_coefficients_csv", "Multiclass coefficients CSV"),
                   br(), br(),
                   downloadButton("download_multi_confusion_csv", "Multiclass confusion matrix CSV"),
                   br(), br(),
                   downloadButton("download_multi_class_distribution_csv", "Class distribution CSV"),
                   br(), br(),
                   downloadButton("download_multi_auc_csv", "One-vs-rest AUC CSV"),
                   br(), br(),
                   downloadButton("download_multi_cv_auc_csv", "K-fold one-vs-rest AUC CSV")
            )
          ),
          hr(),
          h4("LaTeX equation files"),
          p("These files contain publication-ready LaTeX equations with discussion and interpretation for the fitted binary and multiclass REEBLRA models."),
          fluidRow(
            column(4,
                   downloadButton("download_binary_equation_tex", "Binary equation LaTeX (.tex)"),
                   br(), br(),
                   downloadButton("download_multi_equation_tex", "Multiclass equations LaTeX (.tex)")
            )
          ),
          hr(),
          h4("PNG outputs"),
          fluidRow(
            column(4,
                   downloadButton("download_binary_roc_png", "Binary ROC PNG"),
                   br(), br(),
                   downloadButton("download_binary_probability_png", "Binary probability plot PNG")
            ),
            column(4,
                   downloadButton("download_multi_prediction_png", "Multiclass prediction plot PNG"),
                   br(), br(),
                   downloadButton("download_multi_class_distribution_png", "Class distribution plot PNG"),
                   br(), br(),
                   downloadButton("download_multi_auc_png", "One-vs-rest AUC plot PNG"),
                   br(), br(),
                   downloadButton("download_multi_cv_auc_png", "K-fold one-vs-rest AUC plot PNG")
            ),
            column(4,
                   downloadButton("download_sensitivity_png", "Parameter sensitivity plot PNG"),
                   br(), br(),
                   downloadButton("download_vif_png2", "Collinearity VIF PNG"),
                   br(), br(),
                   downloadButton("download_corr_matrix_png2", "Parameter correlation matrix PNG")
            )
          )
        )
      )
    )
  )
)

# =========================================================
# SHINY SERVER
# =========================================================

server <- function(input, output, session) {

  param_ranges <- reactiveVal(model_parameter_table("Closed SIR"))
  sim_data <- reactiveVal(NULL)
  binary_fit <- reactiveVal(NULL)
  multi_fit <- reactiveVal(NULL)

  data_ready_message <- reactive({
    if (identical(input$data_mode, "upload")) {
      "No uploaded R0 table is available yet. Please upload a CSV and click 'Load uploaded R0 table' first."
    } else if (identical(input$data_mode, "custom_ode")) {
      "No custom ODE simulation data is available yet. Please enter equations and click 'Run simulations / load uploaded table' first."
    } else {
      "No ODE simulation data is available yet. Please click 'Run ODE simulation' first."
    }
  })

  binary_ready_message <- reactive({
    "No binary Bayesian model has been fitted yet. Please prepare data first, then click 'Fit binary Bayesian model'."
  })

  multi_ready_message <- reactive({
    "No multiclass R0 bracket model has been fitted yet. Please prepare data first, then click 'Fit multiclass R0 bracket model'."
  })

  note_table <- function(message) {
    data.frame(Note = message, check.names = FALSE)
  }

  note_plotly <- function(message) {
    plotly_empty(type = "scatter", mode = "markers") %>%
      layout(title = message)
  }

  write_note_csv <- function(file, message) {
    write.csv(note_table(message), file, row.names = FALSE)
  }



  # -----------------------------
  # CSV generator tab server
  # -----------------------------
  csvgen_ts_data <- reactive({
    req(input$csvgen_ts_file)
    csvgen_read_uploaded_table(input$csvgen_ts_file)
  })

  csvgen_par_data <- reactive({
    req(input$csvgen_par_file)
    csvgen_read_uploaded_table(input$csvgen_par_file)
  })

  output$csvgen_ts_preview <- renderDT({
    req(csvgen_ts_data())
    DT::datatable(head(csvgen_ts_data(), 20), options = list(scrollX = TRUE, pageLength = 10))
  })

  output$csvgen_par_preview <- renderDT({
    req(csvgen_par_data())
    DT::datatable(head(csvgen_par_data(), 20), options = list(scrollX = TRUE, pageLength = 10))
  })

  output$csvgen_rate_column_controls <- renderUI({
    req(csvgen_par_data())
    pars <- csvgen_par_data()
    num_cols <- names(pars)[sapply(pars, is.numeric)]
    if (length(num_cols) == 0) return(helpText("No numeric parameter columns detected."))
    tagList(
      selectInput("csvgen_gamma_col", "Column for gamma_removal", choices = c("None" = "", num_cols), selected = if ("gamma" %in% num_cols) "gamma" else num_cols[1]),
      selectInput("csvgen_sigma_col", "Column for sigma_latency", choices = c("None" = "", num_cols), selected = if ("sigma" %in% num_cols) "sigma" else ""),
      selectInput("csvgen_tgen_col", "Column for generation time Tgen", choices = c("None" = "", num_cols), selected = "")
    )
  })

  csvgen_generated <- eventReactive(input$csvgen_run, {
    ts <- csvgen_ts_data()
    pars <- csvgen_par_data()

    validate(
      need(ncol(ts) >= 2, "The time-series file must have at least two columns: time and one trajectory."),
      need(nrow(pars) == ncol(ts) - 1,
           paste0("Mismatch: parameter file has ", nrow(pars), " rows, but the time-series file has ",
                  ncol(ts) - 1, " trajectory columns. These must match."))
    )

    time <- csvgen_safe_numeric(ts[[1]])
    validate(need(all(is.finite(time)), "The first column of the time-series file must be numeric time values."))

    traj_names <- names(ts)[-1]
    diag_list <- vector("list", length(traj_names))
    r0_values <- rep(NA_real_, length(traj_names))

    for (i in seq_along(traj_names)) {
      fit_diag <- csvgen_fit_growth_one_series(
        time = time,
        y = ts[[traj_names[i]]],
        min_positive_value = input$csvgen_min_positive_value,
        window_mode = input$csvgen_window_mode,
        first_n = input$csvgen_first_n,
        time_min = input$csvgen_time_min,
        time_max = input$csvgen_time_max,
        peak_fraction = input$csvgen_peak_fraction,
        min_points = input$csvgen_min_points
      )

      if (identical(input$csvgen_rate_source, "columns")) {
        gamma_i <- if (!is.null(input$csvgen_gamma_col) && nzchar(input$csvgen_gamma_col) && input$csvgen_gamma_col %in% names(pars)) csvgen_safe_numeric(pars[[input$csvgen_gamma_col]][i]) else NA_real_
        sigma_i <- if (!is.null(input$csvgen_sigma_col) && nzchar(input$csvgen_sigma_col) && input$csvgen_sigma_col %in% names(pars)) csvgen_safe_numeric(pars[[input$csvgen_sigma_col]][i]) else NA_real_
        tgen_i <- if (!is.null(input$csvgen_tgen_col) && nzchar(input$csvgen_tgen_col) && input$csvgen_tgen_col %in% names(pars)) csvgen_safe_numeric(pars[[input$csvgen_tgen_col]][i]) else NA_real_
      } else {
        gamma_i <- input$csvgen_gamma_removal
        sigma_i <- input$csvgen_sigma_latency
        tgen_i <- input$csvgen_generation_time
      }

      r0_values[i] <- csvgen_convert_r_to_R0(
        r = fit_diag$r,
        conversion_method = input$csvgen_conversion_method,
        gamma_removal = gamma_i,
        sigma_latency = sigma_i,
        generation_time = tgen_i,
        nonpositive_rule = input$csvgen_nonpositive_rule
      )

      diag_list[[i]] <- fit_diag %>%
        dplyr::mutate(
          series = traj_names[i],
          R0estimate = r0_values[i],
          conversion_method = input$csvgen_conversion_method,
          gamma_removal_used = gamma_i,
          sigma_latency_used = sigma_i,
          generation_time_used = tgen_i,
          .before = 1
        )
    }

    diagnostics <- dplyr::bind_rows(diag_list)

    # Keep the downloadable/generated CSV clean and upload-ready:
    # R0estimate, param1, param2, ..., paramp.
    # Internal REEBLRA columns are added only when the user clicks
    # "Use generated CSV in REEBLRA app".
    output_csv <- dplyr::bind_cols(tibble::tibble(R0estimate = r0_values), pars)

    list(output_csv = output_csv, diagnostics = diagnostics)
  })

  output$csvgen_output_table <- renderDT({
    req(csvgen_generated())
    DT::datatable(csvgen_generated()$output_csv, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$csvgen_diagnostics_table <- renderDT({
    req(csvgen_generated())
    DT::datatable(csvgen_generated()$diagnostics, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$csvgen_warnings_text <- renderText({
    req(csvgen_generated())
    diag <- csvgen_generated()$diagnostics
    bad <- diag %>% dplyr::filter(warning != "OK" | !valid_growth_fit | !is.finite(R0estimate))
    if (nrow(bad) == 0) return("No warnings. All trajectories produced finite R0estimate values.")
    paste(
      apply(bad, 1, function(row) {
        paste0("Series: ", row[["series"]],
               " | R0estimate: ", row[["R0estimate"]],
               " | r: ", row[["r"]],
               " | warning: ", row[["warning"]])
      }),
      collapse = "\n"
    )
  })

  observeEvent(input$csvgen_use_in_app, {
    req(csvgen_generated())
    df <- csvgen_generated()$output_csv
    df <- df %>%
      dplyr::mutate(
        R0estimate = suppressWarnings(as.numeric(R0estimate)),
        analytic_R0 = NA_real_,
        numeric_secondary_R0 = R0estimate,
        cases_increase = as.integer(R0estimate > 1),
        max_I = NA_real_,
        final_size = NA_real_,
        threshold_binary = as.integer(R0estimate > input$r0_cutoff),
        R0_bracket_analytic = NA,
        R0_bracket_numeric = make_r0_brackets(R0estimate, step = as.numeric(input$multiclass_step))
      )
    sim_data(df)
    binary_fit(NULL)
    multi_fit(NULL)
    updateRadioButtons(session, "data_mode", selected = "upload")
    showNotification("Generated CSV is now loaded as the active uploaded REEBLRA dataset. You can fit the binary or multiclass model.", type = "message", duration = 8)
  })

  output$csvgen_download_csv <- downloadHandler(
    filename = function() paste0("REEBLRA_upload_ready_R0estimate_parameters_", Sys.Date(), ".csv"),
    content = function(file) {
      req(csvgen_generated())
      readr::write_csv(csvgen_generated()$output_csv, file)
    }
  )

  output$csvgen_download_diagnostics <- downloadHandler(
    filename = function() paste0("growth_rate_R0_diagnostics_", Sys.Date(), ".csv"),
    content = function(file) {
      req(csvgen_generated())
      readr::write_csv(csvgen_generated()$diagnostics, file)
    }
  )


  custom_model_inputs <- reactive({
    eq_info <- custom_parse_equations(input$custom_equations, max_equations = 30)
    aux_info <- custom_parse_auxiliary_equations(input$custom_auxiliary_equations, max_equations = 30)
    init_state <- custom_parse_named_numeric_lines(input$custom_initial_values, value_name = "initial value")
    param_table <- custom_parse_parameter_ranges(input$custom_parameter_ranges, max_params = 50)
    custom_validate_model_inputs(eq_info, aux_info, init_state, param_table)
  })

  output$custom_nonnegative_state_controls <- renderUI({
    out <- tryCatch(custom_parse_equations(input$custom_equations, max_equations = 30), error = function(e) NULL)
    if (is.null(out) || length(out$states) == 0) {
      return(helpText("Enter valid custom state equations to choose which states are forced nonnegative."))
    }
    checkboxGroupInput(
      "custom_nonnegative_states",
      "State variables forced to be nonnegative",
      choices = out$states,
      selected = out$states,
      inline = TRUE
    )
  })

  custom_param_ranges <- reactive({
    custom_model_inputs()$param_table
  })

  output$custom_lhs_recommendation <- renderText({
    out <- tryCatch({
      npar <- nrow(custom_param_ranges())
      paste0("Detected ", npar, " custom parameter(s). Suggested for stable screening: 400 × parameters = ", 400 * npar,
             " simulations. Default shown: 400 for quick runs.")
    }, error = function(e) {
      paste0("Custom ODE input note: ", e$message)
    })
    out
  })

  uploaded_raw <- reactive({
    req(input$upload_r0_file)
    readr::read_csv(input$upload_r0_file$datapath, show_col_types = FALSE)
  })

  output$upload_column_controls <- renderUI({
    req(input$upload_r0_file)
    df <- uploaded_raw()
    num_cols <- names(df)[sapply(df, is.numeric)]
    tagList(
      selectInput("upload_r0_col", "R0 estimate column", choices = names(df), selected = if ("R0estimate" %in% names(df)) "R0estimate" else names(df)[1]),
      helpText(paste0("Detected ", length(setdiff(num_cols, input$upload_r0_col)), " numeric candidate parameter columns. Keep this at <= 50 for stable fitting and interpretation."))
    )
  })

  current_param_ranges <- reactive({
    base <- param_ranges()
    out <- base
    for (i in seq_len(nrow(out))) {
      p <- out$parameter[i]
      min_id <- paste0("min_", p)
      max_id <- paste0("max_", p)
      min_val <- input[[min_id]]
      max_val <- input[[max_id]]
      if (!is.null(min_val) && !is.na(min_val)) out$min[i] <- as.numeric(min_val)
      if (!is.null(max_val) && !is.na(max_val)) out$max[i] <- as.numeric(max_val)
    }
    # Guard against accidental min >= max.
    bad <- which(out$min >= out$max)
    if (length(bad) > 0) {
      out$max[bad] <- out$min[bad] + 1e-6
    }
    out
  })

  observeEvent(input$model_name, {
    param_ranges(model_parameter_table(input$model_name))
    sim_data(NULL)
    binary_fit(NULL)
    multi_fit(NULL)
  })

  output$parameter_range_controls <- renderUI({
    df <- param_ranges()
    tagList(lapply(seq_len(nrow(df)), function(i) {
      p <- df$parameter[i]
      fluidRow(
        column(12, tags$b(p)),
        column(6, numericInput(paste0("min_", p), "min", value = df$min[i], step = signif((df$max[i] - df$min[i]) / 20, 2))),
        column(6, numericInput(paste0("max_", p), "max", value = df$max[i], step = signif((df$max[i] - df$min[i]) / 20, 2)))
      )
    }))
  })

  output$formula_text <- renderText({
    if (identical(input$data_mode, "upload")) {
      "Uploaded dataset: analytic R0/Re is not available. The app uses the uploaded R0estimate column."
    } else if (identical(input$data_mode, "custom_ode")) {
      "Custom ODE dataset: analytic R0/Re is not derived automatically. The app uses the simulated numerical R0estimate from early growth of the selected state or auxiliary/output equation."
    } else {
      r0_formula_text(input$model_name)
    }
  })

  output$formula_text2 <- renderText({
    if (identical(input$data_mode, "upload")) {
      "Uploaded dataset: analytic R0/Re is not available."
    } else if (identical(input$data_mode, "custom_ode")) {
      "Custom ODE dataset: analytic R0/Re is not available automatically. Non-ODE auxiliary/output equations can be used for growth fitting."
    } else {
      r0_formula_text(input$model_name)
    }
  })

  output$ode_text <- renderText({
    if (identical(input$data_mode, "upload")) {
      "Uploaded dataset mode: no built-in ODE is used. The uploaded table should contain R0estimate and parameter columns."
    } else if (identical(input$data_mode, "custom_ode")) {
      out <- tryCatch(custom_model_inputs(), error = function(e) NULL)
      if (is.null(out)) {
        custom_ode_equation_text()
      } else {
        nn_states <- input$custom_nonnegative_states
        if (is.null(nn_states)) nn_states <- out$eq_info$states
        custom_ode_equation_text(out$eq_info, out$aux_info, input$custom_target_state, input$custom_conversion_method, nn_states)
      }
    } else {
      ode_equation_text(input$model_name)
    }
  })

  output$param_def_table <- renderDT({
    if (identical(input$data_mode, "upload")) {
      return(DT::datatable(data.frame(Note = "Parameter definitions are not available for uploaded datasets; use your own variable documentation."), rownames = FALSE))
    }
    if (identical(input$data_mode, "custom_ode")) {
      return(DT::datatable(data.frame(Note = "Parameter definitions are user-defined in custom ODE mode. Use the parameter ranges box as your model documentation."), rownames = FALSE))
    }
    DT::datatable(parameter_definition_table(input$model_name),
                  rownames = FALSE,
                  options = list(pageLength = 100, dom = "t", scrollX = TRUE))
  })

  output$lhs_recommendation <- renderText({
    npar <- nrow(current_param_ranges())
    paste0("Suggested for stable screening: 400 × number of parameters = ", 400 * npar,
           " simulations for this model (", npar, " parameters). Default shown: 400 for quick runs.")
  })

  output$param_table <- renderDT({
    if (identical(input$data_mode, "custom_ode")) {
      df <- sim_data()
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        custom_tbl <- tryCatch(custom_param_ranges(), error = function(e) data.frame(Note = e$message, check.names = FALSE))
        return(DT::datatable(custom_tbl, rownames = FALSE, options = list(pageLength = 100, dom = "t", scrollX = TRUE)))
      }
      predictors <- setdiff(names(df), c(
        "R0estimate", "analytic_R0", "numeric_secondary_R0", "cases_increase", "max_I", "final_size",
        "early_growth_r", "early_growth_r_lwr95", "early_growth_r_upr95", "early_growth_R2",
        "early_growth_n", "early_growth_threshold_used", "early_growth_window_end", "early_growth_valid", "early_growth_warning",
        "R0_bracket_analytic", "R0_bracket_numeric", "threshold_binary"
      ))
      summary_tbl <- data.frame(
        parameter = predictors,
        min = sapply(df[predictors], min, na.rm = TRUE),
        max = sapply(df[predictors], max, na.rm = TRUE),
        mean = sapply(df[predictors], mean, na.rm = TRUE),
        sd = sapply(df[predictors], stats::sd, na.rm = TRUE),
        row.names = NULL,
        check.names = FALSE
      )
      return(datatable(summary_tbl, rownames = FALSE,
                       options = list(pageLength = 100, dom = "t", scrollX = TRUE)))
    }
    if (identical(input$data_mode, "upload")) {
      df <- sim_data()
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(DT::datatable(note_table(data_ready_message()), rownames = FALSE, options = list(dom = "t")))
      }
      predictors <- setdiff(names(df), c(
        "R0estimate", "analytic_R0", "numeric_secondary_R0", "cases_increase", "max_I", "final_size",
        "early_growth_r", "early_growth_r_lwr95", "early_growth_r_upr95", "early_growth_R2",
        "early_growth_n", "early_growth_threshold_used", "early_growth_window_end", "early_growth_valid", "early_growth_warning",
        "R0_bracket_analytic", "R0_bracket_numeric", "threshold_binary"
      ))
      summary_tbl <- data.frame(
        parameter = predictors,
        min = sapply(df[predictors], min, na.rm = TRUE),
        max = sapply(df[predictors], max, na.rm = TRUE),
        mean = sapply(df[predictors], mean, na.rm = TRUE),
        sd = sapply(df[predictors], stats::sd, na.rm = TRUE),
        row.names = NULL,
        check.names = FALSE
      )
      datatable(summary_tbl, rownames = FALSE,
                options = list(pageLength = 100, dom = "t", scrollX = TRUE))
    } else {
      datatable(
        current_param_ranges(),
        editable = TRUE,
        rownames = FALSE,
        options = list(pageLength = 100, dom = "t")
      )
    }
  })

  observeEvent(input$param_table_cell_edit, {
    if (identical(input$data_mode, "upload") || identical(input$data_mode, "custom_ode")) return(NULL)
    info <- input$param_table_cell_edit
    df <- param_ranges()
    i <- info$row
    j <- info$col + 1
    if (names(df)[j] %in% c("min", "max")) {
      df[i, j] <- as.numeric(info$value)
      param_ranges(df)
    }
  })

  output$param_note <- renderText({
    if (identical(input$data_mode, "upload")) {
      if (is.null(sim_data())) return(data_ready_message())
      paste(
        "Uploaded R0estimate table mode. There are no ODE parameter ranges to edit.\n",
        "The table above summarizes the uploaded numeric parameter columns used as predictors.\n",
        "Analytic R0/Re is not available because the uploaded data may come from any model or simulator.\n",
        "Binary threshold uses R0estimate > cutoff; multiclass uses R0estimate brackets."
      )
    } else if (identical(input$data_mode, "custom_ode")) {
      out <- tryCatch(custom_model_inputs(), error = function(e) NULL)
      if (is.null(out)) return(data_ready_message())
      nn_states <- input$custom_nonnegative_states
      if (is.null(nn_states)) nn_states <- out$eq_info$states
      nn_states <- intersect(nn_states, out$eq_info$states)
      paste(
        "Custom user-defined ODE mode.\n",
        "States detected:", paste(out$eq_info$states, collapse = ", "), "\n",
        "Auxiliary/output variables detected:", ifelse(length(out$aux_info$variables) > 0, paste(out$aux_info$variables, collapse = ", "), "None"), "\n",
        "Parameters detected:", paste(out$param_table$parameter, collapse = ", "), "\n",
        "Target infectious/output state:", input$custom_target_state, "\n",
        "State variables forced to be nonnegative:", ifelse(length(nn_states) > 0, paste(nn_states, collapse = ", "), "None"), "\n",
        "Analytic R0/Re is not available automatically; binary and multiclass targets use numerical R0estimate from early growth."
      )
    } else {
      paste(
        "Model:", input$model_name, "\n",
        "Formula:", r0_formula_text(input$model_name), "\n",
        "Parameter ranges currently used by the simulator:\n",
        paste(capture.output(print(current_param_ranges())), collapse = "\n"), "\n",
        "Tip: choose wider ranges if all runs fall into the same threshold class."
      )
    }
  })

  observeEvent(input$run_sim, {
    set.seed(input$seed)

    if (identical(input$data_mode, "upload")) {
      req(input$upload_r0_file)
      r0_col <- if (!is.null(input$upload_r0_col)) input$upload_r0_col else "R0estimate"
      df <- tryCatch({
        prepare_uploaded_r0_data(input$upload_r0_file$datapath, r0_col = r0_col, max_params = 50, bracket_step = as.numeric(input$multiclass_step))
      }, error = function(e) {
        showNotification(paste("Upload error:", e$message), type = "error", duration = 12)
        NULL
      })
      req(df)
      df <- df %>%
        mutate(
          threshold_binary = as.integer(R0estimate > input$r0_cutoff),
          R0_bracket_numeric = make_r0_brackets(R0estimate, step = as.numeric(input$multiclass_step))
        )
      showNotification("Uploaded R0estimate table loaded and prepared.", type = "message", duration = 5)
    } else if (identical(input$data_mode, "custom_ode")) {
      custom_inputs <- tryCatch(custom_model_inputs(), error = function(e) {
        showNotification(paste("Custom ODE input error:", e$message), type = "error", duration = 12)
        NULL
      })
      req(custom_inputs)
      selected_nonnegative_states <- input$custom_nonnegative_states
      if (is.null(selected_nonnegative_states)) selected_nonnegative_states <- custom_inputs$eq_info$states
      selected_nonnegative_states <- intersect(selected_nonnegative_states, custom_inputs$eq_info$states)
      withProgress(message = "Running custom ODE simulations", value = 0, {
        n_sims <- as.integer(round(input$custom_n_samples))
        if (is.na(n_sims) || n_sims < 1) {
          stop("Number of LHS simulations must be a positive integer.")
        }
        df <- simulate_custom_many(
          param_table = custom_inputs$param_table,
          eq_info = custom_inputs$eq_info,
          aux_info = custom_inputs$aux_info,
          init_state = custom_inputs$init_state,
          n = n_sims,
          tmax = input$custom_tmax,
          dt = input$custom_dt,
          target_state = input$custom_target_state,
          nonnegative_states = selected_nonnegative_states,
          growth_epsilon = input$custom_growth_epsilon,
          conversion_method = input$custom_conversion_method,
          gamma_removal = input$custom_gamma_removal,
          sigma_latency = input$custom_sigma_latency,
          generation_time = input$custom_generation_time,
          nonpositive_rule = input$custom_nonpositive_rule,
          min_positive_value = input$custom_min_positive_value,
          growth_window_first_n = input$custom_growth_window_first_n,
          min_points = 5,
          show_progress = TRUE
        )
      })
      df <- df %>%
        mutate(
          threshold_binary = as.integer(R0estimate > input$r0_cutoff),
          R0_bracket_analytic = NA,
          R0_bracket_numeric = make_r0_brackets(R0estimate, step = as.numeric(input$multiclass_step))
        )
      showNotification("Custom ODE simulations finished and prepared for REEBLRA fitting.", type = "message", duration = 6)
    } else {
      withProgress(message = "Running ODE simulations", value = 0, {
        n_sims <- as.integer(round(input$n_samples))
        if (is.na(n_sims) || n_sims < 1) {
          stop("Number of LHS simulations must be a positive integer.")
        }
        df <- simulate_many(
          model_name = input$model_name,
          param_table = current_param_ranges(),
          n = n_sims,
          tmax = input$tmax,
          dt = input$dt,
          I0 = input$I0,
          growth_epsilon = input$growth_epsilon,
          susceptible_threshold = input$susceptible_threshold,
          show_progress = TRUE
        )
      })

      # Create threshold outcome selected by user.
      df <- df %>%
        mutate(
          threshold_binary = dplyr::case_when(
            input$threshold_source == "cases_increase" ~ as.integer(cases_increase),
            input$threshold_source == "analytic" ~ as.integer(analytic_R0 > input$r0_cutoff),
            input$threshold_source == "numeric" ~ as.integer(numeric_secondary_R0 > input$r0_cutoff),
            TRUE ~ as.integer(cases_increase)
          ),
          R0_bracket_analytic = make_r0_brackets(analytic_R0, step = as.numeric(input$multiclass_step)),
          R0_bracket_numeric = make_r0_brackets(numeric_secondary_R0, step = as.numeric(input$multiclass_step))
        )
    }

    sim_data(df)
    binary_fit(NULL)
    multi_fit(NULL)
  })

  sim_with_current_threshold <- reactive({
    df <- sim_data()
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      return(NULL)
    }
    if (identical(input$data_mode, "upload") || identical(input$data_mode, "custom_ode")) {
      df$R0estimate <- suppressWarnings(as.numeric(df$R0estimate))
      df$numeric_secondary_R0 <- df$R0estimate
      
      df %>%
        mutate(
          threshold_binary = as.integer(R0estimate > input$r0_cutoff),
          R0_bracket_numeric = make_r0_brackets(
            R0estimate,
            step = as.numeric(input$multiclass_step)
          )
        )
    } else {
      df %>%
        mutate(
          threshold_binary = dplyr::case_when(
            input$threshold_source == "cases_increase" ~ as.integer(cases_increase),
            input$threshold_source == "analytic" ~ as.integer(analytic_R0 > input$r0_cutoff),
            input$threshold_source == "numeric" ~ as.integer(numeric_secondary_R0 > input$r0_cutoff),
            TRUE ~ as.integer(cases_increase)
          ),
          R0_bracket_analytic = make_r0_brackets(analytic_R0, step = as.numeric(input$multiclass_step)),
          R0_bracket_numeric = make_r0_brackets(numeric_secondary_R0, step = as.numeric(input$multiclass_step))
        )
    }
  })

  current_predictors <- reactive({
    if (identical(input$data_mode, "upload") || identical(input$data_mode, "custom_ode")) {
      df <- sim_with_current_threshold()
      if (is.null(df)) return(character(0))
      setdiff(names(df), c(
        "R0estimate", "analytic_R0", "numeric_secondary_R0", "cases_increase", "max_I", "final_size",
        "early_growth_r", "early_growth_r_lwr95", "early_growth_r_upr95", "early_growth_R2",
        "early_growth_n", "early_growth_threshold_used", "early_growth_window_end", "early_growth_valid", "early_growth_warning",
        "R0_bracket_analytic", "R0_bracket_numeric", "threshold_binary"
      ))
    } else {
      current_param_ranges()$parameter
    }
  })

  output$sim_summary <- renderTable({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(note_table(data_ready_message()))
    tibble::tibble(
      Quantity = c("Simulation runs", "Parameters", "Threshold positives", "Threshold negatives",
                   "Mean analytic R0/Re", "Mean numeric growth-rate R0estimate"),
      Value = c(
        nrow(df),
        paste(current_predictors(), collapse = ", "),
        sum(df$threshold_binary == 1, na.rm = TRUE),
        sum(df$threshold_binary == 0, na.rm = TRUE),
        round(mean(df$analytic_R0, na.rm = TRUE), 4),
        round(mean(df$numeric_secondary_R0, na.rm = TRUE), 4)
      )
    )
  })

  output$r0_compare <- renderTable({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(note_table(data_ready_message()))
    tibble::tibble(
      Metric = c("Correlation: analytic vs numeric growth-rate R0estimate",
                 "MAE: analytic vs numeric growth-rate R0estimate",
                 "AUC: analytic R0 predicts case increase",
                 "AUC: numeric growth-rate R0estimate predicts case increase"),
      Value = c(
        ifelse(all(is.na(df$analytic_R0)), NA_real_, round(cor(df$analytic_R0, df$numeric_secondary_R0, use = "complete.obs"), 4)),
        ifelse(all(is.na(df$analytic_R0)), NA_real_, round(mean(abs(df$analytic_R0 - df$numeric_secondary_R0), na.rm = TRUE), 4)),
        round(tryCatch(as.numeric(pROC::auc(pROC::roc(df$cases_increase, df$analytic_R0, quiet = TRUE))), error = function(e) NA_real_), 4),
        round(tryCatch(as.numeric(pROC::auc(pROC::roc(df$cases_increase, df$numeric_secondary_R0, quiet = TRUE))), error = function(e) NA_real_), 4)
      )
    )
  })

  output$sim_table <- renderDT({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(DT::datatable(note_table(data_ready_message()), rownames = FALSE, options = list(dom = "t")))
    datatable(df, options = list(scrollX = TRUE, pageLength = 100))
  })

  output$vif_table <- renderDT({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(DT::datatable(note_table(data_ready_message()), rownames = FALSE, options = list(dom = "t")))
    DT::datatable(
      vif_table_from_data(df, current_predictors()),
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 100)
    )
  })

  output$corr_matrix_table <- renderDT({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(DT::datatable(note_table(data_ready_message()), rownames = FALSE, options = list(dom = "t")))
    cor_mat <- correlation_matrix_from_data(df, current_predictors())
    if (is.data.frame(cor_mat) && "Note" %in% names(cor_mat)) {
      return(DT::datatable(cor_mat, rownames = FALSE))
    }
    DT::datatable(round(as.data.frame(cor_mat), 4),
                  options = list(scrollX = TRUE, pageLength = 100))
  })

  output$corr_matrix_plot <- renderPlotly({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(note_plotly(data_ready_message()))
    cor_long <- correlation_long_from_data(df, current_predictors())
    if ("Note" %in% names(cor_long)) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = cor_long$Note[1]))
    }
    plot_ly(
      cor_long,
      x = ~Parameter_1,
      y = ~Parameter_2,
      z = ~Correlation,
      type = "heatmap",
      colorscale = "RdBu",
      zmin = -1,
      zmax = 1,
      text = ~paste0(Parameter_1, " vs ", Parameter_2, "<br>Correlation = ", round(Correlation, 4)),
      hoverinfo = "text"
    ) %>%
      layout(
        title = "Parameter correlation matrix",
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "")
      )
  })

  observeEvent(input$fit_binary, {
    df <- sim_with_current_threshold()
    if (is.null(df)) {
      showNotification(data_ready_message(), type = "warning", duration = 8)
      return(NULL)
    }
    predictors <- current_predictors()
    if (length(predictors) == 0) {
      showNotification("No numeric parameter predictors are available for fitting.", type = "warning", duration = 8)
      return(NULL)
    }

    binary_fit(NULL)
    wait_id <- showNotification(
      "Fitting binary Bayesian logistic R0 model. Please wait; results may appear with a short delay, especially with MCMC or k-fold validation.",
      type = "message", duration = NULL
    )
    on.exit(removeNotification(wait_id), add = TRUE)

    fit <- tryCatch({
      withProgress(message = "Fitting binary Bayesian logistic model — please wait", value = 0.1, {
        fit_binary_bayes(
          df = df,
          predictors = predictors,
          outcome = "threshold_binary",
          chains = input$chains,
          iter = input$iter,
          seed = input$seed,
          algorithm = input$stan_algorithm,
          k_folds = input$k_folds
        )
      })
    }, error = function(e) {
      showNotification(paste("Bayesian fitting error:", e$message), type = "error", duration = 12)
      NULL
    })

    if (!is.null(fit)) {
      binary_fit(fit)
      showNotification("Binary Bayesian logistic R0 model finished.", type = "message", duration = 5)
    }
  })

  output$binary_equation <- renderText({
    if (is.null(binary_fit())) return(binary_ready_message())
    fit <- binary_fit()
    eq <- fit$equation
    formula_text <- safe_deparse_formula(fit$formula)

    paste(
      "Outcome:", fit$outcome, "\n",
      "Bayesian algorithm:", fit$algorithm, "\n",
      "K-fold validation setting:", fit$k_folds, "\n",
      "Model formula used internally:", formula_text, "\n\n",
      "Why the equation may be long:\n",
      "If the uploaded dataset has many numeric parameters, REEBLRA fits all selected numeric predictors together. Long equations can look visually overwhelming, especially on the original scale. The standardized equation is usually easier to interpret because each coefficient is comparable across parameters.\n\n",
      "TOP STANDARDIZED TERMS BY ABSOLUTE POSTERIOR MEAN:\n",
      top_binary_terms_text(fit, top_n = 10), "\n\n",
      "STANDARDIZED VERSION USED FOR FITTING:\n",
      "Predictors are standardized internally: z_parameter = (parameter - mean) / sd.\n",
      collapse_long_equation(eq$logit_standardized), "\n",
      collapse_long_equation(eq$probability_standardized), "\n",
      collapse_long_equation(eq$threshold_standardized), "\n\n",
      "ORIGINAL-SCALE EQUATION:\n",
      "This is mathematically useful for direct substitution, but coefficients can look very large when parameters have very small scales or very different units. Prefer the standardized equation and coefficient table for interpretation.\n",
      collapse_long_equation(eq$logit_original), "\n",
      collapse_long_equation(eq$probability_original), "\n",
      collapse_long_equation(eq$threshold_original), "\n\n",
      "Interpretation: this is the in silico Bayesian logistic R0-threshold equation. For R0 cutoff = 1, p > 0.5 means the app estimates R0 > 1 for that parameter set. Positive coefficients increase the probability of crossing the threshold; negative coefficients decrease it."
    )
  })

  output$binary_metrics <- renderTable({
    if (is.null(binary_fit())) return(note_table(binary_ready_message()))
    binary_fit()$metrics
  })

  output$mcmc_diagnostics <- renderDT({
    if (is.null(binary_fit())) return(DT::datatable(note_table(binary_ready_message()), rownames = FALSE, options = list(dom = "t")))
    DT::datatable(mcmc_diagnostics_table(binary_fit()), options = list(scrollX = TRUE, pageLength = 100))
  })

  output$binary_trial_inputs <- renderUI({
    if (is.null(binary_fit())) return(helpText(binary_ready_message()))
    bf <- binary_fit()
    defaults <- if (!is.null(bf$centers)) bf$centers else NULL
    make_trial_inputs("binary_trial_", bf$predictors, defaults)
  })

  output$binary_trial_output <- renderText({
    if (is.null(binary_fit())) return(binary_ready_message())
    bf <- binary_fit()
    newd <- trial_newdata(input, "binary_trial_", bf$predictors)
    if (any(!stats::complete.cases(newd))) return("Please enter complete parameter values.")
    pred <- predict_binary_trial(bf, newd)
    analytic <- safe_analytic_for_trial(input$model_name, newd, input$data_mode)
    pred_text <- paste0(
      "Predicted probability of selected threshold being TRUE: ", round(pred$Predicted_probability_threshold_TRUE[1], 6), "\n",
      "Predicted binary interpretation: ", pred$Predicted_class[1], "\n",
      "Binary cutoff used: ", input$r0_cutoff
    )
    format_trial_comparison(pred_text, analytic)
  })

  output$binary_coef <- renderDT({
    if (is.null(binary_fit())) return(DT::datatable(note_table(binary_ready_message()), rownames = FALSE, options = list(dom = "t")))
    bf <- binary_fit()
    coefs <- posterior_coef_table(bf$fit, z_predictors = bf$z_predictors)
    datatable(coefs, options = list(scrollX = TRUE, pageLength = 100))
  })

  output$binary_roc <- renderPlot({
    if (is.null(binary_fit())) {
      plot.new()
      title("Binary ROC curve unavailable")
      text(0.5, 0.5, binary_ready_message())
      return(invisible(NULL))
    }
    fit <- binary_fit()

    y_app <- fit$data[[fit$outcome]]
    roc_app <- tryCatch(pROC::roc(y_app, fit$prob, quiet = TRUE), error = function(e) NULL)

    if (is.null(roc_app)) {
      plot.new()
      title("ROC curve unavailable")
      text(0.5, 0.5, "ROC could not be computed. Check whether both binary classes are present.")
      return(invisible(NULL))
    }

    plot(
      roc_app,
      main = "Binary ROC curve: apparent/in-sample and k-fold",
      legacy.axes = TRUE,
      lwd = 2
    )
    legend_labels <- paste0("Apparent / in-sample AUC = ", round(as.numeric(pROC::auc(roc_app)), 4))
    legend_lty <- 1
    legend_lwd <- 2

    if (!is.null(fit$cv_predictions) && nrow(fit$cv_predictions) > 0 &&
        ".observed" %in% names(fit$cv_predictions) &&
        ".cv_predicted_prob" %in% names(fit$cv_predictions) &&
        length(unique(fit$cv_predictions$.observed)) == 2) {
      roc_cv <- tryCatch(
        pROC::roc(fit$cv_predictions$.observed, fit$cv_predictions$.cv_predicted_prob, quiet = TRUE),
        error = function(e) NULL
      )
      if (!is.null(roc_cv)) {
        plot(roc_cv, add = TRUE, legacy.axes = TRUE, lwd = 2, lty = 2)
        legend_labels <- c(
          legend_labels,
          paste0(fit$k_folds, "-fold CV AUC = ", round(as.numeric(pROC::auc(roc_cv)), 4))
        )
        legend_lty <- c(legend_lty, 2)
        legend_lwd <- c(legend_lwd, 2)
      }
    }

    legend(
      "bottomright",
      legend = legend_labels,
      lty = legend_lty,
      lwd = legend_lwd,
      bty = "n"
    )
  })


  output$binary_axis_controls <- renderUI({
    if (is.null(binary_fit())) return(helpText(binary_ready_message()))
    bf <- binary_fit()
    choices <- bf$predictors
    if (length(choices) == 0) return(NULL)
    fluidRow(
      column(6, selectInput("binary_x_axis", "X-axis parameter", choices = choices, selected = choices[1])),
      column(6, selectInput("binary_y_axis", "Y-axis parameter", choices = choices, selected = ifelse(length(choices) >= 2, choices[2], choices[1])))
    )
  })

  output$binary_plotly <- renderPlotly({
    if (is.null(binary_fit())) return(note_plotly(binary_ready_message()))
    fit <- binary_fit()
    d <- fit$data
    d$predicted_prob <- fit$prob
    d$predicted_class <- ifelse(d$predicted_prob >= 0.5, "above cutoff", "below cutoff")
    predictors <- fit$predictors
    xvar <- if (!is.null(input$binary_x_axis) && input$binary_x_axis %in% predictors) input$binary_x_axis else predictors[1]
    yvar <- if (!is.null(input$binary_y_axis) && input$binary_y_axis %in% predictors) input$binary_y_axis else ifelse(length(predictors) >= 2, predictors[2], predictors[1])

    hover_text <- apply(d[, fit$predictors, drop = FALSE], 1, function(row) {
      paste(paste(names(row), round(as.numeric(row), 5), sep = " = "), collapse = "<br>")
    })
    d$hover_text <- paste0(
      hover_text,
      "<br>Predicted probability = ", round(d$predicted_prob, 5),
      "<br>Predicted class = ", d$predicted_class
    )

    if (length(predictors) >= 2 && !identical(xvar, yvar)) {
      plot_ly(
        d,
        x = as.formula(paste0("~", xvar)),
        y = as.formula(paste0("~", yvar)),
        color = ~predicted_prob,
        text = ~hover_text,
        hoverinfo = "text",
        type = "scatter",
        mode = "markers",
        marker = list(size = 8)
      ) %>%
        layout(
          title = "2D predicted threshold probability scatter: color = predicted probability",
          xaxis = list(title = xvar),
          yaxis = list(title = yvar)
        )
    } else {
      plot_ly(d, x = as.formula(paste0("~", xvar)), y = ~predicted_prob,
              text = ~hover_text, hoverinfo = "text",
              type = "scatter", mode = "markers") %>%
        layout(title = "Predicted threshold probability", xaxis = list(title = xvar), yaxis = list(title = "Predicted probability"))
    }
  })


  output$binary_ordinary_equation <- renderText({
    if (is.null(binary_fit())) return(binary_ready_message())
    bf <- binary_fit()
    if (is.null(bf$ordinary_fit) || is.null(bf$ordinary_equation)) {
      return("Ordinary binary logistic regression could not be fitted for the current data.")
    }
    eq <- bf$ordinary_equation
    paste(
      "ORDINARY LOGISTIC REGRESSION ON STANDARDIZED PREDICTORS:\n",
      collapse_long_equation(eq$eta_standardized), "\n",
      collapse_long_equation(eq$probability_standardized), "\n\n",
      "ORDINARY LOGISTIC REGRESSION ON ORIGINAL PARAMETER SCALE:\n",
      collapse_long_equation(eq$eta_original), "\n",
      collapse_long_equation(eq$probability_original), "\n\n",
      "Classification rule: p > 0.5 is equivalent to eta > 0."
    )
  })

  observeEvent(input$fit_multi, {
    df <- sim_with_current_threshold()
    if (is.null(df)) {
      showNotification(data_ready_message(), type = "warning", duration = 8)
      return(NULL)
    }
    predictors <- current_predictors()
    if (length(predictors) == 0) {
      showNotification("No numeric parameter predictors are available for fitting.", type = "warning", duration = 8)
      return(NULL)
    }

    multi_fit(NULL)
    wait_id <- showNotification(
      "Fitting multiclass R0 bracket model. Please wait; results may appear with a short delay when many thresholds or k-fold validation are used.",
      type = "message", duration = NULL
    )
    on.exit(removeNotification(wait_id), add = TRUE)

    fit <- tryCatch({
      withProgress(message = "Fitting multiclass R0 bracket model — please wait", value = 0.1, {
        fit_multiclass_model(
          df = df,
          predictors = predictors,
          outcome_class = if (identical(input$data_mode, "upload") || identical(input$data_mode, "custom_ode")) "R0_bracket_numeric" else input$multiclass_source,
          bracket_step = as.numeric(input$multiclass_step),
          seed = input$seed,
          k_folds = input$k_folds
        )
      })
    }, error = function(e) {
      showNotification(paste("Multiclass fitting error:", e$message), type = "error", duration = 12)
      NULL
    })

    if (!is.null(fit)) {
      multi_fit(fit)
      showNotification("Multiclass R0 bracket model finished.", type = "message", duration = 5)
    }
  })

  output$multi_method <- renderText({
    if (is.null(multi_fit())) return(multi_ready_message())
    paste(
      multi_fit()$method,
      "\nOutcome:", multi_fit()$outcome_class,
      "\nPredictors are standardized internally.",
      "\nK-fold validation setting:", multi_fit()$k_folds,
      "\nObserved R0estimate range used for fitting:",
      paste0(
        round(min(multi_fit()$data$R0_for_successive, na.rm = TRUE), 6),
        " to ",
        round(max(multi_fit()$data$R0_for_successive, na.rm = TRUE), 6)
      ),
      "\nNote: this app uses successive binary logistic equations because R0 brackets are ordered. The model fits R0 > 1, then the next selected step-size threshold, and so on until R0 > 10. A finite equation is estimable only when the simulated data contain observations on both sides of that threshold. The cumulative probabilities are converted into bracket probabilities, and the predicted bracket is the one with the largest probability."
    )
  })

  output$multi_equation <- renderText({
    if (is.null(multi_fit())) return(multi_ready_message())
    multiclass_equation_text(multi_fit())
  })


  output$multi_coef_heatmap <- renderPlotly({
    if (is.null(multi_fit())) return(note_plotly(multi_ready_message()))
    heat <- multiclass_coefficient_heatmap_data(multi_fit())
    if (nrow(heat) == 0) {
      return(note_plotly("No multiclass coefficient matrix is available."))
    }

    threshold_order <- paste0(
      "R0 > ",
      vapply(as.character(multi_fit()$thresholds), format_bracket_number, character(1))
    )
    plot_ly(
      heat,
      x = ~Parameter,
      y = ~R0_level,
      z = ~Row_normalized,
      type = "heatmap",
      colorscale = list(
        list(0, "#1F4E79"),
        list(0.5, "#FFFFFF"),
        list(1, "#7A0019")
      ),
      zmin = -1,
      zmax = 1,
      xgap = 1,
      ygap = 1,
      text = ~paste0(
        "R0 level: ", R0_level,
        "<br>Parameter: ", Parameter,
        "<br>Standardized coefficient: ", ifelse(is.finite(Estimate), round(Estimate, 6), "not estimable"),
        "<br>Row-normalized coefficient: ", ifelse(is.finite(Row_normalized), round(Row_normalized, 6), "not estimable"),
        "<br>N at or below threshold: ", N_at_or_below,
        "<br>N above threshold: ", N_above,
        "<br>Status: ", Status
      ),
      hoverinfo = "text",
      colorbar = list(title = "Row-normalized<br>coefficient")
    ) %>%
      layout(
        title = "Standardized coefficient matrix across successive R0 thresholds",
        xaxis = list(title = "Parameter", tickangle = -45),
        yaxis = list(
          title = "R0 threshold level",
          categoryorder = "array",
          categoryarray = rev(threshold_order)
        ),
        margin = list(l = 110, b = 120)
      )
  })

  output$multi_coef <- renderDT({
    if (is.null(multi_fit())) return(DT::datatable(note_table(multi_ready_message()), rownames = FALSE, options = list(dom = "t")))
    datatable(multiclass_coef_table(multi_fit()), options = list(scrollX = TRUE, pageLength = 100))
  })

  output$multi_trial_inputs <- renderUI({
    if (is.null(multi_fit())) return(helpText(multi_ready_message()))
    mf <- multi_fit()
    defaults <- if (!is.null(mf$centers)) mf$centers else NULL
    make_trial_inputs("multi_trial_", mf$predictors, defaults)
  })

  output$multi_trial_output <- renderText({
    if (is.null(multi_fit())) return(multi_ready_message())
    mf <- multi_fit()
    newd <- trial_newdata(input, "multi_trial_", mf$predictors)
    if (any(!stats::complete.cases(newd))) return("Please enter complete parameter values.")
    pr <- predict_multiclass_trial(mf, newd)
    analytic <- safe_analytic_for_trial(input$model_name, newd, input$data_mode)
    probs_text <- paste(capture.output(print(round(pr$probabilities, 6))), collapse = "\n")
    pred_text <- paste0(
      "Predicted R0 bracket: ", pr$predicted_class[1], "\n\n",
      "Predicted probabilities by bracket:\n", probs_text
    )
    if (!is.na(analytic)) {
      pred_text <- paste0(pred_text, "\n\nAnalytic R0/Re for these parameter values: ", round(analytic, 6),
                          "\nAnalytic R0 bracket: ", as.character(make_r0_brackets(analytic, step = as.numeric(input$multiclass_step))))
    } else {
      pred_text <- paste0(pred_text, "\n\nAnalytic R0/Re: not available for uploaded datasets or this input.")
    }
    pred_text
  })

  output$multi_metrics <- renderTable({
    if (is.null(multi_fit())) return(note_table(multi_ready_message()))
    multiclass_metrics(multi_fit())
  })

  output$multi_class_distribution_plot <- renderPlotly({
    if (is.null(multi_fit())) return(note_plotly(multi_ready_message()))
    fit <- multi_fit()
    observed <- if (!is.null(fit$observed_bracket)) {
      factor(as.character(fit$observed_bracket), levels = fit$labels)
    } else {
      factor(as.character(fit$data[[fit$outcome_class]]))
    }
    predicted <- factor(as.character(fit$pred_class), levels = levels(observed))

    dist_df <- dplyr::bind_rows(
      as.data.frame(table(Bracket = observed), stringsAsFactors = FALSE) %>%
        dplyr::mutate(Distribution = "Observed"),
      as.data.frame(table(Bracket = predicted), stringsAsFactors = FALSE) %>%
        dplyr::mutate(Distribution = "Predicted")
    ) %>%
      dplyr::rename(Count = Freq) %>%
      dplyr::mutate(
        Bracket = as.character(Bracket),
        hover_text = paste0(Distribution, " bracket = ", Bracket, "<br>Count = ", Count)
      )

    plot_ly(
      dist_df,
      x = ~Bracket,
      y = ~Count,
      color = ~Distribution,
      type = "bar",
      text = ~hover_text,
      hoverinfo = "text"
    ) %>%
      layout(
        title = "Observed and predicted R0 bracket distribution",
        barmode = "group",
        xaxis = list(title = "R0 bracket", tickangle = -45),
        yaxis = list(title = "Number of simulations")
      )
  })

  output$multi_auc_plot <- renderPlotly({
    if (is.null(multi_fit())) return(note_plotly(multi_ready_message()))
    auc_df <- multiclass_metrics(multi_fit()) %>%
      dplyr::filter(Metric == "One-vs-rest AUC", Validation == "Apparent / in-sample", !Class %in% c("Overall", "Macro average", "Weighted average")) %>%
      dplyr::mutate(
        Class = as.character(Class),
        Value = as.numeric(Value),
        hover_text = paste0("R0 bracket = ", Class, "<br>One-vs-rest AUC = ", round(Value, 4))
      )

    auc_df_valid <- auc_df %>% dplyr::filter(!is.na(Value), is.finite(Value))
    if (nrow(auc_df_valid) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = "One-vs-rest AUC unavailable: some classes may have only one outcome level."))
    }

    plot_ly(
      auc_df_valid,
      x = ~Class,
      y = ~Value,
      type = "bar",
      text = ~hover_text,
      hoverinfo = "text"
    ) %>%
      layout(
        title = "One-vs-rest AUC per R0 bracket",
        xaxis = list(title = "R0 bracket", tickangle = -45),
        yaxis = list(title = "One-vs-rest AUC", range = c(0, 1)),
        shapes = list(
          list(
            type = "line",
            x0 = -0.5,
            x1 = nrow(auc_df_valid) - 0.5,
            y0 = 0.5,
            y1 = 0.5,
            line = list(dash = "dash")
          )
        )
      )
  })


  output$multi_cv_auc_plot <- renderPlotly({
    if (is.null(multi_fit())) return(note_plotly(multi_ready_message()))
    mf <- multi_fit()
    cv_label <- paste0(mf$k_folds, "-fold cross-validation")
    auc_df <- multiclass_metrics(mf) %>%
      dplyr::filter(Metric == "One-vs-rest AUC", Validation == cv_label, !Class %in% c("Overall", "Macro average", "Weighted average")) %>%
      dplyr::mutate(
        Class = as.character(Class),
        Value = as.numeric(Value),
        hover_text = paste0("R0 bracket = ", Class, "<br>K-fold one-vs-rest AUC = ", round(Value, 4))
      )

    auc_df_valid <- auc_df %>% dplyr::filter(!is.na(Value), is.finite(Value))
    if (nrow(auc_df_valid) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = "K-fold one-vs-rest AUC unavailable: set k >= 2 and ensure held-out predictions contain usable classes."))
    }

    plot_ly(
      auc_df_valid,
      x = ~Class,
      y = ~Value,
      type = "bar",
      text = ~hover_text,
      hoverinfo = "text"
    ) %>%
      layout(
        title = "K-fold one-vs-rest AUC per R0 bracket",
        xaxis = list(title = "R0 bracket", tickangle = -45),
        yaxis = list(title = "K-fold one-vs-rest AUC", range = c(0, 1)),
        shapes = list(
          list(
            type = "line",
            x0 = -0.5,
            x1 = nrow(auc_df_valid) - 0.5,
            y0 = 0.5,
            y1 = 0.5,
            line = list(dash = "dash")
          )
        )
      )
  })

  output$multi_confusion <- renderTable({
    if (is.null(multi_fit())) return(note_table(multi_ready_message()))
    fit <- multi_fit()
    observed_raw <- if (!is.null(fit$observed_bracket)) {
      as.character(fit$observed_bracket)
    } else {
      as.character(fit$data[[fit$outcome_class]])
    }
    observed <- factor(observed_raw, levels = fit$labels)
    predicted <- factor(as.character(fit$pred_class), levels = fit$labels)
    table(Observed = observed, Predicted = predicted)
  })


  output$multi_axis_controls <- renderUI({
    if (is.null(multi_fit())) return(helpText(multi_ready_message()))
    mf <- multi_fit()
    choices <- mf$predictors
    if (length(choices) == 0) return(NULL)
    fluidRow(
      column(6, selectInput("multi_x_axis", "X-axis parameter", choices = choices, selected = choices[1])),
      column(6, selectInput("multi_y_axis", "Y-axis parameter", choices = choices, selected = ifelse(length(choices) >= 2, choices[2], choices[1])))
    )
  })

  output$multi_plotly <- renderPlotly({
    if (is.null(multi_fit())) return(note_plotly(multi_ready_message()))
    fit <- multi_fit()
    d <- fit$data
    d$predicted_class <- fit$pred_class
    predictors <- fit$predictors
    xvar <- if (!is.null(input$multi_x_axis) && input$multi_x_axis %in% predictors) input$multi_x_axis else predictors[1]
    yvar <- if (!is.null(input$multi_y_axis) && input$multi_y_axis %in% predictors) input$multi_y_axis else ifelse(length(predictors) >= 2, predictors[2], predictors[1])

    hover_text <- apply(d[, fit$predictors, drop = FALSE], 1, function(row) {
      paste(paste(names(row), round(as.numeric(row), 5), sep = " = "), collapse = "<br>")
    })
    obs <- if (!is.null(fit$observed_bracket)) as.character(fit$observed_bracket) else as.character(fit$data[[fit$outcome_class]])
    d$hover_text <- paste0(
      hover_text,
      "<br>Observed bracket = ", obs,
      "<br>Predicted bracket = ", d$predicted_class
    )

    if (length(predictors) >= 2 && !identical(xvar, yvar)) {
      plot_ly(
        d,
        x = as.formula(paste0("~", xvar)),
        y = as.formula(paste0("~", yvar)),
        color = ~predicted_class,
        text = ~hover_text,
        hoverinfo = "text",
        type = "scatter",
        mode = "markers",
        marker = list(size = 8)
      ) %>%
        layout(
          title = "Predicted R0 bracket",
          xaxis = list(title = xvar),
          yaxis = list(title = yvar)
        )
    } else {
      plot_ly(
        d,
        x = as.formula(paste0("~", xvar)),
        color = ~predicted_class,
        text = ~hover_text,
        hoverinfo = "text",
        type = "histogram"
      ) %>%
        layout(title = "Predicted R0 bracket", xaxis = list(title = xvar))
    }
  })

  output$compare_plot <- renderPlotly({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(note_plotly(data_ready_message()))
    if (!is.null(binary_fit())) {
      bf <- binary_fit()
      # Match rows used by binary fit.
      d <- bf$data
      d$predicted_prob <- bf$prob
      if (all(is.na(d$analytic_R0))) {
        d$row_id <- seq_len(nrow(d))
        plot_ly(
          d,
          x = ~row_id,
          y = ~numeric_secondary_R0,
          z = ~predicted_prob,
          color = ~predicted_prob,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 3)
        ) %>%
          layout(
            title = "Uploaded R0estimate vs learned probability",
            scene = list(
              xaxis = list(title = "Row"),
              yaxis = list(title = "Uploaded R0estimate"),
              zaxis = list(title = "Learned threshold probability")
            )
          )
      } else {
        plot_ly(
          d,
          x = ~analytic_R0,
          y = ~numeric_secondary_R0,
          z = ~predicted_prob,
          color = ~predicted_prob,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 3)
        ) %>%
          layout(
            title = "Analytic R0 vs numeric growth-rate R0estimate vs learned probability",
            scene = list(
              xaxis = list(title = "Analytic R0/Re"),
              yaxis = list(title = "Numeric growth-rate R0estimate"),
              zaxis = list(title = "Learned threshold probability")
            )
          )
      }
    } else {
      if (all(is.na(df$analytic_R0))) {
        df$row_id <- seq_len(nrow(df))
        plot_ly(df, x = ~row_id, y = ~numeric_secondary_R0, color = ~factor(threshold_binary),
                type = "scatter", mode = "markers") %>%
          layout(title = "Uploaded R0estimate by row", xaxis = list(title = "Row"), yaxis = list(title = "R0estimate"))
      } else {
        plot_ly(
          df,
          x = ~analytic_R0,
          y = ~numeric_secondary_R0,
          color = ~factor(cases_increase),
          type = "scatter",
          mode = "markers"
        ) %>%
          layout(title = "Analytic R0/Re vs numerical growth-rate R0estimate")
      }
    }
  })

  output$compare_metrics <- renderTable({
    df <- sim_with_current_threshold()
    if (is.null(df)) return(note_table(data_ready_message()))
    out <- tibble::tibble(
      Metric = c(
        "Correlation analytic vs numeric growth-rate R0estimate",
        "MAE analytic vs numeric growth-rate R0estimate",
        "AUC learned probability vs selected threshold"
      ),
      Value = c(
        ifelse(all(is.na(df$analytic_R0)), NA_real_, round(cor(df$analytic_R0, df$numeric_secondary_R0, use = "complete.obs"), 4)),
        ifelse(all(is.na(df$analytic_R0)), NA_real_, round(mean(abs(df$analytic_R0 - df$numeric_secondary_R0), na.rm = TRUE), 4)),
        NA_real_
      )
    )
    if (!is.null(binary_fit())) {
      bf <- binary_fit()
      out$Value[3] <- round(tryCatch(as.numeric(pROC::auc(pROC::roc(bf$data[[bf$outcome]], bf$prob, quiet = TRUE))), error = function(e) NA_real_), 4)
    }
    out
  })


  output$sensitivity_table <- renderDT({
    if (is.null(binary_fit())) return(DT::datatable(note_table(binary_ready_message()), rownames = FALSE, options = list(dom = "t")))
    datatable(sensitivity_from_binary_fit(binary_fit()),
              rownames = FALSE,
              options = list(pageLength = 100, scrollX = TRUE))
  })

  output$sensitivity_plot <- renderPlotly({
    if (is.null(binary_fit())) return(note_plotly(binary_ready_message()))
    s <- sensitivity_from_binary_fit(binary_fit())
    score_col <- "Normalized_Bayesian_sensitivity_coefficient_median"
    lwr_col <- "Normalized_Bayesian_sensitivity_coefficient_lwr95"
    upr_col <- "Normalized_Bayesian_sensitivity_coefficient_upr95"
    if (!score_col %in% names(s)) {
      return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = "No sensitivity table available"))
    }
    s$Parameter <- factor(s$Parameter, levels = rev(s$Parameter))
    plot_ly(
      s,
      x = ~Normalized_Bayesian_sensitivity_coefficient_median,
      y = ~Parameter,
      type = "bar",
      orientation = "h",
      error_x = list(
        type = "data",
        symmetric = FALSE,
        array = ~pmax(0, Normalized_Bayesian_sensitivity_coefficient_upr95 - Normalized_Bayesian_sensitivity_coefficient_median),
        arrayminus = ~pmax(0, Normalized_Bayesian_sensitivity_coefficient_median - Normalized_Bayesian_sensitivity_coefficient_lwr95)
      ),
      text = ~paste0(
        "Parameter = ", Parameter,
        "<br>Median normalized coefficient = ", round(Normalized_Bayesian_sensitivity_coefficient_median, 5),
        "<br>95% CrI = [", round(Normalized_Bayesian_sensitivity_coefficient_lwr95, 5), ", ", round(Normalized_Bayesian_sensitivity_coefficient_upr95, 5), "]",
        "<br>Mean normalized coefficient = ", round(Normalized_Bayesian_sensitivity_coefficient_mean, 5),
        "<br>Transformation = tanh(standardized Bayesian coefficient)"
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        title = "Normalized Bayesian sensitivity coefficient",
        xaxis = list(title = "Median normalized Bayesian sensitivity coefficient", range = c(-1, 1), zeroline = TRUE),
        yaxis = list(title = "")
      )
  })

  # ---------------- Downloads ----------------

  output$download_sim <- downloadHandler(
    filename = function() paste0("r0_simulation_", gsub(" ", "_", ifelse(identical(input$data_mode, "upload"), "Uploaded_R0_Table", ifelse(identical(input$data_mode, "custom_ode"), "Custom_ODE_Model", input$model_name))), ".csv"),
    content = function(file) write.csv(sim_with_current_threshold(), file, row.names = FALSE)
  )

  output$download_sim2 <- downloadHandler(
    filename = function() paste0("r0_simulation_", gsub(" ", "_", ifelse(identical(input$data_mode, "upload"), "Uploaded_R0_Table", ifelse(identical(input$data_mode, "custom_ode"), "Custom_ODE_Model", input$model_name))), ".csv"),
    content = function(file) write.csv(sim_with_current_threshold(), file, row.names = FALSE)
  )

  output$download_binary_metrics <- downloadHandler(
    filename = function() "binary_bayesian_r0_metrics.csv",
    content = function(file) {
      req(binary_fit())
      write.csv(binary_fit()$metrics, file, row.names = FALSE)
    }
  )

  output$download_binary_metrics2 <- downloadHandler(
    filename = function() "binary_bayesian_r0_metrics.csv",
    content = function(file) {
      req(binary_fit())
      write.csv(binary_fit()$metrics, file, row.names = FALSE)
    }
  )

  output$download_binary_predictions <- downloadHandler(
    filename = function() "binary_bayesian_r0_predictions.csv",
    content = function(file) {
      req(binary_fit())
      bf <- binary_fit()
      d <- bf$data
      d$predicted_prob <- bf$prob
      d$predicted_class <- as.integer(d$predicted_prob >= 0.5)
      write.csv(d, file, row.names = FALSE)
    }
  )

  output$download_binary_predictions2 <- downloadHandler(
    filename = function() "binary_bayesian_r0_predictions.csv",
    content = function(file) {
      req(binary_fit())
      bf <- binary_fit()
      d <- bf$data
      d$predicted_prob <- bf$prob
      d$predicted_class <- as.integer(d$predicted_prob >= 0.5)
      write.csv(d, file, row.names = FALSE)
    }
  )

  output$download_binary_cv_predictions_csv <- downloadHandler(
    filename = function() "binary_bayesian_kfold_predictions.csv",
    content = function(file) {
      req(binary_fit())
      bf <- binary_fit()
      if (is.null(bf$cv_predictions) || nrow(bf$cv_predictions) == 0) {
        write.csv(data.frame(Note = "No k-fold predictions available. Set k_folds >= 2 and refit the binary model."), file, row.names = FALSE)
      } else {
        write.csv(bf$cv_predictions, file, row.names = FALSE)
      }
    }
  )

  output$download_binary_roc_csv <- downloadHandler(
    filename = function() "binary_roc_data.csv",
    content = function(file) {
      req(binary_fit())
      write.csv(binary_roc_df(binary_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_predictions <- downloadHandler(
    filename = function() "multiclass_r0_bracket_predictions.csv",
    content = function(file) {
      req(multi_fit())
      mf <- multi_fit()
      d <- mf$data
      d$predicted_class <- mf$pred_class
      write.csv(d, file, row.names = FALSE)
    }
  )
  output$download_multi_predictions2 <- downloadHandler(
    filename = function() "multiclass_r0_bracket_predictions.csv",
    content = function(file) {
      req(multi_fit())
      mf <- multi_fit()
      d <- mf$data
      d$predicted_class <- mf$pred_class
      write.csv(d, file, row.names = FALSE)
    }
  )

  output$download_multi_cv_predictions_csv <- downloadHandler(
    filename = function() "multiclass_r0_kfold_predictions.csv",
    content = function(file) {
      req(multi_fit())
      mf <- multi_fit()
      if (is.null(mf$cv_predictions) || nrow(mf$cv_predictions) == 0) {
        write.csv(data.frame(Note = "No k-fold predictions available. Set k_folds >= 2 and refit the multiclass model."), file, row.names = FALSE)
      } else {
        write.csv(mf$cv_predictions, file, row.names = FALSE)
      }
    }
  )

  output$download_param_ranges_csv <- downloadHandler(
    filename = function() {
      if (identical(input$data_mode, "upload")) {
        "uploaded_parameter_summary.csv"
      } else if (identical(input$data_mode, "custom_ode")) {
        "custom_ode_parameter_ranges.csv"
      } else {
        "ode_parameter_ranges.csv"
      }
    },
    content = function(file) {
      if (identical(input$data_mode, "upload")) {
        req(sim_with_current_threshold())
        params <- current_predictors()
        out <- data.frame(
          parameter = params,
          min = sapply(sim_with_current_threshold()[params], min, na.rm = TRUE),
          max = sapply(sim_with_current_threshold()[params], max, na.rm = TRUE),
          mean = sapply(sim_with_current_threshold()[params], mean, na.rm = TRUE),
          sd = sapply(sim_with_current_threshold()[params], stats::sd, na.rm = TRUE),
          row.names = NULL
        )
        write.csv(out, file, row.names = FALSE)
      } else if (identical(input$data_mode, "custom_ode")) {
        out <- tryCatch(custom_param_ranges(), error = function(e) data.frame(Note = e$message, check.names = FALSE))
        write.csv(out, file, row.names = FALSE)
      } else {
        write.csv(current_param_ranges(), file, row.names = FALSE)
      }
    }
  )

  output$download_param_definitions_csv <- downloadHandler(
    filename = function() "parameter_definitions.csv",
    content = function(file) {
      if (identical(input$data_mode, "upload")) {
        out <- data.frame(
          parameter = current_predictors(),
          definition = "Uploaded numeric parameter column; user-defined meaning.",
          row.names = NULL
        )
      } else {
        out <- if (identical(input$data_mode, "custom_ode")) {
          data.frame(Note = "Parameter definitions are user-defined in custom ODE mode.", check.names = FALSE)
        } else if (identical(input$data_mode, "upload")) {
          data.frame(Note = "Parameter definitions are not available for uploaded datasets.", check.names = FALSE)
        } else {
          parameter_definition_table(input$model_name)
        }
      }
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_sim_summary_csv <- downloadHandler(
    filename = function() "simulation_summary.csv",
    content = function(file) {
      req(sim_with_current_threshold())
      df <- sim_with_current_threshold()
      out <- tibble::tibble(
        Quantity = c("Simulation runs", "Parameters", "Threshold positives", "Threshold negatives",
                     "Mean analytic R0/Re", "Mean numeric growth-rate R0estimate"),
        Value = c(
          nrow(df),
          paste(current_predictors(), collapse = ", "),
          sum(df$threshold_binary == 1, na.rm = TRUE),
          sum(df$threshold_binary == 0, na.rm = TRUE),
          round(mean(df$analytic_R0, na.rm = TRUE), 4),
          round(mean(df$numeric_secondary_R0, na.rm = TRUE), 4)
        )
      )
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_r0_compare_csv <- downloadHandler(
    filename = function() "r0_source_comparison.csv",
    content = function(file) {
      req(sim_with_current_threshold())
      df <- sim_with_current_threshold()
      out <- tibble::tibble(
        Metric = c("Correlation: analytic vs numeric growth-rate R0estimate",
                   "MAE: analytic vs numeric growth-rate R0estimate",
                   "AUC: analytic R0 predicts case increase",
                   "AUC: numeric growth-rate R0estimate predicts case increase"),
        Value = c(
          ifelse(all(is.na(df$analytic_R0)), NA_real_, round(cor(df$analytic_R0, df$numeric_secondary_R0, use = "complete.obs"), 4)),
          ifelse(all(is.na(df$analytic_R0)), NA_real_, round(mean(abs(df$analytic_R0 - df$numeric_secondary_R0), na.rm = TRUE), 4)),
          round(tryCatch(as.numeric(pROC::auc(pROC::roc(df$cases_increase, df$analytic_R0, quiet = TRUE))), error = function(e) NA_real_), 4),
          round(tryCatch(as.numeric(pROC::auc(pROC::roc(df$cases_increase, df$numeric_secondary_R0, quiet = TRUE))), error = function(e) NA_real_), 4)
        )
      )
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_vif_csv <- downloadHandler(
    filename = function() "collinearity_vif.csv",
    content = function(file) {
      req(sim_with_current_threshold())
      write.csv(vif_table_from_data(sim_with_current_threshold(), current_predictors()), file, row.names = FALSE)
    }
  )

  output$download_vif_csv2 <- downloadHandler(
    filename = function() "collinearity_vif.csv",
    content = function(file) {
      req(sim_with_current_threshold())
      write.csv(vif_table_from_data(sim_with_current_threshold(), current_predictors()), file, row.names = FALSE)
    }
  )

  output$download_vif_png <- downloadHandler(
    filename = function() "collinearity_vif.png",
    content = function(file) {
      req(sim_with_current_threshold())
      ggplot2::ggsave(file, plot = plot_collinearity_vif_gg(sim_with_current_threshold(), current_predictors()), width = 8, height = 6, dpi = 300)
    }
  )

  output$download_vif_png2 <- downloadHandler(
    filename = function() "collinearity_vif.png",
    content = function(file) {
      req(sim_with_current_threshold())
      ggplot2::ggsave(file, plot = plot_collinearity_vif_gg(sim_with_current_threshold(), current_predictors()), width = 8, height = 6, dpi = 300)
    }
  )

  output$download_corr_matrix_csv <- downloadHandler(
    filename = function() "parameter_correlation_matrix.csv",
    content = function(file) {
      req(sim_with_current_threshold())
      cor_mat <- correlation_matrix_from_data(sim_with_current_threshold(), current_predictors())
      write.csv(cor_mat, file, row.names = TRUE)
    }
  )

  output$download_corr_matrix_csv2 <- downloadHandler(
    filename = function() "parameter_correlation_matrix.csv",
    content = function(file) {
      req(sim_with_current_threshold())
      cor_mat <- correlation_matrix_from_data(sim_with_current_threshold(), current_predictors())
      write.csv(cor_mat, file, row.names = TRUE)
    }
  )

  output$download_corr_matrix_png <- downloadHandler(
    filename = function() "parameter_correlation_matrix.png",
    content = function(file) {
      req(sim_with_current_threshold())
      ggplot2::ggsave(file, plot = plot_collinearity_corr_gg(sim_with_current_threshold(), current_predictors()), width = 8, height = 7, dpi = 300)
    }
  )

  output$download_corr_matrix_png2 <- downloadHandler(
    filename = function() "parameter_correlation_matrix.png",
    content = function(file) {
      req(sim_with_current_threshold())
      ggplot2::ggsave(file, plot = plot_collinearity_corr_gg(sim_with_current_threshold(), current_predictors()), width = 8, height = 7, dpi = 300)
    }
  )

  output$download_binary_coefficients_csv <- downloadHandler(
    filename = function() "binary_bayesian_posterior_coefficients.csv",
    content = function(file) {
      req(binary_fit())
      bf <- binary_fit()
      write.csv(posterior_coef_table(bf$fit, z_predictors = bf$z_predictors), file, row.names = FALSE)
    }
  )

  output$download_mcmc_diagnostics_csv <- downloadHandler(
    filename = function() "binary_bayesian_mcmc_diagnostics.csv",
    content = function(file) {
      req(binary_fit())
      write.csv(mcmc_diagnostics_table(binary_fit()), file, row.names = FALSE)
    }
  )

  output$download_sensitivity_csv <- downloadHandler(
    filename = function() "parameter_sensitivity_bayesian_logistic.csv",
    content = function(file) {
      req(binary_fit())
      write.csv(sensitivity_from_binary_fit(binary_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_metrics_csv <- downloadHandler(
    filename = function() "multiclass_r0_metrics.csv",
    content = function(file) {
      req(multi_fit())
      write.csv(multiclass_metrics(multi_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_coefficients_csv <- downloadHandler(
    filename = function() "multiclass_successive_binary_coefficients.csv",
    content = function(file) {
      req(multi_fit())
      write.csv(multiclass_coef_table(multi_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_confusion_csv <- downloadHandler(
    filename = function() "multiclass_confusion_matrix.csv",
    content = function(file) {
      req(multi_fit())
      write.csv(multiclass_confusion_df(multi_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_class_distribution_csv <- downloadHandler(
    filename = function() "multiclass_class_distribution.csv",
    content = function(file) {
      req(multi_fit())
      write.csv(multiclass_distribution_df(multi_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_auc_csv <- downloadHandler(
    filename = function() "multiclass_one_vs_rest_auc.csv",
    content = function(file) {
      req(multi_fit())
      write.csv(multiclass_all_auc_df(multi_fit()), file, row.names = FALSE)
    }
  )

  output$download_multi_cv_auc_csv <- downloadHandler(
    filename = function() "multiclass_kfold_one_vs_rest_auc.csv",
    content = function(file) {
      req(multi_fit())
      mf <- multi_fit()
      cv_label <- paste0(mf$k_folds, "-fold cross-validation")
      out <- multiclass_auc_df(mf, validation_label = cv_label)
      if (nrow(out) == 0) out <- data.frame(Note = "No k-fold AUC available. Set k >= 2 and refit multiclass model.")
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_binary_equation_tex <- downloadHandler(
    filename = function() "binary_REEBLRA_equation.tex",
    content = function(file) {
      req(binary_fit())
      writeLines(
        binary_latex_document(
          binary_fit(),
          model_name = ifelse(identical(input$data_mode, "upload"), "Uploaded R0 table", ifelse(identical(input$data_mode, "custom_ode"), "Custom ODE model", input$model_name)),
          data_mode = input$data_mode,
          cutoff = input$r0_cutoff
        ),
        con = file,
        useBytes = TRUE
      )
    }
  )

  output$download_multi_equation_tex <- downloadHandler(
    filename = function() "multiclass_REEBLRA_equations.tex",
    content = function(file) {
      req(multi_fit())
      writeLines(
        multiclass_latex_document(
          multi_fit(),
          model_name = ifelse(identical(input$data_mode, "upload"), "Uploaded R0 table", ifelse(identical(input$data_mode, "custom_ode"), "Custom ODE model", input$model_name)),
          data_mode = input$data_mode
        ),
        con = file,
        useBytes = TRUE
      )
    }
  )

  output$download_binary_roc_png <- downloadHandler(
    filename = function() "binary_roc_curve.png",
    content = function(file) {
      req(binary_fit())
      ggplot2::ggsave(file, plot = plot_binary_roc_gg(binary_fit()), width = 7, height = 5, dpi = 300)
    }
  )

  output$download_binary_probability_png <- downloadHandler(
    filename = function() "binary_predicted_probability_plot.png",
    content = function(file) {
      req(binary_fit())
      ggplot2::ggsave(file, plot = plot_binary_probability_gg(binary_fit()), width = 8, height = 6, dpi = 300)
    }
  )

  output$download_multi_prediction_png <- downloadHandler(
    filename = function() "multiclass_prediction_plot.png",
    content = function(file) {
      req(multi_fit())
      ggplot2::ggsave(file, plot = plot_multiclass_prediction_gg(multi_fit()), width = 8, height = 6, dpi = 300)
    }
  )

  output$download_multi_class_distribution_png <- downloadHandler(
    filename = function() "multiclass_class_distribution.png",
    content = function(file) {
      req(multi_fit())
      ggplot2::ggsave(file, plot = plot_multiclass_distribution_gg(multi_fit()), width = 9, height = 6, dpi = 300)
    }
  )

  output$download_multi_auc_png <- downloadHandler(
    filename = function() "multiclass_one_vs_rest_auc.png",
    content = function(file) {
      req(multi_fit())
      ggplot2::ggsave(file, plot = plot_multiclass_auc_gg(multi_fit(), validation_label = "Apparent / in-sample"), width = 9, height = 6, dpi = 300)
    }
  )

  output$download_multi_cv_auc_png <- downloadHandler(
    filename = function() "multiclass_kfold_one_vs_rest_auc.png",
    content = function(file) {
      req(multi_fit())
      mf <- multi_fit()
      cv_label <- paste0(mf$k_folds, "-fold cross-validation")
      plot_obj <- tryCatch(
        plot_multiclass_auc_gg(mf, validation_label = cv_label),
        error = function(e) {
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, label = "No k-fold AUC available. Set k >= 2 and refit multiclass model.") +
            ggplot2::theme_void()
        }
      )
      ggplot2::ggsave(file, plot = plot_obj, width = 9, height = 6, dpi = 300)
    }
  )

  output$download_sensitivity_png <- downloadHandler(
    filename = function() "parameter_sensitivity_bayesian_logistic.png",
    content = function(file) {
      req(binary_fit())
      ggplot2::ggsave(file, plot = plot_sensitivity_gg(binary_fit()), width = 8, height = 6, dpi = 300)
    }
  )

  output$download_all_outputs_zip <- downloadHandler(
    filename = function() paste0("REEBLRA_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      tmpdir <- file.path(tempdir(), paste0("REEBLRA_outputs_", as.integer(Sys.time())))
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      oldwd <- getwd()
      on.exit(setwd(oldwd), add = TRUE)

      add_csv <- function(obj, filename) {
        tryCatch(write.csv(obj, file.path(tmpdir, filename), row.names = FALSE), error = function(e) NULL)
      }
      add_png <- function(plot_obj, filename, width = 8, height = 6) {
        tryCatch(ggplot2::ggsave(file.path(tmpdir, filename), plot = plot_obj, width = width, height = height, dpi = 300), error = function(e) NULL)
      }
      add_tex <- function(tex_text, filename) {
        tryCatch(writeLines(tex_text, file.path(tmpdir, filename), useBytes = TRUE), error = function(e) NULL)
      }

      if (!is.null(sim_data())) {
        df <- sim_with_current_threshold()
        add_csv(df, "simulation_or_uploaded_data.csv")
        add_csv(tibble::tibble(
          Quantity = c("Simulation runs", "Parameters", "Threshold positives", "Threshold negatives",
                       "Mean analytic R0/Re", "Mean numeric growth-rate R0estimate"),
          Value = c(nrow(df), paste(current_predictors(), collapse = ", "),
                    sum(df$threshold_binary == 1, na.rm = TRUE),
                    sum(df$threshold_binary == 0, na.rm = TRUE),
                    round(mean(df$analytic_R0, na.rm = TRUE), 4),
                    round(mean(df$numeric_secondary_R0, na.rm = TRUE), 4))
        ), "simulation_summary.csv")
        add_csv(vif_table_from_data(df, current_predictors()), "collinearity_vif.csv")
        add_png(plot_collinearity_vif_gg(df, current_predictors()), "collinearity_vif.png", 8, 6)
        add_csv(correlation_matrix_from_data(df, current_predictors()), "parameter_correlation_matrix.csv")
        add_png(plot_collinearity_corr_gg(df, current_predictors()), "parameter_correlation_matrix.png", 8, 7)
      }

      if (identical(input$data_mode, "upload") && !is.null(sim_data())) {
        params <- current_predictors()
        add_csv(data.frame(
          parameter = params,
          min = sapply(sim_with_current_threshold()[params], min, na.rm = TRUE),
          max = sapply(sim_with_current_threshold()[params], max, na.rm = TRUE),
          mean = sapply(sim_with_current_threshold()[params], mean, na.rm = TRUE),
          sd = sapply(sim_with_current_threshold()[params], stats::sd, na.rm = TRUE),
          row.names = NULL
        ), "uploaded_parameter_summary.csv")
        add_csv(data.frame(parameter = params, definition = "Uploaded numeric parameter column; user-defined meaning.", row.names = NULL), "parameter_definitions.csv")
      } else if (identical(input$data_mode, "custom_ode")) {
        add_csv(tryCatch(custom_param_ranges(), error = function(e) data.frame(Note = e$message, check.names = FALSE)), "custom_ode_parameter_ranges.csv")
        add_csv(data.frame(Note = "Parameter definitions are user-defined in custom ODE mode.", check.names = FALSE), "parameter_definitions.csv")
      } else {
        add_csv(current_param_ranges(), "ode_parameter_ranges.csv")
        add_csv(parameter_definition_table(input$model_name), "parameter_definitions.csv")
      }

      if (!is.null(binary_fit())) {
        bf <- binary_fit()
        add_csv(bf$metrics, "binary_metrics.csv")
        add_csv(binary_predictions_df(bf), "binary_predictions.csv")
        add_csv(binary_roc_df(bf), "binary_roc_data.csv")
        if (!is.null(bf$cv_predictions) && nrow(bf$cv_predictions) > 0) add_csv(bf$cv_predictions, "binary_kfold_predictions.csv")
        add_csv(posterior_coef_table(bf$fit, z_predictors = bf$z_predictors), "binary_posterior_coefficients.csv")
        add_csv(mcmc_diagnostics_table(bf), "binary_mcmc_diagnostics.csv")
        add_csv(sensitivity_from_binary_fit(bf), "parameter_sensitivity.csv")
        add_tex(binary_latex_document(
          bf,
          model_name = ifelse(identical(input$data_mode, "upload"), "Uploaded R0 table", ifelse(identical(input$data_mode, "custom_ode"), "Custom ODE model", input$model_name)),
          data_mode = input$data_mode,
          cutoff = input$r0_cutoff
        ), "binary_REEBLRA_equation.tex")
        add_png(plot_binary_roc_gg(bf), "binary_roc_curve.png", 7, 5)
        add_png(plot_binary_probability_gg(bf), "binary_predicted_probability_plot.png", 8, 6)
        add_png(plot_sensitivity_gg(bf), "parameter_sensitivity_plot.png", 8, 6)
      }

      if (!is.null(multi_fit())) {
        mf <- multi_fit()
        add_csv(multiclass_predictions_df(mf), "multiclass_predictions.csv")
        if (!is.null(mf$cv_predictions) && nrow(mf$cv_predictions) > 0) add_csv(mf$cv_predictions, "multiclass_kfold_predictions.csv")
        add_csv(multiclass_metrics(mf), "multiclass_metrics.csv")
        add_csv(multiclass_coef_table(mf), "multiclass_coefficients.csv")
        add_csv(multiclass_confusion_df(mf), "multiclass_confusion_matrix.csv")
        add_csv(multiclass_distribution_df(mf), "multiclass_class_distribution.csv")
        add_csv(multiclass_all_auc_df(mf), "multiclass_one_vs_rest_auc.csv")
        add_csv(multiclass_auc_df(mf, validation_label = paste0(mf$k_folds, "-fold cross-validation")), "multiclass_kfold_one_vs_rest_auc.csv")
        add_tex(multiclass_latex_document(
          mf,
          model_name = ifelse(identical(input$data_mode, "upload"), "Uploaded R0 table", ifelse(identical(input$data_mode, "custom_ode"), "Custom ODE model", input$model_name)),
          data_mode = input$data_mode
        ), "multiclass_REEBLRA_equations.tex")
        add_png(plot_multiclass_prediction_gg(mf), "multiclass_prediction_plot.png", 8, 6)
        add_png(plot_multiclass_distribution_gg(mf), "multiclass_class_distribution.png", 9, 6)
        add_png(plot_multiclass_auc_gg(mf, validation_label = "Apparent / in-sample"), "multiclass_one_vs_rest_auc.png", 9, 6)
        tryCatch(add_png(plot_multiclass_auc_gg(mf, validation_label = paste0(mf$k_folds, "-fold cross-validation")), "multiclass_kfold_one_vs_rest_auc.png", 9, 6), error = function(e) NULL)
      }

      files <- list.files(tmpdir, full.names = TRUE)
      if (length(files) == 0) {
        writeLines("No outputs were available. Run/load data and fit models first.", file.path(tmpdir, "README.txt"))
        files <- list.files(tmpdir, full.names = TRUE)
      }
      setwd(tmpdir)
      utils::zip(zipfile = file, files = basename(files))
    }
  )

}

shinyApp(ui, server)
