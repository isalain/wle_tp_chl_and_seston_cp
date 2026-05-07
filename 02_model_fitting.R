################################################################################
# Bayesian Hierarchical models for TP-Chlorophyll relationships and Seston C:P in
# western Lake Erie based on GLERL- CIGLR water quality data available at NCEI
################################################################################

# SNIPPET 2: BAYESIAN MODEL FITTING

library(brms)
library(dplyr)

work_dir <- "C:/Users/aisab/Desktop/tp_chl"
setwd(work_dir)

options(mc.cores = parallel::detectCores())

# ==============================================================================
# MODEL 1: TP-CHLOROPHYLL BY MONTH
# ==============================================================================

if (file.exists("model_month.rds")) {
  model_month <- readRDS("model_month.rds")
} else {
  model_month <- brm(
    formula = log_Chla ~ log_TP + (1 + log_TP | Month_factor),
    data = model_data_tpchl,
    family = gaussian(),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    seed = 12345,
    file = "model_month"
  )
  saveRDS(model_month, "model_month.rds")
}

summary(model_month)

# ==============================================================================
# MODEL 2: TP-CHLOROPHYLL BY SITE 
# ==============================================================================

if (file.exists("model_site.rds")) {
  model_site <- readRDS("model_site.rds")
} else {
  model_site <- brm(
    formula = log_Chla ~ log_TP + (1 + log_TP | Site_factor),
    data = model_data_tpchl,
    family = gaussian(),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    control = list(adapt_delta = 0.97),
    seed = 12345,
    file = "model_site"
  )
  saveRDS(model_site, "model_site.rds")
}

summary(model_site)

# ==============================================================================
# MODEL 3: SESTON C:P HIERARCHICAL MODEL (MONTH + SITE)
# ==============================================================================

if (file.exists("seston_model_month_site.rds")) {
  seston_model <- readRDS("seston_model_month_site.rds")
} else {
  seston_model <- brm(
    formula = log_CP ~ (1 | Month_factor) + (1 | Site_factor),
    data = seston_data,
    family = gaussian(),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    control = list(adapt_delta = 0.99),
    seed = 12345,
    file = "seston_model_month_site"
  )
  saveRDS(seston_model, "seston_model_month_site.rds")
}

summary(seston_model)

# ==============================================================================
# MODEL DIAGNOSTICS
# ==============================================================================

get_diagnostics <- function(model) {
  model_summary <- summary(model)
  
  rhat_fixed <- model_summary$fixed[, "Rhat"]
  rhat_random <- model_summary$random
  
  rhat_vals <- rhat_fixed
  if (length(rhat_random) > 0) {
    for (re in names(rhat_random)) {
      rhat_vals <- c(rhat_vals, rhat_random[[re]][, "Rhat"])
    }
  }
  
  ess_fixed <- model_summary$fixed[, "Bulk_ESS"]
  ess_vals <- ess_fixed
  if (length(rhat_random) > 0) {
    for (re in names(rhat_random)) {
      ess_vals <- c(ess_vals, rhat_random[[re]][, "Bulk_ESS"])
    }
  }
  
  r2 <- bayes_R2(model)
  
  list(
    rhat_range = c(min(rhat_vals, na.rm = TRUE), max(rhat_vals, na.rm = TRUE)),
    ess_range = c(min(ess_vals, na.rm = TRUE), max(ess_vals, na.rm = TRUE)),
    bayes_r2 = r2[1],
    bayes_r2_sd = r2[2]
  )
}

diag_month <- get_diagnostics(model_month)
diag_site  <- get_diagnostics(model_site)
diag_seston <- get_diagnostics(seston_model)