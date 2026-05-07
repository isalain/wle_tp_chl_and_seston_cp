################################################################################
# Bayesian Hierarchical models for TP-Chlorophyll relationships and Seston C:P in
# western Lake Erie based on GLERL- CIGLR water quality data available at NCEI
################################################################################


# SNIPPET 3: MODEL PARAMETER EXTRACTION AND PREDICTIONS

library(brms)
library(dplyr)
library(tidyverse)
library(posterior)
library(ellipse)

work_dir <- "C:/Users/aisab/Desktop/tp_chl"
setwd(work_dir)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

site_levels <- c("WE2", "WE4", "WE6", "WE8", "WE9", 
                 "WE12", "WE13", "WE14", "WE15", "WE16")
month_order <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

extract_coef_simple <- function(model, group_var) {
  coef_summary <- coef(model, summary = TRUE, robust = TRUE)[[group_var]]
  
  intercepts <- data.frame(
    group = rownames(coef_summary[, , "Intercept"]),
    intercept = coef_summary[, "Estimate", "Intercept"],
    intercept_lower = coef_summary[, "Q2.5", "Intercept"],
    intercept_upper = coef_summary[, "Q97.5", "Intercept"]
  )
  
  slopes <- data.frame(
    group = rownames(coef_summary[, , "log_TP"]),
    slope = coef_summary[, "Estimate", "log_TP"],
    slope_lower = coef_summary[, "Q2.5", "log_TP"],
    slope_upper = coef_summary[, "Q97.5", "log_TP"]
  )
  
  result <- merge(intercepts, slopes, by = "group")
  rownames(result) <- NULL
  return(result)
}

extract_joint_posterior_ellipses <- function(model, group_var, factor_levels = NULL) {
  draws <- as_draws_df(model)
  
  intercept_cols <- grep(paste0("r_", group_var, "\\[.+,Intercept\\]"), names(draws), value = TRUE)
  slope_cols <- grep(paste0("r_", group_var, "\\[.+,log_TP\\]"), names(draws), value = TRUE)
  
  pop_intercept <- draws$b_Intercept
  pop_slope <- draws$b_log_TP
  
  group_names <- gsub(paste0("r_", group_var, "\\[(.+),Intercept\\]"), "\\1", intercept_cols)
  
  ellipse_list <- list()
  
  for (i in seq_along(group_names)) {
    grp <- group_names[i]
    
    group_intercept_dev <- draws[[intercept_cols[i]]]
    group_slope_dev <- draws[[slope_cols[i]]]
    
    group_intercept_abs <- pop_intercept + group_intercept_dev
    group_slope_abs <- pop_slope + group_slope_dev
    
    cov_mat <- cov(cbind(group_intercept_abs, group_slope_abs))
    mean_vec <- c(median(group_intercept_abs), median(group_slope_abs))
    
    ell <- ellipse(cov_mat, centre = mean_vec, level = 0.95, npoints = 100)
    
    ellipse_df <- data.frame(
      intercept = ell[,1],
      slope = ell[,2],
      group = grp
    )
    
    ellipse_list[[grp]] <- ellipse_df
  }
  
  all_ellipses <- bind_rows(ellipse_list)
  
  point_estimates <- data.frame(
    group = group_names,
    intercept = sapply(group_names, function(g) {
      median(pop_intercept + draws[[paste0("r_", group_var, "[", g, ",Intercept]")]])
    }),
    slope = sapply(group_names, function(g) {
      median(pop_slope + draws[[paste0("r_", group_var, "[", g, ",log_TP]")]])
    })
  )
  
  if (!is.null(factor_levels)) {
    all_ellipses$group <- factor(all_ellipses$group, levels = factor_levels)
    point_estimates$group <- factor(point_estimates$group, levels = factor_levels)
  }
  
  return(list(ellipses = all_ellipses, points = point_estimates))
}

extract_hierarchical_estimates <- function(model, factor_name, levels_order) {
  ranef_names <- names(ranef(model))
  
  if (factor_name %in% ranef_names) {
    group_coef <- coef(model, robust=TRUE)[[factor_name]]
    
    estimates_df <- data.frame(
      group = rownames(group_coef[,,"Intercept"]),
      estimate = group_coef[,"Estimate","Intercept"],
      lower = group_coef[,"Q2.5","Intercept"],
      upper = group_coef[,"Q97.5","Intercept"]
    )
    
  } else {
    fixef_summary <- fixef(model, robust=TRUE)
    
    estimates_df <- data.frame(
      group = levels_order,
      estimate = NA,
      lower = NA,
      upper = NA
    )
    
    estimates_df$estimate[1] <- fixef_summary["Intercept", "Estimate"]
    estimates_df$lower[1] <- fixef_summary["Intercept", "Q2.5"]
    estimates_df$upper[1] <- fixef_summary["Intercept", "Q97.5"]
    
    for (i in 2:length(levels_order)) {
      level <- levels_order[i]
      param_name <- paste0(factor_name, level)
      
      if (param_name %in% rownames(fixef_summary)) {
        estimates_df$estimate[i] <- fixef_summary["Intercept", "Estimate"] + 
          fixef_summary[param_name, "Estimate"]
        
        int_se <- (fixef_summary["Intercept", "Q97.5"] - 
                     fixef_summary["Intercept", "Q2.5"]) / (2 * 1.96)
        contrast_se <- (fixef_summary[param_name, "Q97.5"] - 
                          fixef_summary[param_name, "Q2.5"]) / (2 * 1.96)
        total_se <- sqrt(int_se^2 + contrast_se^2)
        
        estimates_df$lower[i] <- estimates_df$estimate[i] - 1.96 * total_se
        estimates_df$upper[i] <- estimates_df$estimate[i] + 1.96 * total_se
      }
    }
  }
  
  estimates_df$group <- factor(estimates_df$group, levels = levels_order)
  return(estimates_df)
}

get_population_intercept <- function(model) {
  fixef_summary <- fixef(model, robust=TRUE)
  return(fixef_summary["Intercept", "Estimate"])
}

get_population_slope <- function(model) {
  fixef_summary <- fixef(model)
  return(fixef_summary["log_TP", "Estimate"])
}

# ==============================================================================
# EXTRACT TP-CHLOROPHYLL MODEL COEFFICIENTS
# ==============================================================================

coef_site <- extract_coef_simple(model_site, "Site_factor")
coef_site$group <- factor(coef_site$group, levels = site_levels)

coef_month <- extract_coef_simple(model_month, "Month_factor")
coef_month$group <- factor(coef_month$group, levels = month_order)

# ==============================================================================
# POPULATION-LEVEL ESTIMATES
# ==============================================================================

pop_estimates <- data.frame(
  Model = c("Site", "Month"),
  Intercept = c(
    get_population_intercept(model_site),
    get_population_intercept(model_month)
  ),
  Slope = c(
    get_population_slope(model_site),
    get_population_slope(model_month)
  )
)

# ==============================================================================
# EXTRACT JOINT POSTERIOR ELLIPSES
# ==============================================================================

site_posterior <- extract_joint_posterior_ellipses(model_site, "Site_factor", site_levels)
month_posterior <- extract_joint_posterior_ellipses(model_month, "Month_factor", month_order)

# ==============================================================================
# EXTRACT SESTON C:P MODEL ESTIMATES
# ==============================================================================

cp_site_effects <- extract_hierarchical_estimates(seston_model, "Site_factor", site_levels)
names(cp_site_effects) <- c("group", "baseline_intercept", "baseline_lower", "baseline_upper")

cp_month_effects <- extract_hierarchical_estimates(seston_model, "Month_factor", month_order)
names(cp_month_effects) <- c("group", "baseline_intercept", "baseline_lower", "baseline_upper")

cp_population_mean <- get_population_intercept(seston_model)

# ==============================================================================
# COMBINE TP-CHL SLOPES WITH C:P ESTIMATES
# ==============================================================================

combined_site <- data.frame(
  group = cp_site_effects$group,
  cp_baseline = cp_site_effects$baseline_intercept,
  cp_lower = cp_site_effects$baseline_lower,
  cp_upper = cp_site_effects$baseline_upper,
  tpchl_slope = coef_site$slope[match(cp_site_effects$group, coef_site$group)],
  tpchl_slope_lower = coef_site$slope_lower[match(cp_site_effects$group, coef_site$group)],
  tpchl_slope_upper = coef_site$slope_upper[match(cp_site_effects$group, coef_site$group)]
)
combined_site$group <- factor(combined_site$group, levels = site_levels)

combined_month <- data.frame(
  group = cp_month_effects$group,
  cp_baseline = cp_month_effects$baseline_intercept,
  cp_lower = cp_month_effects$baseline_lower,
  cp_upper = cp_month_effects$baseline_upper,
  tpchl_slope = coef_month$slope[match(cp_month_effects$group, coef_month$group)],
  tpchl_slope_lower = coef_month$slope_lower[match(cp_month_effects$group, coef_month$group)],
  tpchl_slope_upper = coef_month$slope_upper[match(cp_month_effects$group, coef_month$group)]
)
combined_month$group <- factor(combined_month$group, levels = month_order)

# Correlation tests
cor_site <- cor.test(combined_site$cp_baseline, combined_site$tpchl_slope)
cor_month <- cor.test(combined_month$cp_baseline, combined_month$tpchl_slope)

# ==============================================================================
# CALCULATE AXIS LIMITS
# ==============================================================================

all_intercepts <- c(
  site_posterior$ellipses$intercept,
  month_posterior$ellipses$intercept
)
all_slopes <- c(
  site_posterior$ellipses$slope,
  month_posterior$ellipses$slope
)

x_lim_fig2 <- range(all_intercepts, na.rm = TRUE) + c(-0.5, 0.5)
y_lim_fig2 <- range(all_slopes, na.rm = TRUE) + c(-0.1, 0.1)

log_TP_range <- range(model_data_tpchl$log_TP, na.rm = TRUE)
log_Chla_range <- range(model_data_tpchl$log_Chla, na.rm = TRUE)
x_lim_pred <- log_TP_range + c(-0.2, 0.5)
y_lim_pred <- log_Chla_range + c(-0.5, 1.5)

all_cp <- c(
  cp_site_effects$baseline_intercept, cp_site_effects$baseline_lower, cp_site_effects$baseline_upper,
  cp_month_effects$baseline_intercept, cp_month_effects$baseline_lower, cp_month_effects$baseline_upper
)
cp_axis_limits <- range(all_cp, na.rm = TRUE) + c(-0.05, 0.05)

all_tpchl_slopes <- c(
  combined_site$tpchl_slope, combined_site$tpchl_slope_lower, combined_site$tpchl_slope_upper,
  combined_month$tpchl_slope, combined_month$tpchl_slope_lower, combined_month$tpchl_slope_upper
)
slope_axis_limits <- range(all_tpchl_slopes, na.rm = TRUE) + c(-0.1, 0.1)

# ==============================================================================
# PREDICTION TABLES
# ==============================================================================

create_prediction_table <- function(model, group_var, factor_levels, 
                                    TP_values = c(10, 50, 100)) {
  
  log_TP_values <- log(TP_values)
  
  results_list <- list()
  
  for (grp in factor_levels) {
    newdata <- data.frame(log_TP = log_TP_values)
    newdata[[group_var]] <- grp
    
    pred_draws <- posterior_epred(model, newdata = newdata, re_formula = NULL)
    
    pred_natural <- exp(pred_draws)
    
    result_df <- data.frame(
      Group = grp,
      TP_ugL = TP_values,
      Chla_mean = apply(pred_natural, 2, mean),
      Chla_median = apply(pred_natural, 2, median),
      Chla_lower95 = apply(pred_natural, 2, quantile, probs = 0.025),
      Chla_upper95 = apply(pred_natural, 2, quantile, probs = 0.975)
    )
    
    results_list[[as.character(grp)]] <- result_df
  }
  
  all_results <- bind_rows(results_list)
  return(all_results)
}

table_site <- create_prediction_table(model_site, "Site_factor", site_levels)
table_month <- create_prediction_table(model_month, "Month_factor", month_order)

write.csv(table_site, "Table_Predictions_by_Site.csv", row.names = FALSE)
write.csv(table_month, "Table_Predictions_by_Month.csv", row.names = FALSE)