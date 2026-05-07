################################################################################
# Bayesian Hierarchical models for TP-Chlorophyll relationships and Seston C:P in
# western Lake Erie based on GLERL- CIGLR water quality data available at NCEI
################################################################################

# SNIPPET 4: VISUALIZING RESULTS IN FIGURES and Tables

library(dplyr)
library(ggplot2)
library(patchwork)
library(brms)
library(tidyverse)
library(posterior)
library(ggrepel)
library(viridis)
library(bayesplot)
library(ggExtra)

work_dir <- "C:/Users/aisab/Desktop/tp_chl"
setwd(work_dir)

site_levels <- c("WE2", "WE4", "WE6", "WE8", "WE9", 
                 "WE12", "WE13", "WE14", "WE15", "WE16")
month_order <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

# ==============================================================================
# JOURNAL DIMENSIONS
# ==============================================================================

col_width  <- 3.54   # 90 mm single column
full_width <- 7.01   # 178 mm full page width
full_sq    <- 7.01   # 178 x 178 mm square

# ==============================================================================
# FIGURE 1: RAW DATA BOXPLOTS (4 rows x 2 columns: Site | Month)
# 90 mm wide, compact height (~150 mm)
# ==============================================================================

create_improved_boxplot <- function(data_input, x_var, y_var, x_label, y_label, 
                                    colors, show_x_label = TRUE, show_y_label = TRUE,
                                    show_x_ticks = TRUE, show_y_ticks = TRUE,
                                    log_scale = FALSE) {
  
  mean_data <- data_input %>%
    filter(!is.na(!!sym(y_var))) %>%
    group_by(!!sym(x_var)) %>%
    summarise(mean_val = mean(!!sym(y_var), na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(data_input, aes(x = !!sym(x_var), y = !!sym(y_var), color = !!sym(x_var))) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "white", linewidth = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.15, size = 0.5) +
    geom_line(data = mean_data, aes(x = !!sym(x_var), y = mean_val, group = 1),
              color = "black", linewidth = 0.6) +
    geom_point(data = mean_data, aes(x = !!sym(x_var), y = mean_val),
               color = "black", size = 1.5, shape = 18) +
    scale_color_manual(values = colors, guide = "none") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.margin = margin(1, 2, 1, 2))
  
  if (!show_x_label) {
    p <- p + labs(x = NULL)
  } else {
    p <- p + labs(x = x_label) + theme(axis.title.x = element_text(size = 8, face = "bold"))
  }
  
  if (!show_y_label) {
    p <- p + labs(y = NULL) + theme(axis.title.y = element_blank())
  } else {
    p <- p + labs(y = y_label) + theme(axis.title.y = element_text(size = 8, face = "bold"))
  }
  
  if (!show_x_ticks) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  }
  
  if (!show_y_ticks) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    p <- p + theme(axis.text.y = element_text(size = 7))
  }
  
  if (log_scale) p <- p + scale_y_log10()
  
  return(p)
}

chla_site <- create_improved_boxplot(RawData_clean, "Site_factor", "Chla", "Site", "Chlorophyll-a (µg/L)", 
                                     site_colors, FALSE, TRUE, FALSE, TRUE, TRUE)
chla_month <- create_improved_boxplot(RawData_clean, "Month_factor", "Chla", "Month", "Chlorophyll-a (µg/L)", 
                                      month_colors, FALSE, FALSE, FALSE, FALSE, TRUE)

tp_site <- create_improved_boxplot(RawData_clean, "Site_factor", "TP", "Site", "Total Phosphorus (µg/L)", 
                                   site_colors, FALSE, TRUE, FALSE, TRUE, TRUE)
tp_month <- create_improved_boxplot(RawData_clean, "Month_factor", "TP", "Month", "Total Phosphorus (µg/L)", 
                                    month_colors, FALSE, FALSE, FALSE, FALSE, TRUE)

poc_site <- create_improved_boxplot(RawData_clean, "Site_factor", "POC", "Site", "POC (mg/L)", 
                                    site_colors, FALSE, TRUE, FALSE, TRUE, TRUE)
poc_month <- create_improved_boxplot(RawData_clean, "Month_factor", "POC", "Month", "POC (mg/L)", 
                                     month_colors, FALSE, FALSE, FALSE, FALSE, TRUE)

cp_site <- create_improved_boxplot(seston_data, "Site_factor", "Seston_C_to_P", "Site", "Seston C:P (molar)", 
                                   site_colors, TRUE, TRUE, TRUE, TRUE, TRUE)
cp_month <- create_improved_boxplot(seston_data, "Month_factor", "Seston_C_to_P", "Month", "Seston C:P (molar)", 
                                    month_colors, TRUE, FALSE, TRUE, FALSE, TRUE)

figure1 <- (chla_site | chla_month) /
  (tp_site | tp_month) /
  (poc_site | poc_month) /
  (cp_site | cp_month)

ggsave("Figure1_Raw_Data_Boxplots.pdf", figure1, 
       width = col_width, height = 5.91, dpi = 600)  # 90 x 150 mm

# ==============================================================================
# FIGURE 2: COMBINED INTERCEPTS, SLOPES & POSTERIOR PREDICTIONS
# 178 x 178 mm full square — larger text for readability
# ==============================================================================

all_intercepts <- c(
  coef_site$intercept, coef_site$intercept_lower, coef_site$intercept_upper,
  coef_month$intercept, coef_month$intercept_lower, coef_month$intercept_upper
)

all_slopes <- c(
  coef_site$slope, coef_site$slope_lower, coef_site$slope_upper,
  coef_month$slope, coef_month$slope_lower, coef_month$slope_upper
)

intercept_limits <- c(min(all_intercepts, na.rm = TRUE) - 0.5, 
                      max(all_intercepts, na.rm = TRUE) + 0.5)
slope_limits <- c(min(all_slopes, na.rm = TRUE) - 0.1, 
                  max(all_slopes, na.rm = TRUE) + 0.1)

# Panel A-B: Joint posterior ellipses
create_ellipse_panel <- function(ellipse_data, point_data, group_name, colors) {
  ggplot() +
    geom_path(data = ellipse_data, aes(x = intercept, y = slope, color = group, group = group),
              linewidth = 0.8, alpha = 0.5) +
    geom_point(data = point_data, aes(x = intercept, y = slope, color = group), size = 2, alpha = 0.8) +
    geom_text_repel(data = point_data, aes(x = intercept, y = slope, label = group, color = group),
                    size = 2, fontface = "bold", 
                    box.padding = 1.0,
                    min.segment.length = 0,
                    max.overlaps = 50, seed = 42) +
    scale_color_manual(values = colors, guide = "none") +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-1, 2)) +
    labs(x = "Intercept", y = "Slope", title = group_name) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5), 
          panel.grid.minor = element_blank())
}

top_site <- create_ellipse_panel(site_posterior$ellipses, site_posterior$points, "Site", site_colors)
top_month <- create_ellipse_panel(month_posterior$ellipses, month_posterior$points, "Month", month_colors)

# Panel C-F: Intercept and Slope panels
create_intercept_panel <- function(coef_data, pop_value, colors) {
  ggplot(coef_data, aes(x = intercept, y = group, color = group)) +
    geom_vline(xintercept = pop_value, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = intercept_lower, xmax = intercept_upper), height = 0.3, linewidth = 0.6) +
    geom_point(size = 2) +
    scale_color_manual(values = colors, guide = "none") +
    scale_y_discrete(limits = rev(levels(coef_data$group))) +
    scale_x_continuous(limits = intercept_limits) +
    labs(x = "Intercept (log scale)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9))
}

create_slope_panel <- function(coef_data, pop_value, colors) {
  ggplot(coef_data, aes(x = slope, y = group, color = group)) +
    geom_vline(xintercept = pop_value, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = slope_lower, xmax = slope_upper), height = 0.3, linewidth = 0.6) +
    geom_point(size = 2) +
    scale_color_manual(values = colors, guide = "none") +
    scale_y_discrete(limits = rev(levels(coef_data$group))) +
    scale_x_continuous(limits = slope_limits) +
    labs(x = "Slope", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          axis.text.y = element_blank())
}

mid_site_int <- create_intercept_panel(coef_site, pop_estimates$Intercept[1], site_colors)
mid_site_slope <- create_slope_panel(coef_site, pop_estimates$Slope[1], site_colors)
mid_month_int <- create_intercept_panel(coef_month, pop_estimates$Intercept[2], month_colors)
mid_month_slope <- create_slope_panel(coef_month, pop_estimates$Slope[2], month_colors)

# Panel G-H: Prediction ribbons restricted to each group's observed TP range
create_prediction_panel_marginal <- function(model, group_var, data_input, colors, 
                                             factor_levels, n_pred = 100) {
  
  pred_list <- list()
  for (grp in factor_levels) {
    grp_data <- data_input %>% filter(!!sym(group_var) == grp)
    grp_tp_range <- range(grp_data$log_TP, na.rm = TRUE)
    log_TP_pred <- seq(grp_tp_range[1], grp_tp_range[2], length.out = n_pred)
    
    newdata <- data.frame(log_TP = log_TP_pred)
    newdata[[group_var]] <- grp
    pred_draws <- posterior_epred(model, newdata = newdata, re_formula = NULL)
    
    pred_list[[as.character(grp)]] <- data.frame(
      log_TP = log_TP_pred, group = grp,
      mean = apply(pred_draws, 2, mean),
      lower_50 = apply(pred_draws, 2, quantile, probs = 0.25),
      upper_50 = apply(pred_draws, 2, quantile, probs = 0.75),
      lower_95 = apply(pred_draws, 2, quantile, probs = 0.025),
      upper_95 = apply(pred_draws, 2, quantile, probs = 0.975)
    )
  }
  
  all_preds <- bind_rows(pred_list)
  all_preds$group <- factor(all_preds$group, levels = factor_levels)
  
  overall_median_x <- median(data_input$log_TP, na.rm = TRUE)
  overall_median_y <- median(data_input$log_Chla, na.rm = TRUE)
  
  x_labels_natural <- c(1, 10, 100, 1000)
  y_labels_natural <- c(1, 10, 100, 1000)
  x_breaks <- log(x_labels_natural)
  y_breaks <- log(y_labels_natural)
  
  p_main <- ggplot(data_input, aes(x = log_TP, y = log_Chla)) +
    geom_point(alpha = 0) +
    geom_vline(xintercept = overall_median_x, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_hline(yintercept = overall_median_y, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_ribbon(data = all_preds, aes(x = log_TP, ymin = lower_95, ymax = upper_95, fill = group), 
                alpha = 0.2, inherit.aes = FALSE) +
    geom_ribbon(data = all_preds, aes(x = log_TP, ymin = lower_50, ymax = upper_50, fill = group), 
                alpha = 0.4, inherit.aes = FALSE) +
    geom_line(data = all_preds, aes(x = log_TP, y = mean, color = group), 
              linewidth = 1.2, inherit.aes = FALSE) +
    scale_color_manual(values = colors, guide = "none") +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels_natural) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels_natural) +
    coord_cartesian(xlim = range(data_input$log_TP, na.rm = TRUE) + c(-0.2, 0.5), 
                    ylim = range(data_input$log_Chla, na.rm = TRUE) + c(-0.5, 1.5)) +
    labs(x = "Total Phosphorus (µg/L)", y = "Chlorophyll-a (µg/L)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  
  p_with_marginals <- ggMarginal(p_main, type = "density", fill = "gray70", color = "gray40", alpha = 0.7)
  
  return(p_with_marginals)
}

bottom_site <- create_prediction_panel_marginal(model_site, "Site_factor", model_data_tpchl, site_colors, site_levels)
bottom_month <- create_prediction_panel_marginal(model_month, "Month_factor", model_data_tpchl, month_colors, month_order)

figure2 <- wrap_plots(
  A = top_site, B = top_month,
  C = mid_site_int, D = mid_site_slope, E = mid_month_int, F = mid_month_slope,
  G = wrap_elements(full = bottom_site), 
  H = wrap_elements(full = bottom_month),
  design = "AABB\nAABB\nCDEF\nCDEF\nGGHH\nGGHH"
)

ggsave("Figure2_Combined.pdf", figure2, 
       width = full_sq, height = full_sq, dpi = 600)  # 178 x 178 mm

# ==============================================================================
# FIGURE 3: SESTON C:P EFFECTS (stacked: Site / Month)
# 90 mm wide
# ==============================================================================

create_cp_panel <- function(effects_df, pop_mean, colors, group_name, levels_order) {
  
  plot_data <- effects_df %>%
    mutate(
      baseline_nat = 10^baseline_intercept,
      lower_nat = 10^baseline_lower,
      upper_nat = 10^baseline_upper
    )
  
  pop_mean_nat <- 10^pop_mean
  
  ggplot(plot_data, aes(x = baseline_nat, y = reorder(group, baseline_nat), color = group)) +
    geom_vline(xintercept = pop_mean_nat, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lower_nat, xmax = upper_nat), height = 0.3, linewidth = 0.4) +
    geom_point(size = 1.8) +
    scale_color_manual(values = colors, guide = "none") +
    scale_y_discrete(limits = rev(levels_order)) +
    labs(x = "Seston C:P (molar)", y = NULL, title = paste("By", group_name)) +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 7),
          plot.margin = margin(2, 4, 2, 2))
}

p_site_cp <- create_cp_panel(cp_site_effects, cp_population_mean, site_colors, "Site", site_levels)
p_month_cp <- create_cp_panel(cp_month_effects, cp_population_mean, month_colors, "Month", month_order)

figure3 <- p_site_cp / p_month_cp

ggsave("Figure3_Seston_CP.pdf", figure3, 
       width = col_width, height = col_width * 1.4, dpi = 600)

# ==============================================================================
# FIGURE 4: C:P VS TP-CHL SLOPE (stacked: Site / Month)
# 90 mm wide
# ==============================================================================

create_combined_panel <- function(combined_df, colors, group_name, cor_result) {
  
  plot_data <- combined_df %>%
    mutate(
      cp_nat = 10^cp_baseline,
      cp_lower_nat = 10^cp_lower,
      cp_upper_nat = 10^cp_upper
    )
  
  cp_xlim_nat <- 10^cp_axis_limits
  
  ggplot(plot_data, aes(x = cp_nat, y = tpchl_slope, color = group)) +
    geom_errorbar(aes(ymin = tpchl_slope_lower, ymax = tpchl_slope_upper), width = 0, linewidth = 0.4, alpha = 0.6) +
    geom_errorbarh(aes(xmin = cp_lower_nat, xmax = cp_upper_nat), height = 0, linewidth = 0.4, alpha = 0.6) +
    geom_point(size = 1.8) +
    geom_text_repel(aes(label = group), size = 2.2, fontface = "bold", max.overlaps = 20, seed = 42) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.5, alpha = 0.3, 
                inherit.aes = FALSE, aes(x = cp_nat, y = tpchl_slope)) +
    scale_color_manual(values = colors, guide = "none") +
    scale_x_continuous(limits = cp_xlim_nat) +
    scale_y_continuous(limits = slope_axis_limits) +
    labs(x = "Seston C:P (molar)", y = "TP-Chl Slope", title = paste("By", group_name)) +
    annotate("text", x = cp_xlim_nat[1] + 5, y = slope_axis_limits[2] - 0.1,
             label = sprintf("r = %.3f\np = %.4f", cor_result$estimate, cor_result$p.value), 
             hjust = 0, vjust = 1, size = 2.5) +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(2, 2, 2, 2))
}

p_combined_site <- create_combined_panel(combined_site, site_colors, "Site", cor_site)
p_combined_month <- create_combined_panel(combined_month, month_colors, "Month", cor_month)

figure4 <- p_combined_site / p_combined_month

ggsave("Figure4_Combined_CP_TPChl.pdf", figure4, 
       width = col_width, height = col_width * 1.6, dpi = 600)

# ==============================================================================
# SUPPLEMENTARY FIGURES
# ==============================================================================

# Figure S1: MCMC Trace Plots
trace_month <- mcmc_trace(as_draws_array(model_month), 
                          pars = c("b_Intercept", "b_log_TP", "sigma"), 
                          facet_args = list(ncol = 1)) +
  labs(title = "Monthly Model") + theme_bw(base_size = 8)

trace_site <- mcmc_trace(as_draws_array(model_site), 
                         pars = c("b_Intercept", "b_log_TP", "sigma"), 
                         facet_args = list(ncol = 1)) +
  labs(title = "Site Model") + theme_bw(base_size = 8)

trace_seston <- mcmc_trace(as_draws_array(seston_model), 
                           pars = c("b_Intercept", "sigma"), 
                           facet_args = list(ncol = 1)) +
  labs(title = "Seston C:P Model") + theme_bw(base_size = 8)

figureS1 <- trace_month / trace_site / trace_seston
ggsave("FigureS1_MCMC_trace.pdf", figureS1, 
       width = col_width, height = col_width * 2.5, dpi = 600)

# Figure S2: Natural Scale TP-Chl (stacked)
create_natural_plot <- function(data, group_var, group_name, colors) {
  ggplot(data, aes(x = TP, y = Chla, color = !!sym(group_var))) +
    geom_point(alpha = 0.5, size = 0.8) + 
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.5) +
    scale_color_manual(values = colors, name = group_name) +
    labs(x = "Total Phosphorus (µg/L)", y = "Chlorophyll-a (µg/L)", 
         title = paste("TP-Chl:", group_name)) +
    theme_bw(base_size = 8) + 
    theme(legend.position = "right",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          legend.key.size = unit(0.25, "cm"))
}

figureS2 <- create_natural_plot(model_data_tpchl, "Month_factor", "Month", month_colors) /
  create_natural_plot(model_data_tpchl, "Site_factor", "Site", site_colors)

ggsave("FigureS2_TP_Chl_natural_scale.pdf", figureS2, 
       width = col_width, height = col_width * 1.8, dpi = 600)

# ==============================================================================
# Figure S3: Population-level posterior ellipses (Month vs Site models)
# ==============================================================================

# Extract population-level joint posteriors for both TP-Chl models
extract_population_ellipse <- function(model, model_name, level = 0.95) {
  draws <- as_draws_df(model)
  pop_int <- draws$b_Intercept
  pop_slope <- draws$b_log_TP
  
  cov_mat <- cov(cbind(pop_int, pop_slope))
  mean_vec <- c(mean(pop_int), mean(pop_slope))
  
  ell <- ellipse(cov_mat, centre = mean_vec, level = level, npoints = 200)
  
  ellipse_df <- data.frame(
    intercept = ell[, 1],
    slope = ell[, 2],
    Model = model_name
  )
  
  point_df <- data.frame(
    intercept = mean_vec[1],
    slope = mean_vec[2],
    Model = model_name
  )
  
  list(ellipse = ellipse_df, point = point_df)
}

pop_month_ell <- extract_population_ellipse(model_month, "Month")
pop_site_ell  <- extract_population_ellipse(model_site, "Site")

pop_ellipses <- bind_rows(pop_month_ell$ellipse, pop_site_ell$ellipse)
pop_points   <- bind_rows(pop_month_ell$point, pop_site_ell$point)

model_colors <- c("Month" = "#E41A1C", "Site" = "#377EB8")

figureS3 <- ggplot() +
  geom_polygon(data = pop_ellipses, aes(x = intercept, y = slope, fill = Model), 
               alpha = 0.25) +
  geom_path(data = pop_ellipses, aes(x = intercept, y = slope, color = Model), 
            linewidth = 0.6) +
  geom_point(data = pop_points, aes(x = intercept, y = slope, color = Model), 
             size = 2.5) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(x = "Population-level intercept (log Chl-a)",
       y = "Population-level slope (log TP effect)") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))

ggsave("FigureS3_Population_Ellipses.pdf", figureS3, 
       width = col_width, height = col_width, dpi = 600)

# ==============================================================================
# Figure S4: Chlorophyll yield (µg Chla / µg TP) at reference TP levels
# ==============================================================================

compute_yield <- function(model, group_var, factor_levels, 
                          TP_values = c(11, 45, 242)) {
  
  log_TP_values <- log(TP_values)
  results_list <- list()
  
  # Population-level predictions
  for (tp_i in seq_along(TP_values)) {
    newdata_pop <- data.frame(log_TP = log_TP_values[tp_i])
    newdata_pop[[group_var]] <- factor_levels[1]
    
    pred_pop <- posterior_epred(model, newdata = newdata_pop, re_formula = NA)
    chla_pop <- exp(pred_pop)
    yield_pop <- chla_pop / TP_values[tp_i]
    
    results_list[[paste0("Pop_", TP_values[tp_i])]] <- data.frame(
      group = "Pop.",
      TP_level = paste0(TP_values[tp_i], " µg/L"),
      TP_value = TP_values[tp_i],
      yield_mean = mean(yield_pop),
      yield_lower = quantile(yield_pop, 0.025),
      yield_upper = quantile(yield_pop, 0.975)
    )
  }
  
  # Group-level predictions
  for (grp in factor_levels) {
    for (tp_i in seq_along(TP_values)) {
      newdata <- data.frame(log_TP = log_TP_values[tp_i])
      newdata[[group_var]] <- grp
      
      pred <- posterior_epred(model, newdata = newdata, re_formula = NULL)
      chla_pred <- exp(pred)
      yield_pred <- chla_pred / TP_values[tp_i]
      
      results_list[[paste0(grp, "_", TP_values[tp_i])]] <- data.frame(
        group = as.character(grp),
        TP_level = paste0(TP_values[tp_i], " µg/L"),
        TP_value = TP_values[tp_i],
        yield_mean = mean(yield_pred),
        yield_lower = quantile(yield_pred, 0.025),
        yield_upper = quantile(yield_pred, 0.975)
      )
    }
  }
  
  bind_rows(results_list)
}

yield_site  <- compute_yield(model_site, "Site_factor", site_levels)
yield_month <- compute_yield(model_month, "Month_factor", month_order)

# Set factor levels for plotting order
yield_site$group <- factor(yield_site$group, 
                           levels = c("Pop.", site_levels))
yield_month$group <- factor(yield_month$group, 
                            levels = c("Pop.", month_order))

# TP level as factor for faceting
yield_site$TP_level <- factor(yield_site$TP_level, 
                              levels = c("11 µg/L", "45 µg/L", "242 µg/L"))
yield_month$TP_level <- factor(yield_month$TP_level, 
                               levels = c("11 µg/L", "45 µg/L", "242 µg/L"))

# Color palettes including Pop. as black
site_colors_with_pop <- c("Pop." = "black", site_colors)
month_colors_with_pop <- c("Pop." = "black", month_colors)

create_yield_panel <- function(yield_data, colors, group_name) {
  ggplot(yield_data, aes(x = group, y = yield_mean, fill = group)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = yield_lower, ymax = yield_upper), 
                  width = 0.3, linewidth = 0.3) +
    facet_wrap(~ TP_level, nrow = 1, scales = "free_x",
               strip.position = "bottom") +
    scale_fill_manual(values = colors, guide = "none") +
    labs(y = expression(paste("yield, ", frac(µg~Chl-italic(a), µg~TP))),
         x = NULL, title = paste("By", group_name)) +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 7),
          strip.placement = "outside",
          plot.margin = margin(2, 2, 2, 2))
}

p_yield_site  <- create_yield_panel(yield_site, site_colors_with_pop, "Site")
p_yield_month <- create_yield_panel(yield_month, month_colors_with_pop, "Month")

figureS4 <- p_yield_site / p_yield_month

ggsave("FigureS4_Yield_BarPlots.pdf", figureS4, 
       width = col_width, height = col_width * 1.8, dpi = 600)