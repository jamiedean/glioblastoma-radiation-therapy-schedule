setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(survival)
library(survminer)
library(swimplot)

source('code/custom_plot_theme.R')

################################################################################################################

perform_survival_analysis <- function(experimental_data, control_data, covariates, metric) {
  
  if (metric == 'OS') {
    covariates_outcome <- append(covariates, c('time_reirradiation_to_death', 'death'))
  } else if (metric == 'PFS') {
    covariates_outcome <- append(covariates, c('time_reirradiation_to_progression', 'progression'))
  } 
  
  rownames(control_data) <- paste0('c', seq(1, nrow(control_data)))
  rownames(experimental_data) <- paste0('r', seq(1, nrow(experimental_data)))
  
  data <- rbind(control_data[, covariates_outcome], experimental_data[, covariates_outcome])
  data$experimental <- c(rep(0, nrow(control_data)), rep(1, nrow(experimental_data)))
  data <- data[complete.cases(data),]
  
  # Univariable Cox regression
  if (metric == 'OS') {
    print(coxph(Surv(time_reirradiation_to_death, death) ~ experimental, data = data))
  } else if (metric == 'PFS') {
    print(coxph(Surv(time_reirradiation_to_progression, progression) ~ experimental, data = data))
  } 
  
  # Multivariable Cox regression
  if (metric == 'OS') {
    cox_model <- coxph(Surv(time_reirradiation_to_death, death) ~ ., data = data)
  } else if (metric == 'PFS') {
    cox_model <- coxph(Surv(time_reirradiation_to_progression, progression) ~ ., data = data)
  } 
  print(summary(cox_model))
  temp <- cox.zph(cox_model)
  plot(temp)
  
  cox_model <-
    data.frame(
      covariate = c('Male', 'Age', 'KPS', 'Progression number', 'Tumor size', 'RT BED', 'Bevacizumab',
                    'Experimental'),
      hr = summary(cox_model)[[8]][, 1],
      lower_ci = summary(cox_model)[[8]][, 3],
      upper_ci = summary(cox_model)[[8]][, 4],
      p_value = paste('p =', as.character(signif(summary(cox_model)[[7]][, 5], 2))))
  
  if (metric == 'PFS') {
    x_axis_label <- 'PFS hazard ratio'
    color <- plot_colors[1]
  } else if (metric == 'OS') {
    x_axis_label <- 'OS hazard ratio'
    color <- plot_colors[2]
  }
  
  forest_plot <- 
    ggplot(cox_model) +
    geom_point(aes(x = hr, y = covariate), size = 2, color = color) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                   size = 1, height = 0.3, color = color) +
    geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
    geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
    xlab(x_axis_label) + ylab('')
  
  return(forest_plot)
}

################################################################################################################

#lrfs_data <- read.csv('data/lrfs.csv')
#lrfs_object <- Surv(time = lrfs_data$months, event = lrfs_data$event)
#lrfs_fit <- survfit(lrfs_object ~ 1, data = lrfs_data, conf.type = 'log-log')
#lrfs_plot <- 
#  ggsurvplot(lrfs_fit, data = lrfs_data, break.time.by = 6, risk.table = FALSE, palette = plot_colors[1],
#             xlab = 'Time (months)', ylab = 'Local recurrence-free survival', legend = 'none',
#             ggtheme = custom_plot_theme)

#drfs_data <- read.csv('data/drfs.csv')
#drfs_object <- Surv(time = drfs_data$months, event = drfs_data$event)
#drfs_fit <- survfit(drfs_object ~ 1, data = drfs_data, conf.type = 'log-log')
#drfs_plot <- 
#  ggsurvplot(drfs_fit, data = drfs_data, break.time.by = 6, risk.table = FALSE, palette = plot_colors[2],
#             xlab = 'Time (months)', ylab = 'Distant recurrence-free survival', legend = 'none',
#             ggtheme = custom_plot_theme)

pfs_data <- read.csv('data/pfs.csv')
pfs_object <- Surv(time = pfs_data$months, event = pfs_data$event)
pfs_fit <- survfit(pfs_object ~ 1, data = pfs_data, conf.type = 'log-log')
pfs_fit
pfs_plot <- 
  ggsurvplot(pfs_fit, data = pfs_data, break.time.by = 6, risk.table = FALSE, palette = plot_colors[1],
             xlab = 'Time (months)', ylab = 'Progression-free survival', legend = 'none',
             ggtheme = custom_plot_theme)

os_data <- read.csv('data/os.csv')
os_object <- Surv(time = os_data$months, event = os_data$event)
os_fit <- survfit(os_object ~ 1, data = os_data, conf.type = 'log-log')
os_fit
os_plot <- 
  ggsurvplot(os_fit, data = os_data, break.time.by = 6, risk.table = FALSE, palette = plot_colors[2],
             xlab = 'Time (months)', ylab = 'Overall survival', legend = 'none',
             ggtheme = custom_plot_theme)

################################################################################################################

external_control <- read.csv('data/external_control_data.csv')
trial <- read.csv('data/trial_data.csv')

# Radiation therapy biologically effective dose
trial$rt_bed <- 47.25

alpha_beta_ratio <- 10
external_control$rt_bed <-
  external_control$rt_dose*(1 + (external_control$rt_dose/external_control$rt_fractions)/alpha_beta_ratio)

# Impute missing KPS data
external_control$kps[is.na(external_control$kps)] <- median(external_control$kps, na.rm = TRUE)

# Bevacizumab
trial$any_bevacizumab <- 
  rowMaxs(as.matrix(trial[, c('concurrent_bevacizumab', 'prior_bevacizumab')]))
external_control$any_bevacizumab <- 
  rowMaxs(as.matrix(external_control[, c('concurrent_bevacizumab', 'prior_bevacizumab')]))

covariates <- c('male', 'age', 'kps', 'progression_number', 'volume', 'rt_bed', 'any_bevacizumab')

p_os_cox_model <- perform_survival_analysis(trial, external_control, covariates, 'OS')
p_pfs_cox_model <- perform_survival_analysis(trial, external_control, covariates, 'PFS')

swimmer_os <- read.csv('data/swimmer_os.csv')
swimmer_resp <- read.csv('data/swimmer_resp.csv')
swimmer_plot <- swimmer_plot(swimmer_os, id = 'casenum', end = 'months', name_fill = 'Best_response',
                             label = 'Best Overall Response') +
  xlab('Patient') + ylab('Time (months)') +
  swimmer_points(swimmer_resp, id = 'casenum', time = 'time', name_shape = 'Response') +
  swimmer_arrows(swimmer_os, id = 'casenum', arrow_start = 'months', cont = 'censor') +
  annotate('text', x = 3, y = 36, label = 'Alive', size = 6) +
  annotate('text', x = 2, y = 36, label = sprintf('\u2192'), size = 12) +
  coord_flip(clip = 'off', ylim = c(0, 40)) +
  scale_fill_manual(values = c(plot_colors[3], plot_colors[1], plot_colors[2])) +
  theme_set(theme_bw(base_size = 16) +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.ticks.length = unit(-0.25, 'cm'),
                    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), size = 16),
                    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), size = 16),
                    plot.title = element_text(hjust = 0.5, size = 16),
                    legend.position = 'right',
                    legend.text = element_text(size = 16)))

marginal_progression_example <- ggdraw() + 
  draw_image('figures/marginal_progression_example_patient_6.png',
             x = 1, width = 1, height = 1, hjust = 1)

distant_progression_example <- ggdraw() + 
  draw_image('figures/distant_progression_example_patient_1.png',
             x = 1, width = 1, height = 1, hjust = 1)

###############################################################################################################

figure <- 
  plot_grid(
    plot_grid(
      plot_grid(pfs_plot$plot, os_plot$plot,
                p_pfs_cox_model, p_os_cox_model,
                labels = 'AUTO'),
      swimmer_plot,
      ncol = 1,
      rel_heights = c(2, 1.5),
      labels = c('', 'E')),
    plot_grid(marginal_progression_example,
              distant_progression_example,
              ncol = 2,
              rel_widths = c(3, 2),
              labels = c('F', 'G')),
    ncol = 1,
    rel_heights = c(3, 1)) +
  theme(plot.background = element_rect(fill = 'white', color = NA))
save_plot('figures/figure_2.png', figure, ncol = 2.5, nrow = 5.5)
save_plot('figures/figure_2.tiff', figure, ncol = 2.5, nrow = 5.5, dpi = 150)
