setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(MatrixGenerics)
library(survival)
library(survminer)

source('code/custom_plot_theme.R')

###############################################################################################################

plot_covariate_distributions <- function(experimental_data, control_data) {
  
  covariates <- c('male', 'age', 'kps', 'time_initial_reirradiation', 'volume', 'rt_bed', 'any_bevacizumab')
  
  rownames(control_data) <- paste0('c', seq(1, nrow(control_data)))
  rownames(experimental_data) <- paste0('r', seq(1, nrow(experimental_data)))
  data <- rbind(control_data[, covariates], experimental_data[, covariates])
  data$experimental <- as.factor(c(rep('External control', nrow(control_data)),
                                   rep('Experimental', nrow(experimental_data))))
  data <- data[complete.cases(data),]
  
  p_male <- 
    ggplot(data, aes(x = male, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 1) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('Male') + ylab('Density') +
    scale_x_continuous(breaks = c(0, 1)) +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_age <- 
    ggplot(data, aes(x = age, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 5) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('Age (y)') + ylab('Density') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_kps <- 
    ggplot(data, aes(x = kps, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 10) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('KPS') + ylab('Density') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_time_initial_reirradiation <- 
    ggplot(data, aes(x = time_initial_reirradiation, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 2) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('Time diagnosis - re-RT (months)') + ylab('Density') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_volume <- 
    ggplot(data, aes(x = (volume), y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 200) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab(bquote('Tumor size'~(mm^2))) + ylab('Density') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_rt_bed <- 
    ggplot(data, aes(x = rt_bed, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 5) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('RT BED (Gy)') + ylab('Density') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_any_bevacizumab <- ggplot(data, aes(x = any_bevacizumab, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 1) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('Bevacizumab') + ylab('Density') +
    scale_x_continuous(breaks = c(0, 1)) +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_covariates <- 
    plot_grid(p_male, p_age, p_kps, p_time_initial_reirradiation, p_volume, p_rt_bed, p_any_bevacizumab,
              ncol = 4)
  
  return(p_covariates)
}

perform_survival_analysis <- function(experimental_data, control_data, covariates) {
  
  covariates_outcome <- append(covariates, c('survival_following_reirradiation', 'death'))
  
  rownames(control_data) <- paste0('c', seq(1, nrow(control_data)))
  rownames(experimental_data) <- paste0('r', seq(1, nrow(experimental_data)))
  data <- rbind(control_data[, covariates_outcome], experimental_data[, covariates_outcome])
  data$experimental <- c(rep(0, nrow(control_data)), rep(1, nrow(experimental_data)))
  data <- data[complete.cases(data),]
  
  # Univariable Cox regression
  print(coxph(Surv(survival_following_reirradiation, death) ~ experimental, data = data))
  # Multivariable Cox regression
  cox_model_overall_survival <- coxph(Surv(survival_following_reirradiation, death) ~ ., data = data)
  print(summary(cox_model_overall_survival))
  temp <- cox.zph(cox_model_overall_survival)
  plot(temp)
  
  cox_model_overall_survival <-
    data.frame(
      covariate = c('Male', 'Age', 'KPS', 'Time diagnosis - re-RT', 'Tumor size', 'RT BED', 'Bevacizumab',
                    'Experimental'),
      hr = summary(cox_model_overall_survival)[[8]][, 1],
      lower_ci = summary(cox_model_overall_survival)[[8]][, 3],
      upper_ci = summary(cox_model_overall_survival)[[8]][, 4],
      p_value = paste('p =', as.character(signif(summary(cox_model_overall_survival)[[7]][, 5], 2))))
  
  forest_plot_overall_survival <- 
    ggplot(cox_model_overall_survival) +
    geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                   size = 1, height = 0.3, color = plot_colors[1]) +
    geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
    geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
    xlab('Hazard ratio') + ylab('')
  
  return(forest_plot_overall_survival)
}

perfrom_progression_location_analysis <- function(experimental_data, control_data, covariates) {
  
  covariates_outcome <- append(covariates, 'distal_progression')
  
  rownames(control_data) <- paste0('c', seq(1, nrow(control_data)))
  rownames(experimental_data) <- paste0('r', seq(1, nrow(experimental_data)))
  data <- rbind(control_data[, covariates_outcome], experimental_data[, covariates_outcome])
  data$experimental <- c(rep(0, nrow(control_data)), rep(1, nrow(experimental_data)))
  data <- data[complete.cases(data),]
  
  # Univariable logistic regression
  print(summary(glm(distal_progression ~ experimental, data = data, family = 'binomial')))
  # Multivariable logistic regression
  logistic_regression <- glm(distal_progression ~ ., data = data, family = 'binomial')
  print(summary(logistic_regression))
  
  logistic_regression_odds_ratios <-
    data.frame(covariate =
                 c('Male', 'Age', 'KPS', 'Time diagnosis - re-RT', 'Tumor size', 'RT BED', 'Bevacizumab',
                   'Experimental'),
               or = exp(coef(logistic_regression)[-1]),
               lower_ci = exp(confint(logistic_regression)[-1, 1]),
               upper_ci = exp(confint(logistic_regression)[-1, 2]),
               p_value = paste('p =',
                               as.character(signif(summary(logistic_regression)[['coefficients']][-1, 4], 2))))
  
  forest_plot_distal_progression <- 
    ggplot(logistic_regression_odds_ratios) +
    geom_point(aes(x = or, y = covariate), size = 2, color = plot_colors[2]) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                   size = 1, height = 0.3, color = plot_colors[2]) +
    geom_text_repel(aes(x = or, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
    geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
    xlab('Odds ratio') + ylab('')
  
  return(forest_plot_distal_progression)
}

###############################################################################################################

external_control <- 
  read.csv('data/external_control_data.csv')
trial <- read.csv('data/trial_data.csv')

# Radiation therapy biologically effective dose
trial$rt_bed <- 47.25

alpha_beta_ratio <- 10
external_control$rt_bed <-
  external_control$rt_dose*(1 + (external_control$rt_dose/external_control$rt_fractions)/alpha_beta_ratio)

#trial$time_between_progressions <-
#  trial$time_initial_reirradiation/trial$progression_number
#external_control$time_between_progressions <-
#  external_control$time_initial_reirradiation/external_control$progression_number

# Impute missing KPS data
external_control$kps[is.na(external_control$kps)] <- median(external_control$kps, na.rm = TRUE)

# Bevacizumab
trial$any_bevacizumab <- 
  rowMaxs(as.matrix(trial[, c('concurrent_bevacizumab', 'prior_bevacizumab')]))
external_control$any_bevacizumab <- 
  rowMaxs(as.matrix(external_control[, c('concurrent_bevacizumab', 'prior_bevacizumab')]))

covariates <- c('male', 'age', 'kps', 'time_initial_reirradiation', 'volume', 'rt_bed', 'any_bevacizumab')

p_covariate_distributions <- plot_covariate_distributions(trial, external_control)

p_survival <- perform_survival_analysis(trial, external_control, covariates)
p_progression_location <- perfrom_progression_location_analysis(trial, external_control, covariates)

###############################################################################################################

figure <- 
  plot_grid(
    p_covariate_distributions,
    plot_grid(
      p_survival,
      p_progression_location,
      ncol = 2,
      labels = c('B', 'C')),
    ncol = 1,
    labels = c('A', NULL),
    rel_heights = c(2, 1))

save_plot('figures/figure_s3.pdf', figure, ncol = 3, nrow = 5, device = cairo_pdf)
