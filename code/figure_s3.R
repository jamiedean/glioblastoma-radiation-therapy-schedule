setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(MatrixGenerics)
library(survival)
library(survminer)

source('code/custom_plot_theme.R')

###############################################################################################################

plot_covariate_distributions <- function(experimental_data, control_data) {
  
  covariates <- c('male', 'age', 'kps', 'progression_number', 'volume', 'rt_bed', 'any_bevacizumab')
  
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
  
  p_progression_number <- 
    ggplot(data, aes(x = progression_number, y = ..density.., fill = experimental)) +
    geom_histogram(color = 'black', binwidth = 1) +
    facet_wrap(~ experimental, nrow = 2) +
    xlab('Progression number') + ylab('Density') +
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
    plot_grid(p_male, p_age, p_kps, p_progression_number, p_volume, p_rt_bed, p_any_bevacizumab, ncol = 4)
  
  return(p_covariates)
}

perform_time_to_progression_ratio_analysis <- function(experimental_data, control_data, comparison) {
  
  plot_data <- data.frame(
    ttp_rert = c(control_data$time_reirradiation_to_progression,
             experimental_data$time_reirradiation_to_progression),
    ttp_prior = c(control_data$time_previous_therapy_to_progression,
             experimental_data$time_previous_therapy_to_progression),
    ttp_initial = c(control_data$time_initial_diagnosis_to_first_recurrence,
             experimental_data$time_initial_diagnosis_to_first_recurrence),
    dataset = c(rep('External Control',
                    length(control_data$ratio_reirradiation_progression_initial_progression)),
                rep('Trial',
                    length(experimental_data$ratio_reirradiation_progression_initial_progression))))
  plot_data$ttp_rert_ttp_prior_ratio <- plot_data$ttp_rert/plot_data$ttp_prior
  plot_data$ttp_rert_ttp_initial_ratio <- plot_data$ttp_rert/plot_data$ttp_initial
  plot_data$log_ttp_initial <- log(plot_data$ttp_initial)
  plot_data$log_ttp_prior <- log(plot_data$ttp_prior)
  plot_data$log_ttp_rert <- log(plot_data$ttp_rert)
  plot_data$log_ttp_rert_ttp_prior_ratio <- log(plot_data$ttp_rert_ttp_prior_ratio)
  plot_data$log_ttp_rert_ttp_initial_ratio <- log(plot_data$ttp_rert_ttp_initial_ratio)
  
  my_comparisons <- list(c('External Control', 'Trial'))
  
  if (comparison == 'prior_therapy') {
    ratio <- 'log_ttp_rert_ttp_prior_ratio'
    y_label <- 'Log (Re-RT TTP / Prior Therapy TTP)'
  } else if (comparison == 'initial_therapy') {
    ratio <- 'log_ttp_rert_ttp_initial_ratio'
    y_label <- 'Log (Re-RT TTP / 1st Line Therapy TTP)'
  }
  
  p_ttp_ratio <- 
    ggviolin(plot_data, x = 'dataset', y = ratio, fill = 'dataset',
             palette = c(plot_colors[1], plot_colors[2]),
             add = 'boxplot', add.params = list(fill = 'white')) +
    stat_compare_means(comparisons = my_comparisons, label.y = 4, method = 't.test', size = 5) + 
    xlab('Cohort') + ylab(y_label) +
    custom_plot_theme + theme(legend.position = 'none')
  
  return(p_ttp_ratio)
}

perfrom_progression_location_analysis <- function(experimental_data, control_data, covariates, metric) {
  
  covariates_outcome <- append(covariates, 'distal_progression')
  colnames(experimental_data)[which(colnames(trial) == metric)] <- 'distal_progression'
  
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
                 c('Male', 'Age', 'KPS', 'Progression number', 'Tumor size', 'RT BED', 'Bevacizumab',
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

p_covariate_distributions <- plot_covariate_distributions(trial, external_control)

p_ttp_ratio_prior_therapy <- 
  perform_time_to_progression_ratio_analysis(trial, external_control, 'prior_therapy')
p_ttp_ratio_initial_therapy <- 
  perform_time_to_progression_ratio_analysis(trial, external_control, 'initial_therapy')

p_progression_location_ptv_1 <- 
  perfrom_progression_location_analysis(trial, external_control, covariates, 'distal_progression_ptv_1')
p_progression_location_ptv_2 <- 
  perfrom_progression_location_analysis(trial, external_control, covariates, 'distal_progression_ptv_2')

###############################################################################################################

figure <- 
  plot_grid(
    p_covariate_distributions,
    plot_grid(
      p_ttp_ratio_prior_therapy,
      p_ttp_ratio_initial_therapy,
      p_progression_location_ptv_1,
      p_progression_location_ptv_2,
      ncol = 2,
      labels = c('B', 'C', 'D', 'E')),
    ncol = 1,
    labels = c('A', NULL),
    rel_heights = c(2, 2))

save_plot('figures/figure_s3.pdf', figure, ncol = 3, nrow = 6, device = cairo_pdf)

plot(external_control$progression_number, log(external_control$time_initial_diagnosis_to_first_recurrence))
plot(external_control$progression_number, log(external_control$time_reirradiation_to_progression))
plot(external_control$progression_number, log(external_control$ratio_reirradiation_progression_initial_progression))

plot(log(external_control$time_initial_diagnosis_to_first_recurrence),
     log(external_control$time_reirradiation_to_progression))
cor.test(log(external_control$time_initial_diagnosis_to_first_recurrence),
         log(external_control$time_reirradiation_to_progression))

plot(log(external_control$time_previous_therapy_to_progression),
     log(external_control$time_reirradiation_to_progression))
cor.test(log(external_control$time_previous_therapy_to_progression),
         log(external_control$time_reirradiation_to_progression))

summary(lm(log(time_reirradiation_to_progression) ~ 
             log(time_previous_therapy_to_progression), external_control))
summary(lm(log(time_reirradiation_to_progression) ~ 
             progression_number, external_control))
summary(lm(log(time_reirradiation_to_progression) ~ 
             time_initial_diagnosis_to_first_recurrence, external_control))

plot(log(trial$time_initial_diagnosis_to_first_recurrence),
     log(trial$time_reirradiation_to_progression))
cor.test(log(trial$time_initial_diagnosis_to_first_recurrence),
         log(trial$time_reirradiation_to_progression))


mean(external_control$ratio_reirradiation_progression_initial_progression > 1, na.rm = TRUE)
mean(trial$ratio_reirradiation_progression_initial_progression > 1, na.rm = TRUE)

summary(lm(log(time_initial_diagnosis_to_first_recurrence) ~ progression_number, external_control))
summary(lm(log(ratio_reirradiation_progression_initial_progression) ~ progression_number, external_control))
summary(lm(log(time_reirradiation_to_progression) ~ 
             #log(time_initial_diagnosis_to_first_recurrence) +
             #age +
             #male +
             #volume +
             #kps +
             progression_number,
           external_control))

plot_data <- data.frame(
  pfs2 = c(external_control$time_reirradiation_to_progression,
           trial$time_reirradiation_to_progression),
  pfs1 = c(external_control$time_previous_therapy_to_progression,
           trial$time_previous_therapy_to_progression),
  pfs0 = c(external_control$time_initial_diagnosis_to_first_recurrence,
           trial$time_initial_diagnosis_to_first_recurrence),
  dataset = c(rep('External Control', length(external_control$ratio_reirradiation_progression_initial_progression)),
    rep('Trial', length(trial$ratio_reirradiation_progression_initial_progression))))
plot_data$pfs2_pfs1_ratio <- plot_data$pfs2/plot_data$pfs1
plot_data$pfs2_pfs0_ratio <- plot_data$pfs2/plot_data$pfs0
plot_data$log_pfs0 <- log(plot_data$pfs0)
plot_data$log_pfs1 <- log(plot_data$pfs1)
plot_data$log_pfs2 <- log(plot_data$pfs2)
plot_data$log_pfs2_pfs1_ratio <- log(plot_data$pfs2_pfs1_ratio)
plot_data$log_pfs2_pfs0_ratio <- log(plot_data$pfs2_pfs0_ratio)

my_comparisons <- list(c('External Control', 'Trial'))

ggviolin(plot_data, x = 'dataset', y = 'log_pfs2', fill = 'dataset',
         palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons, label.y = 5.5, method = 't.test') + 
  xlab('Cohort') + ylab('Log Reirradiation TTP') +
  custom_plot_theme + theme(legend.position = 'none')

ggviolin(plot_data, x = 'dataset', y = 'log_pfs1', fill = 'dataset',
         palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons, label.y = 5.5, method = 't.test') + 
  xlab('Cohort') + ylab('Log Previous Therapy TTP') +
  custom_plot_theme + theme(legend.position = 'none')

ggviolin(plot_data, x = 'dataset', y = 'log_pfs0', fill = 'dataset',
         palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons, label.y = 5.5, method = 't.test') + 
  xlab('Cohort') + ylab('Log First Line Therapy TTP') +
  custom_plot_theme + theme(legend.position = 'none')

ggviolin(plot_data, x = 'dataset', y = 'log_pfs2_pfs1_ratio', fill = 'dataset',
         palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons, label.y = 3.5, method = 't.test') + 
  xlab('Cohort') + ylab('Log (Re-RT TTP / Prior Therapy TTP)') +
  custom_plot_theme + theme(legend.position = 'none')

ggviolin(plot_data, x = 'dataset', y = 'log_pfs2_pfs0_ratio', fill = 'dataset',
         palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons, label.y = 3, method = 't.test') + 
  xlab('Cohort') + ylab('Log (Re-RT TTP / 1st Line Therapy TTP)') +
  custom_plot_theme + theme(legend.position = 'none')

plot_data$dataset
exp(mean(plot_data[plot_data$dataset == 'trial', 'log_pfs2_pfs1_ratio']))
exp(mean(plot_data[plot_data$dataset == 'external_control', 'log_pfs2_pfs1_ratio'], na.rm = TRUE))




pfs_correlation <- read.csv('data/pfs1_pfs2_correlation.csv')
pfs_correlation$pfs2_pfs1_ratio <- pfs_correlation$pfs2/pfs_correlation$pfs1
pfs_correlation$log_pfs1 <- log(pfs_correlation$pfs1)
pfs_correlation$log_pfs2 <- log(pfs_correlation$pfs2)
pfs_correlation$log_pfs2_pfs1_ratio <- log(pfs_correlation$pfs2/pfs_correlation$pfs1)
pfs_correlation

my_comparisons <- list( c('1', '2'), c('1', '3'), c('1', '4'), c('1', '5'),
                        c('2', '3'), c('2', '4'), c('2', '5'),
                        c('3', '4'), c('3', '5'),
                        c('4', '5'))
ggviolin(pfs_correlation, x = 'progression_number', y = 'log_pfs2_pfs1_ratio', fill = 'progression_number',
         #palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(label.y = 10) +
  xlab('Dataset') + ylab('Log PFS ratio')

ggviolin(pfs_correlation, x = 'progression_number', y = 'log_pfs2', fill = 'progression_number',
         #palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(label.y = 10) +
  xlab('Dataset') + ylab('Log PFS 2')

ggviolin(pfs_correlation, x = 'progression_number', y = 'log_pfs1', fill = 'progression_number',
         #palette = c(plot_colors[1], plot_colors[2]),
         add = 'boxplot', add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(label.y = 10) +
  xlab('Dataset') + ylab('Log PFS 1')

t.test(pfs_correlation[pfs_correlation$progression_number == 1, 'log_pfs2_pfs1_ratio'])
t.test(pfs_correlation[pfs_correlation$progression_number == 2, 'log_pfs2_pfs1_ratio'])
t.test(pfs_correlation[pfs_correlation$progression_number == 3, 'log_pfs2_pfs1_ratio'])
t.test(pfs_correlation[pfs_correlation$progression_number == 4, 'log_pfs2_pfs1_ratio'])
t.test(pfs_correlation[pfs_correlation$progression_number == 5, 'log_pfs2_pfs1_ratio'])

t.test(plot_data[plot_data$dataset == 'external_control', 'log_pfs2_pfs1_ratio'])
t.test(plot_data[plot_data$dataset == 'trial', 'log_pfs2_pfs1_ratio'])
t.test(plot_data[plot_data$dataset == 'external_control', 'log_pfs2_pfs0_ratio'])
t.test(plot_data[plot_data$dataset == 'trial', 'log_pfs2_pfs0_ratio'])

plot(pfs_correlation[pfs_correlation$rt == 1, 'log_pfs1'],
     pfs_correlation[pfs_correlation$rt == 1, 'log_pfs2'])
cor.test(pfs_correlation[pfs_correlation$rt == 1, 'log_pfs1'],
         pfs_correlation[pfs_correlation$rt == 1, 'log_pfs2'])

