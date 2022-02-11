setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(ggalt)
library(ggplot2)
library(ggrepel)
library(grid)
library(gtable)
library(reshape2)
library(survival)
library(survminer)
library(tidyverse)

source('code/custom_plot_theme.R')

################################################################################################################

data <- read.csv('data/interfraction_interval.csv')
data$study_number <- rownames(data)

data <- arrange(data, data$interval_control)

shaded_region <- data.frame(xmin = 1.75, xmax = 4.75, ymin = -Inf, ymax = Inf)

p_interfraction_interval <- ggplot() +
  geom_rect(data = shaded_region, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            alpha = 0.3, fill = plot_colors[4]) +
  geom_dumbbell(data = data[data$experimental_superior == 'superior',],
                aes(x = interval_experimental, xend = interval_control, y = reference_abbreviated,
                    col = species), size_x = 7, size_xend = 1) +
  geom_dumbbell(data = data[data$experimental_superior == 'equivalent',],
                aes(x = interval_experimental, xend = interval_control, y = reference_abbreviated,
                    col = species), size_x = 4, size_xend = 4) +
  xlab('Interfraction interval (h)') + ylab('Study') +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  scale_color_manual(values = c(plot_colors[1], plot_colors[2], plot_colors[3])) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  custom_plot_theme

################################################################################################################

pallini2011 <- read.csv('data/pallini2011.csv')

pallini2011$t_primary_surgery_recurrence_surgery <- 
  pallini2011$overall_survival - pallini2011$t_recurrence_surgery_death
pallini2011$cd133_change <- 
  log(pallini2011$cd133_progression/pallini2011$cd133_initial_surgery)

cox_model_overall_survival <- 
  coxph(Surv(overall_survival) ~ age + male + kps + mgmt_methylated + cd133_change, data = pallini2011)
summary(cox_model_overall_survival)

pallini_cox_model_overall_survival <-
  data.frame(covariate =
               c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ change'),
             hr = summary(cox_model_overall_survival)[[8]][, 1],
             lower_ci = summary(cox_model_overall_survival)[[8]][, 3],
             upper_ci = summary(cox_model_overall_survival)[[8]][, 4],
             p_value = paste('p =', as.character(signif(summary(cox_model_overall_survival)[[7]][, 5], 2))))

pallini_forest_plot_overall_survival <- 
  ggplot(pallini_cox_model_overall_survival) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')
pallini_forest_plot_overall_survival

cox_model_t_recurrence_surgery_death <- 
  coxph(Surv(t_recurrence_surgery_death) ~ age + male + kps + mgmt_methylated + cd133_change, data = pallini2011)
summary(cox_model_t_recurrence_surgery_death)

pallini_cox_model_t_recurrence_surgery_death <-
  data.frame(covariate =
               c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ change'),
             hr = summary(cox_model_t_recurrence_surgery_death)[[8]][, 1],
             lower_ci = summary(cox_model_t_recurrence_surgery_death)[[8]][, 3],
             upper_ci = summary(cox_model_t_recurrence_surgery_death)[[8]][, 4],
             p_value = 
               paste('p =', as.character(signif(summary(cox_model_t_recurrence_surgery_death)[[7]][, 5], 2))))

pallini_forest_plot_t_recurrence_surgery_death <- 
  ggplot(pallini_cox_model_t_recurrence_surgery_death) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')
pallini_forest_plot_t_recurrence_surgery_death

cox_model_t_primary_surgery_recurrence_surgery <- 
  coxph(Surv(t_primary_surgery_recurrence_surgery) ~ age + male + kps + mgmt_methylated + cd133_change,
        data = pallini2011)
summary(cox_model_t_primary_surgery_recurrence_surgery)

pallini_cox_model_t_primary_surgery_recurrence_surgery <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ change'),
    hr = summary(cox_model_t_primary_surgery_recurrence_surgery)[[8]][, 1],
    lower_ci = summary(cox_model_t_primary_surgery_recurrence_surgery)[[8]][, 3],
    upper_ci = summary(cox_model_t_primary_surgery_recurrence_surgery)[[8]][, 4],
    p_value = 
      paste('p =', as.character(signif(summary(cox_model_t_primary_surgery_recurrence_surgery)[[7]][, 5], 2))))

pallini_forest_plot_t_primary_surgery_recurrence_surgery <- 
  ggplot(pallini_cox_model_t_primary_surgery_recurrence_surgery) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')
pallini_forest_plot_t_primary_surgery_recurrence_surgery

################################################################################################################

tamura2013 <- read.csv('data/tamura2013.csv')

tamura2013$min_total_dose <- tamura2013$ebrt_dose + tamura2013$min_srs_dose
tamura2013$max_total_dose <- tamura2013$ebrt_dose + tamura2013$max_srs_dose
tamura2013$interval_between_ebrt_and_srs <- tamura2013$interval_from_ebrt - tamura2013$interval_from_srs
tamura2013$interval_from_end_of_treatment <-
  c(6.3, 9.2, 8.9, 6.7, 26.2, 4.6, 7.8, 6.3, 22.7, 3.0, 8.6, 5.7, 11.0, 16.1, 36.9, 11.1, 10.7, 12.7, 6.7, 1.8)
tamura2013$cd133_change <- (tamura2013$cd133_progression)^(1/5) - (tamura2013$cd133_initial_surgery)^(1/5)

cox_model <- 
  coxph(Surv(interval_from_end_of_treatment) ~ age + male + max_total_dose + chemo + cd133_progression,
        data = tamura2013)
summary(cox_model)

tamura_cox_model <-
  data.frame(covariate =
               c('Age', 'Male', 'Maximum total dose', 'Chemotherapy', 'CD133+ fraction at progression'),
             hr = summary(cox_model)[[8]][, 1],
             lower_ci = summary(cox_model)[[8]][, 3],
             upper_ci = summary(cox_model)[[8]][, 4],
             p_value = paste('p =', as.character(signif(summary(cox_model)[[7]][, 5], 2))))

tamura_forest_plot <- 
  ggplot(tamura_cox_model) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[2]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[2]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

################################################################################################################

figure <- 
  plot_grid(
    p_interfraction_interval,
    plot_grid(
      pallini_forest_plot_overall_survival, pallini_forest_plot_t_recurrence_surgery_death,
      pallini_forest_plot_t_primary_surgery_recurrence_surgery, tamura_forest_plot,
      labels = c('B', 'C', 'D', 'E'),
      ncol = 2),
    nrow = 2,
    rel_heights = c(0.66, 1))

save_plot('figures/figure_3.pdf', figure, ncol = 3, nrow = 4, device = cairo_pdf)
