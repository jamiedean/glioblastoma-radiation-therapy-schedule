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

###############################################################################################################

pallini2011 <- read.csv('data/pallini2011.csv')

pallini2011$t_primary_surgery_recurrence_surgery <- 
  pallini2011$overall_survival - pallini2011$t_recurrence_surgery_death
pallini2011$cd133_change <- 
  log(pallini2011$cd133_progression/pallini2011$cd133_initial_surgery)

cox_model_overall_survival_primary <- 
  coxph(Surv(overall_survival) ~ age + male + kps + mgmt_methylated + cd133_initial_surgery, data = pallini2011)
summary(cox_model_overall_survival_primary)

pallini_cox_model_overall_survival_primary <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ primary surgery'),
    hr = summary(cox_model_overall_survival_primary)[[8]][, 1],
    lower_ci = summary(cox_model_overall_survival_primary)[[8]][, 3],
    upper_ci = summary(cox_model_overall_survival_primary)[[8]][, 4],
    p_value = paste('p =', as.character(signif(summary(cox_model_overall_survival_primary)[[7]][, 5], 2))))

pallini_forest_plot_overall_survival_primary <- 
  ggplot(pallini_cox_model_overall_survival_primary) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_t_recurrence_surgery_death_primary <- 
  coxph(Surv(t_recurrence_surgery_death) ~ age + male + kps + mgmt_methylated + cd133_initial_surgery,
        data = pallini2011)
summary(cox_model_t_recurrence_surgery_death_primary)

pallini_cox_model_t_recurrence_surgery_death_primary <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ primary surgery'),
    hr = summary(cox_model_t_recurrence_surgery_death_primary)[[8]][, 1],
    lower_ci = summary(cox_model_t_recurrence_surgery_death_primary)[[8]][, 3],
    upper_ci = summary(cox_model_t_recurrence_surgery_death_primary)[[8]][, 4],
    p_value = 
      paste('p =', as.character(signif(summary(cox_model_t_recurrence_surgery_death_primary)[[7]][, 5], 2))))

pallini_forest_plot_t_recurrence_surgery_death_primary <- 
  ggplot(pallini_cox_model_t_recurrence_surgery_death_primary) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_t_primary_surgery_recurrence_surgery_primary <- 
  coxph(Surv(t_primary_surgery_recurrence_surgery) ~ age + male + kps + mgmt_methylated + cd133_initial_surgery,
        data = pallini2011)
summary(cox_model_t_primary_surgery_recurrence_surgery_primary)

pallini_cox_model_t_primary_surgery_recurrence_surgery_primary <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ primary surgery'),
    hr = summary(cox_model_t_primary_surgery_recurrence_surgery_primary)[[8]][, 1],
    lower_ci = summary(cox_model_t_primary_surgery_recurrence_surgery_primary)[[8]][, 3],
    upper_ci = summary(cox_model_t_primary_surgery_recurrence_surgery_primary)[[8]][, 4],
    p_value = 
      paste('p =',
            as.character(signif(summary(cox_model_t_primary_surgery_recurrence_surgery_primary)[[7]][, 5], 2))))

pallini_forest_plot_t_primary_surgery_recurrence_surgery_primary <- 
  ggplot(pallini_cox_model_t_primary_surgery_recurrence_surgery_primary) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_overall_survival_recurrence <- 
  coxph(Surv(overall_survival) ~ age + male + kps + mgmt_methylated + cd133_progression, data = pallini2011)
summary(cox_model_overall_survival_recurrence)

pallini_cox_model_overall_survival_recurrence <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ recurrence surgery'),
    hr = summary(cox_model_overall_survival_recurrence)[[8]][, 1],
    lower_ci = summary(cox_model_overall_survival_recurrence)[[8]][, 3],
    upper_ci = summary(cox_model_overall_survival_recurrence)[[8]][, 4],
    p_value = paste('p =', as.character(signif(summary(cox_model_overall_survival_recurrence)[[7]][, 5], 2))))

pallini_forest_plot_overall_survival_recurrence <- 
  ggplot(pallini_cox_model_overall_survival_recurrence) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_t_recurrence_surgery_death_recurrence <- 
  coxph(Surv(t_recurrence_surgery_death) ~ age + male + kps + mgmt_methylated + cd133_progression,
        data = pallini2011)
summary(cox_model_t_recurrence_surgery_death_recurrence)

pallini_cox_model_t_recurrence_surgery_death_recurrence <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ recurrence surgery'),
    hr = summary(cox_model_t_recurrence_surgery_death_recurrence)[[8]][, 1],
    lower_ci = summary(cox_model_t_recurrence_surgery_death_recurrence)[[8]][, 3],
    upper_ci = summary(cox_model_t_recurrence_surgery_death_recurrence)[[8]][, 4],
    p_value = 
      paste('p =', as.character(signif(summary(cox_model_t_recurrence_surgery_death_recurrence)[[7]][, 5], 2))))

pallini_forest_plot_t_recurrence_surgery_death_recurrence <- 
  ggplot(pallini_cox_model_t_recurrence_surgery_death_recurrence) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_t_primary_surgery_recurrence_surgery_recurrence <- 
  coxph(Surv(t_primary_surgery_recurrence_surgery) ~ age + male + kps + mgmt_methylated + cd133_progression,
        data = pallini2011)
summary(cox_model_t_primary_surgery_recurrence_surgery_recurrence)

pallini_cox_model_t_primary_surgery_recurrence_surgery_recurrence <-
  data.frame(
    covariate = c('Age', 'Male', 'KPS', 'MGMT methylated', 'CD133+ recurrence surgery'),
    hr = summary(cox_model_t_primary_surgery_recurrence_surgery_recurrence)[[8]][, 1],
    lower_ci = summary(cox_model_t_primary_surgery_recurrence_surgery_recurrence)[[8]][, 3],
    upper_ci = summary(cox_model_t_primary_surgery_recurrence_surgery_recurrence)[[8]][, 4],
    p_value = 
      paste('p =', as.character(
        signif(summary(cox_model_t_primary_surgery_recurrence_surgery_recurrence)[[7]][, 5], 2))))

pallini_forest_plot_t_primary_surgery_recurrence_surgery_recurrence <- 
  ggplot(pallini_cox_model_t_primary_surgery_recurrence_surgery_recurrence) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[1]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[1]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

###############################################################################################################

tamura2013 <- read.csv('data/tamura2013.csv')

tamura2013$min_total_dose <- tamura2013$ebrt_dose + tamura2013$min_srs_dose
tamura2013$max_total_dose <- tamura2013$ebrt_dose + tamura2013$max_srs_dose
tamura2013$interval_between_ebrt_and_srs <- tamura2013$interval_from_ebrt - tamura2013$interval_from_srs
tamura2013$interval_from_end_of_treatment <-
  c(6.3, 9.2, 8.9, 6.7, 26.2, 4.6, 7.8, 6.3, 22.7, 3.0, 8.6, 5.7, 11.0, 16.1, 36.9, 11.1, 10.7, 12.7, 6.7, 1.8)
tamura2013$cd133_change <- (tamura2013$cd133_progression)^(1/5) - (tamura2013$cd133_initial_surgery)^(1/5)

cox_model_tmz <- 
  coxph(Surv(interval_from_end_of_treatment) ~ age + male + max_total_dose + tmz + cd133_progression,
        data = tamura2013)
summary(cox_model_tmz)

tamura_cox_model_tmz <-
  data.frame(covariate =
               c('Age', 'Male', 'Maximum total dose', 'Temozolomide', 'CD133+ recurrence surgery'),
             hr = summary(cox_model_tmz)[[8]][, 1],
             lower_ci = summary(cox_model_tmz)[[8]][, 3],
             upper_ci = summary(cox_model_tmz)[[8]][, 4],
             p_value = paste('p =', as.character(signif(summary(cox_model_tmz)[[7]][, 5], 2))))

tamura_forest_plot_tmz <- 
  ggplot(tamura_cox_model_tmz) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[2]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[2]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_cd133_change <- 
  coxph(Surv(interval_from_end_of_treatment) ~ age + male + max_total_dose + chemo + cd133_change,
        data = tamura2013)
summary(cox_model_cd133_change)

tamura_cox_model_cd133_change <-
  data.frame(covariate =
               c('Age', 'Male', 'Maximum total dose', 'Chemotherapy', 'CD133+ change'),
             hr = summary(cox_model_cd133_change)[[8]][, 1],
             lower_ci = summary(cox_model_cd133_change)[[8]][, 3],
             upper_ci = summary(cox_model_cd133_change)[[8]][, 4],
             p_value = paste('p =', as.character(signif(summary(cox_model_cd133_change)[[7]][, 5], 2))))

tamura_forest_plot_cd133_change <- 
  ggplot(tamura_cox_model_cd133_change) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[2]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[2]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

cox_model_cd133_primary <- 
  coxph(Surv(interval_from_end_of_treatment) ~ age + male + max_total_dose + chemo + cd133_initial_surgery,
        data = tamura2013)
summary(cox_model_cd133_primary)

tamura_cox_model_cd133_primary <-
  data.frame(covariate =
               c('Age', 'Male', 'Maximum total dose', 'Chemotherapy', 'CD133+ primary surgery'),
             hr = summary(cox_model_cd133_primary)[[8]][, 1],
             lower_ci = summary(cox_model_cd133_primary)[[8]][, 3],
             upper_ci = summary(cox_model_cd133_primary)[[8]][, 4],
             p_value = paste('p =', as.character(signif(summary(cox_model_cd133_primary)[[7]][, 5], 2))))

tamura_forest_plot_cd133_primary <- 
  ggplot(tamura_cox_model_cd133_primary) +
  geom_point(aes(x = hr, y = covariate), size = 2, color = plot_colors[2]) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, y = covariate),
                 size = 1, height = 0.3, color = plot_colors[2]) +
  geom_text_repel(aes(x = hr, y = covariate, label = p_value), nudge_y = -0.3, size = 5) +
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') +
  xlab('Hazard ratio') + ylab('')

###############################################################################################################

figure <- 
  plot_grid(pallini_forest_plot_overall_survival_primary,
            pallini_forest_plot_t_recurrence_surgery_death_primary,
            pallini_forest_plot_t_primary_surgery_recurrence_surgery_primary,
            pallini_forest_plot_overall_survival_recurrence,
            pallini_forest_plot_t_recurrence_surgery_death_recurrence,
            pallini_forest_plot_t_primary_surgery_recurrence_surgery_recurrence,
            tamura_forest_plot_tmz,
            tamura_forest_plot_cd133_change,
            tamura_forest_plot_cd133_primary,
            ncol = 3,
            labels = 'AUTO')

save_plot('figures/figure_s1.pdf', figure, ncol = 4, nrow = 5, device = cairo_pdf)
