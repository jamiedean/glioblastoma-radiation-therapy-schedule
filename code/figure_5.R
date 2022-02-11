setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(survival)
library(survminer)
library(swimplot)

source('code/custom_plot_theme.R')

################################################################################################################

lrfs_data <- read.csv('data/lrfs.csv')
lrfs_object <- Surv(time = lrfs_data$months, event = lrfs_data$event)
lrfs_fit <- survfit(lrfs_object ~ 1, data = lrfs_data, conf.type = 'log-log')
lrfs_plot <- 
  ggsurvplot(lrfs_fit, data = lrfs_data, break.time.by = 2, risk.table = FALSE, palette = plot_colors[1],
             xlab = 'Time (months)', ylab = 'Local recurrence-free survival', legend = 'none',
             ggtheme = custom_plot_theme)

drfs_data <- read.csv('data/drfs.csv')
drfs_object <- Surv(time = drfs_data$months, event = drfs_data$event)
drfs_fit <- survfit(drfs_object ~ 1, data = drfs_data, conf.type = 'log-log')
drfs_plot <- 
  ggsurvplot(drfs_fit, data = drfs_data, break.time.by = 2, risk.table = FALSE, palette = plot_colors[2],
             xlab = 'Time (months)', ylab = 'Distant recurrence-free survival', legend = 'none',
             ggtheme = custom_plot_theme)

pfs_data <- read.csv('data/pfs.csv')
pfs_object <- Surv(time = pfs_data$months, event = pfs_data$event)
pfs_fit <- survfit(pfs_object ~ 1, data = pfs_data, conf.type = 'log-log')
pfs_plot <- 
  ggsurvplot(pfs_fit, data = pfs_data, break.time.by = 2, risk.table = FALSE, palette = plot_colors[3],
             xlab = 'Time (months)', ylab = 'Progression-free survival', legend = 'none',
             ggtheme = custom_plot_theme)

os_data <- read.csv('data/os.csv')
os_object <- Surv(time = os_data$months, event = os_data$event)
os_fit <- survfit(os_object ~ 1, data = os_data, conf.type = 'log-log')
os_plot <- 
  ggsurvplot(os_fit, data = os_data, break.time.by = 2, risk.table = FALSE, palette = plot_colors[4],
             xlab = 'Time (months)', ylab = 'Overall survival', legend = 'none',
             ggtheme = custom_plot_theme)

swimmer_os <- read.csv('data/swimmer_os.csv')
swimmer_resp <- read.csv('data/swimmer_resp.csv')
swimmer_plot <- swimmer_plot(swimmer_os, id = 'casenum', end = 'months', name_fill = 'Best_response',
                             label = 'Best Overall Response') +
  xlab('Patient') + ylab('Time (months)') +
  swimmer_points(swimmer_resp, id = 'casenum', time = 'time', name_shape = 'Response') +
  swimmer_arrows(swimmer_os, id = 'casenum', arrow_start = 'months', cont = 'censor') +
  annotate('text', x = 3, y = 31, label = 'Alive', size = 6) +
  annotate('text', x = 2, y = 31, label = sprintf('\u2192'), size = 12) +
  coord_flip(clip = 'off', ylim = c(0, 33)) +
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

###############################################################################################################

figure <- 
  plot_grid(
    plot_grid(lrfs_plot$plot, drfs_plot$plot, pfs_plot$plot, os_plot$plot,
              labels = 'AUTO'),
    swimmer_plot,
    marginal_progression_example,
    ncol = 1,
    rel_heights = c(2, 1.5, 1),
    labels = c('', 'E', 'F'))

save_plot('figures/figure_5.png', figure, ncol = 2.5, nrow = 5.5)
