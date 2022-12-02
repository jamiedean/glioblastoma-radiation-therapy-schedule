setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(magick)
library(pdftools)

source('code/custom_plot_theme.R')

################################################################################################################

treatment_strategy_diagram <- ggdraw() + 
  draw_image('figures/treatment_strategy_diagram.png',
             x = 1, width = 1, height = 1, hjust = 1)

compare_schedules_tid_interval_volume <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_tid_interval_volume.pdf', density = 600))

compare_schedules_qd_interval_volume <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_qd_interval_volume.pdf', density = 600))

compare_schedules_start_day_volume <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_start_day_volume.pdf', density = 600))

################################################################################################################

figure <-
  plot_grid(
    treatment_strategy_diagram,
    plot_grid(compare_schedules_tid_interval_volume,
              compare_schedules_qd_interval_volume,
              compare_schedules_start_day_volume,
              labels = c('B', 'C', 'D'),
              ncol = 1),
    labels = c('A', ''))

save_plot('figures/figure_s2.pdf', figure, ncol = 2, nrow = 3)
