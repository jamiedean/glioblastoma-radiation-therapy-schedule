setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(magick)
library(pdftools)

source('code/custom_plot_theme.R')

################################################################################################################

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
  plot_grid(compare_schedules_tid_interval_volume,
            compare_schedules_qd_interval_volume,
            compare_schedules_start_day_volume,
            labels = 'AUTO',
            ncol = 3)

save_plot('figures/figure_s2.pdf', figure, ncol = 3, nrow = 1)
