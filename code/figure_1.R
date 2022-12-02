setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(magick)
library(pdftools)

source('code/custom_plot_theme.R')

################################################################################################################

compare_schedules_tid_10d_volume <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_tid_10d_volume.pdf', density = 600))

compare_schedules_hypo_tid_volume <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_hypo_tid_volume.pdf', density = 600))

compare_schedules_tid_10d_stem_cells <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_tid_10d_stem_cell_fraction.pdf', density = 600))

compare_schedules_hypo_tid_stem_cells <- 
  ggdraw() + 
  draw_image(magick::image_read_pdf(
    'figures/compare_schedules_hypo_tid_stem_cell_fraction.pdf', density = 600))

treatment_strategy_diagram <- ggdraw() + 
  draw_image('figures/treatment_strategy_diagram.png',
             x = 1, width = 1, height = 1, hjust = 1)

###############################################################################################################

figure <-
  plot_grid(compare_schedules_tid_10d_volume,
            compare_schedules_hypo_tid_volume,
            compare_schedules_tid_10d_stem_cells,
            compare_schedules_hypo_tid_stem_cells,
            labels = c('A', 'C', 'B', 'D'),
            ncol = 2) +
  theme(plot.background = element_rect(fill = 'white', color = NA))
save_plot('figures/figure_1.pdf', figure, ncol = 2, nrow = 2)
save_plot('figures/figure_1.tiff', figure, ncol = 2, nrow = 2)
