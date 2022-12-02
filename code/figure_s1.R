setwd('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/')

library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)

source('code/custom_plot_theme.R')

#################################################################################################################

data <- read.csv('data/preclinical_clinical_trial_correspondence.csv')

model1 <- lm(log(Human.OS.HR) ~ log(Mouse.OS.HR), data)
summary(model1)

model2 <- lm(log(Human.OS.HR) ~ log(Mouse.OS.HR.Upper.95..CI), data)
summary(model2)

p_fit <- ggplot(data, aes(x = log(Mouse.OS.HR), y = log(Human.OS.HR))) +
  geom_point(aes(x = log(Mouse.OS.HR), y = log(Human.OS.HR))) +
  geom_errorbar(aes(x = log(Mouse.OS.HR),
                    ymin = log(Human.OS.HR.Lower.95..CI), ymax = log(Human.OS.HR.Upper.95..CI))) +
  geom_errorbarh(aes(y = log(Human.OS.HR),
                     xmin = log(Mouse.OS.HR.Lower.95..CI), xmax = log(Mouse.OS.HR.Upper.95..CI))) +
  geom_abline(aes(intercept = coef(model1)[1], slope = coef(model1)[2], color = 'model1'),
              size = 1, show.legend = TRUE) +
  geom_abline(aes(intercept = coef(model2)[1], slope = coef(model2)[2], color = 'model2'),
              size = 1, show.legend = TRUE) +
  xlab(expression(ln(HR[mouse]^'measured'))) + ylab(expression(ln(HR[human]^'measured'))) +
  scale_color_manual(labels = c('Model 1', 'Model 2'),
                     values = c(model1 = plot_colors[1], model2 = plot_colors[2])) +
  theme(legend.position = c(0.03, 0.99), legend.justification = c(0, 1), legend.title = element_blank())
print(p_fit)

predictions <- data[data$Fit.Predict == 'predict',]

predictions$model1_pred <- 
  exp(predict(model1, newdata = data.frame(Mouse.OS.HR = predictions$Mouse.OS.HR)))
predictions$model2_pred <- 
  exp(predict(model2, newdata = data.frame(Mouse.OS.HR.Upper.95..CI = predictions$Mouse.OS.HR.Upper.95..CI)))
predictions$ensemble <- colMeans(rbind(predictions$model1_pred, predictions$model2_pred))
predictions$ensemble

p_predict <- ggplot(predictions) +
  geom_point(aes(x = log(Mouse.OS.HR), y = log(model1_pred), col = 'model1'), size = 3) +
  geom_point(aes(x = log(Mouse.OS.HR.Upper.95..CI), y = log(model2_pred), col = 'model2'), size = 3) +
  geom_text_repel(aes(x = log(Mouse.OS.HR), y = log(model1_pred), label = Label),
            hjust = -0.2, vjust = 0.5, size = 4) +
  geom_text_repel(aes(x = log(Mouse.OS.HR.Upper.95..CI), y = log(model2_pred), label = Label),
            hjust = -0.2, vjust = -0.0, size = 4) +
  xlab(expression(ln(HR[mouse]))) + ylab(expression(ln(HR[human]))) +
  scale_color_manual(labels = c('Model 1', 'Model 2'),
                     values = c(model1 = plot_colors[1], model2 = plot_colors[2])) +
  theme(legend.position = c(0.03, 0.99), legend.justification = c(0, 1), legend.title = element_blank())

p_predict_ensemble <- ggplot(predictions) +
  geom_point(aes(x = log(Mouse.OS.HR), y = log(ensemble), color = log(Mouse.OS.HR.Upper.95..CI)), size = 3) +
  scale_colour_gradient(name = expression(atop(ln(HR[mouse]^'measured'), 'Upper 95% CI')),
                        low = plot_colors[1], high = plot_colors[2]) +
  geom_text_repel(aes(x = log(Mouse.OS.HR), y = log(ensemble), label = Label),
                  hjust = -1, vjust = 0.0, size = 4) +
  xlab(expression(ln(HR[mouse]^'measured'))) + ylab(expression(ln(HR[human]^'predicted'))) +
  xlim(c(-2, 2)) + ylim(c(-0.4, 0.2))

#################################################################################################################

figure <- 
  plot_grid(p_fit,
            p_predict_ensemble + theme(legend.position = 'right'),
            ncol = 2, labels = 'AUTO', rel_widths = c(1, 1.3))
save_plot('figures/figure_s1.pdf', figure, ncol = 2, nrow = 1.5)
