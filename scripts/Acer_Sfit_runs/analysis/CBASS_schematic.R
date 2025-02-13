# script plot CBASS schematic

# load libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)


exp_des <- data.frame(control= rep(0, 6), low = c(0, 3, 3, 0, 0, 0), medium = c(0, 5, 5, 0, 0, 0), high = c(0, 7, 7, 0, 0, 0), hours = c(0, 3, 6, 7, 12, 19), stringsAsFactors = F)

exp_des_long <- exp_des %>% tidyr::pivot_longer(cols = c(control,low,medium,high), names_to = "treatment", values_to = "temperature")

# set the colors
Cont_col="#FED976"
Low_col = "#FD8D3C"
Med_col = "#E31A1C"
High_col = "#800026"
exp_des_long$treatment<- factor(exp_des_long$treatment, levels = c( "control", "low", "medium", "high"))

plot_exp <- ggplot(exp_des_long, aes(x=hours, y = temperature, color = treatment, group= treatment)) + geom_line(linewidth= 1.4) + theme_bw() + scale_color_manual(values = c("control" = Cont_col, "low" = Low_col, "medium" = Med_col, "high" = High_col))


