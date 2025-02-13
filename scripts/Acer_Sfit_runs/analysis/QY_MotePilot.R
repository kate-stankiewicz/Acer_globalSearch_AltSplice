
# Script for processing physiology data for Acer AS
# Fluorescence-CBASS-Mote-2021
# authors K. Gomez-Campo and K. Stankiewicz


setwd("~/alt_splice/GS_Pilot_Acer/scripts/Acer_Sfit_runs/analysis")

library(ggplot2)
library(tidyverse)
library(dplyr)
geom_errorbar()
geom_linerange()
geom_pointrange()
geom_crossbar()
geom_errorbarh()

# read Table
Timepoints = read.csv("F-Yield.csv")
Timepoints$Timepoint1 <- as.factor(Timepoints$Timepoint1)
head(Timepoints)

# mean and the standard deviation
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

#Summarize the data
tpSummary <- data_summary(Timepoints, varname="Yield", 
                    groupnames=c("Treatment", "Timepoint1"))
# Convert dose to a factor variable
tpSummary$Timepoint1=as.factor(tpSummary$Timepoint1)
head(tpSummary)

# add hours as a column
tpSummary <- tpSummary %>%
  mutate(Hours = case_when(
    Timepoint1 == 0 ~ 0,
    Timepoint1 == 1 ~ 3,
    Timepoint1 == 2 ~ 6,
    Timepoint1 == 3 ~ 7,
    Timepoint1 == 4 ~ 12,
    Timepoint1 == 5 ~ 19,
    TRUE ~ NA_real_  # In case there are unexpected values
  ))

#PLOTS Quantum Yield of PSII kinetics as a function of time

#Create the plot
QY_t <- ggplot(tpSummary, aes(x = Hours, y = Yield, group = Treatment, color = Treatment)) + 
  geom_line(size = 1.5) +  # Increase the size of the lines
  geom_point(size = 3) +   # Optionally, increase the size of the points
  geom_errorbar(aes(ymin = Yield - sd, ymax = Yield + sd), width = 0.2, 
                size = 0.5,    # Thinner error bars
                position = position_dodge(0.05)) +
  labs(x = "Hours", y = "Quantum Yield of PSII") +  # Custom axis labels
  scale_color_manual(breaks = c("Control", "Low", "Medium", "High"), values = c("#FED976", "#FD8D3C", "#E31A1C", "#800026")) + theme_bw() +  geom_rect(aes(xmin = 6.5, xmax = 16.5, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.01)


#PLOTS Quantum Yield of PSII kinetics as a function of time


#RELATIVE TO CONTROLS
Timepoints = read.csv("F-Yield.csv")
Timepoints$Timepoint1 <- as.factor(Timepoints$Timepoint1)
head(Timepoints)

# mean and the standard deviation
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

#Summarize the data

tpSRel <- data_summary(Timepoints, varname="Rel.Control", 
                    groupnames=c("Treatment", "Timepoint1"))
# Convert Timepoint to a factor variable
tpSRel$Timepoint1=as.factor(tpSRel$Timepoint1)
head(tpSRel)
tpSRel

# Default line plot 
Rel <- ggplot(tpSRel, aes(x=Timepoint1, y=Rel.Control, group=Treatment, color=Treatment)) + 
  geom_line(size = 1.5) +  # Increase the size of the lines
  geom_point(size = 3) +   # Optionally, increase the size of the points
  ylim(0, 1.2) +
  geom_errorbar(aes(ymin = Rel.Control - sd, ymax = Rel.Control + sd), 
                width = 0.2, 
                size = 0.5,    # Thinner error bars
                position = position_dodge(0.05)) +
  labs(x = "Timepoints", y = "Relative Quantum Yield of PSII") + # Custom axis labels
  scale_color_manual(breaks = c("Control", "Low", "Medium", "High"), 
                         values = c("#757575", "#F9A825", "#E64A19", "#C62828"))



