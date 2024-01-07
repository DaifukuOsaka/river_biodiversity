rm(list=ls())

getwd()
setwd("/Users/fuqing/rhistory/Thur_data") 

library(spam)
library(Rcpp)
library(OCNet)

library(ggplot2)
library(stringr)
library(tidyverse)
library(gridExtra)
sourceCpp("evalConc2.cpp")

data <- read.csv("results_baetis.csv")
data$strategy <- paste(data$strategyE, data$strategyK, sep="_")

head(data, n = 5)
#for strategy
color_mapping <- c("ue10_uk10" = "darkred",
                   "ue10_uk15" = "red",
                   "ue10_uk20" = "lightcoral",
                   "ue15_uk10" = "darkblue",
                   "ue15_uk15" = "blue",
                   "ue15_uk20" = "lightblue",
                   "ue20_uk10" = "darkgreen",
                   "ue20_uk15" = "green",
                   "ue20_uk20" = "lightgreen")

genus <- c("Baetis","Caenis","Rhyacophila","Habroleptoides","Leuctra","Drusus","Ephemera","Perla")
spName = genus[2]

plot <- data %>%
  filter(species == spName) %>%
  ggplot(aes(x = strategy, y = loglik, color = factor(strategy))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_mapping) +
  labs(title = spName, x = "strategy", y = "loglik") +
  theme(plot.title = element_text(hjust = 0.5))

print(plot)

plot_r <- data %>%
  filter(species == spName) %>%
  ggplot(aes(x = strategy, y = loglik, fill = strategy)) +
  geom_boxplot(notch = FALSE) +
  facet_wrap(~ intensity, ncol = 1) +
  #geom_hline(yintercept = 0.5, color = "red", linewidth = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = spName, x = "strategy", y = "loglik") +
  theme(plot.title = element_text(hjust = 0.5))

print(plot_r)
