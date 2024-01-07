rm(list=ls())

getwd()
setwd("/Users/fuqing/rhistory/simulation_optim_test") 

library(ggplot2)
library(stringr)
library(tidyverse)
library(gridExtra)

library(stringr)
library(dplyr)

data <- read.csv("results.csv")

data <- data %>%
  mutate(nsites = as.numeric(str_extract(intensity, "\\d+")))

a1 <- 20
a2 <- 10
B1 <- 1
B2 <- 2

data <- data %>%
  mutate(
    preference_num = as.numeric(str_extract(preference, "\\d+")) / 100,
    nsites_k = nsites * preference_num,
    nsites_e = nsites - nsites_k,
    cost_e = a1 + B1 * (nsites - nsites * preference_num),
    cost_k = a2 + B2 * nsites * preference_num,
    cost_total = cost_e + cost_k
  )

data$strategy <- paste(data$strategyE, data$strategyK, sep="_")
data$disType <- str_extract(data$distrib, "[A-Za-z]")
head(data, 5)

strategy_colors <- c("ue25_uk25" = "darkred",
                     "ue25_uk50" = "red",
                     "ue25_uk75" = "lightcoral",
                     "ue50_uk25" = "darkblue",
                     "ue50_uk50" = "blue",
                     "ue50_uk75" = "lightblue",
                     "ue75_uk25" = "darkgreen",
                     "ue75_uk50" = "green",
                     "ue75_uk75" = "lightgreen")

ggplot(data, aes(x = cost_total, y = PA, color = strategy)) +
  geom_point() +
  scale_color_manual(values = strategy_colors) +
  theme_minimal() +
  labs(x = "Cost Total", y = "PA")

data_filtered <- data %>% filter(strategy == 'ue25_uk25')

ggplot(data_filtered, aes(x = cost_total, y = PA)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Cost Total", y = "PA")

