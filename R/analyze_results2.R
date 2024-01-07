rm(list=ls())

getwd()
setwd("/Users/fuqing/rhistory/simulation_optim_test") 

library(spam)
library(Rcpp)
library(OCNet)

library(ggplot2)
library(stringr)
library(tidyverse)
library(gridExtra)
source("eDITH_functions_RcppT.R")

data <- read.csv("results.csv")
data$strategy <- paste(data$strategyE, data$strategyK, sep="_")
data$disType <- str_extract(data$distrib, "[A-Za-z]")
head(data, n = 5)

######for model variants
data <- read.csv("results_variant_bestseeds.csv")
data$strategy <- paste(data$strategyE, data$strategyK, sep="_")
data$disType <- str_extract(data$distrib, "[A-Za-z]")
head(data, n = 5)
##remove error convergence && fn value >= -150
library(readr)
data <- subset(data, value >= -150)
data <- subset(data, convergence != 52)

#for strategy
color_mapping <- c("ue25_uk25" = "darkred",
                   "ue25_uk50" = "red",
                   "ue25_uk75" = "lightcoral",
                   "ue50_uk25" = "darkblue",
                   "ue50_uk50" = "blue",
                   "ue50_uk75" = "lightblue",
                   "ue75_uk25" = "darkgreen",
                   "ue75_uk50" = "green",
                   "ue75_uk75" = "lightgreen")
plot_r <- data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = factor(model, levels = c("fixed","unknown","me")), y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Scattered", x = "Model Variants", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = factor(model, levels = c("fixed","unknown","me")), y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Evenly", x = "Model Variants", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)

#for preference
#samplingPreference <- c("k25","k50","k75")#k25:25% kicknet,75% eDNA
color_mapping_preference <- c("k25" = "red",
                              "k50" = "blue",
                              "k75" = "green")
plot_r <- data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = factor(model, levels = c("fixed","unknown","me")), y = PA, fill = preference)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping_preference) +
  labs(title = "Scattered", x = "Model Variants", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = factor(model, levels = c("fixed","unknown","me")), y = PA, fill = preference)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping_preference) +
  labs(title = "Evenly", x = "Model Variants", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)



#####for outlierdots
data <- read.csv("results_dist.csv")
data$strategy <- paste(data$strategyE, data$strategyK, sep="_")
data$disType <- str_extract(data$distrib, "[A-Za-z]")
head(data, n = 5)

# 选择满足条件的数据
subset_data <- data[data$intensity == "s48" & grepl("R", data$disType), ]
print(subset_data)

color_mapping <- c("ue25_uk25" = "darkred",
                   "ue25_uk50" = "red",
                   "ue25_uk75" = "lightcoral",
                   "ue50_uk25" = "darkblue",
                   "ue50_uk50" = "blue",
                   "ue50_uk75" = "lightblue",
                   "ue75_uk25" = "darkgreen",
                   "ue75_uk50" = "green",
                   "ue75_uk75" = "lightgreen")

dodge_width <- 0.75

ggplot(subset_data, aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, position = position_dodge(dodge_width)) + 
  geom_jitter(position = position_dodge(dodge_width), size = 1.5, aes(color = strategy), alpha = 0.6) +
  scale_fill_manual(values = color_mapping) + 
  scale_color_manual(values = color_mapping) + 
  theme_minimal() +
  labs(title = "Notched Box Plot of PA by preference", 
       fill = "Strategy", color = "Strategy") +
  ylab("PA values") +
  xlab("Preference")

library(dplyr)

# 计算每个组的四分位数和IQR
outliers_data_dist <- subset_data %>%
  group_by(preference, strategy) %>%
  summarise(Q1 = quantile(PA, 0.25),
            Q3 = quantile(PA, 0.75),
            IQR = Q3 - Q1) %>%
  left_join(subset_data, by = c("preference", "strategy")) %>%
  filter(PA < (Q1 - 1.5 * IQR) | PA > (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR) 

print(outliers_data_dist)
write.csv(outliers_data_dist, "/Users/fuqing/rhistory/simulation_optim_test/outliers.csv", row.names = FALSE)

# 计算distrib中各类数据的数量
counts <- table(outliers_data_dist$distrib)
counts
##plot with selected seeds(R2,R5,R8,R9,R15)
results_data <- read.csv("results.csv") %>%
  filter(distrib %in% c("R2", "R5", "A1", "A2", "A3", "A4", "A5"))

results_dist_data <- read.csv("results_dist.csv") %>%
  filter(distrib %in% c("R8", "R9", "R15"))

combined_data <- bind_rows(results_data, results_dist_data)
combined_data$strategy <- paste(combined_data$strategyE, combined_data$strategyK, sep="_")
combined_data$disType <- str_extract(combined_data$distrib, "[A-Za-z]")
head(combined_data, n = 5)
#box plot
color_mapping <- c("ue25_uk25" = "darkred",
                   "ue25_uk50" = "red",
                   "ue25_uk75" = "lightcoral",
                   "ue50_uk25" = "darkblue",
                   "ue50_uk50" = "blue",
                   "ue50_uk75" = "lightblue",
                   "ue75_uk25" = "darkgreen",
                   "ue75_uk50" = "green",
                   "ue75_uk75" = "lightgreen")

plot_r <- combined_data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Scattered", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- combined_data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Evenly", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)




# 计算所占百分比

percentages <- prop.table(counts) * 100

# 将结果合并为一个数据框
result <- data.frame(Category = names(counts),
                     Count = as.numeric(counts),
                     Percentage = percentages)





# 创建基础的notched box plot，并使用 outlier.shape = NA 来隐藏离群点
p <- ggplot(subset_data, aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Notched Box Plot with Outliers Highlighted",
       y = "PA Value",
       x = "Preference",
       fill = "Strategy")

# 在基础的box plot上叠加单独的离群点
p + geom_point(data = outliers_data, aes(x = preference, y = PA), position = position_dodge(0.8), size = 3, shape = 16, color = "black")

#新表格处理分析
data <- read.csv("results_cov.csv")
data$strategy <- paste(data$strategyE, data$strategyK, sep="_")
data$disType <- str_extract(data$distrib, "[A-Za-z]")
head(data, n = 5)

subset_data <- data[data$intensity == "s48" & grepl("R", data$disType), ]
outliers_data <- subset_data %>%
  group_by(preference, strategy) %>%
  summarise(Q1 = quantile(PA, 0.25),
            Q3 = quantile(PA, 0.75),
            IQR = Q3 - Q1) %>%
  left_join(subset_data, by = c("preference", "strategy")) %>%
  filter(PA < (Q1 - 1.5 * IQR) | PA > (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR) 

print(outliers_data )
write.csv(outliers_data, "/Users/fuqing/rhistory/simulation_optim_test/outlierdots.csv", row.names = FALSE)

summary(data$value)
summary(subset_data$value)
summary(outliers_data$value)
summary(outliers_data$convergence)

count_values <- table(outliers_data$convergence)
values_to_extract <- c(0, 1, 10, 51, 52)
specific_counts <- count_values[as.character(values_to_extract)]
print(specific_counts)
#new plot: remove fn value>150
library(readr)
data <- subset(data, value >= -150)
data <- subset(data, convergence != 52)
color_mapping <- c("ue25_uk25" = "darkred",
                   "ue25_uk50" = "red",
                   "ue25_uk75" = "lightcoral",
                   "ue50_uk25" = "darkblue",
                   "ue50_uk50" = "blue",
                   "ue50_uk75" = "lightblue",
                   "ue75_uk25" = "darkgreen",
                   "ue75_uk50" = "green",
                   "ue75_uk75" = "lightgreen")

plot_r <- data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Scattered", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Evenly", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)





# boxplot PA 
ggplot(data, aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE, alpha = 0.6, outlier.shape = NA) +  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1.5) +  # PA=0.5
  facet_wrap(~intensity, scales = "free_x") +  # by intensity
  scale_fill_brewer(palette = "Set1") +
  theme_bw() + 
  labs(
    title = "Notched Box Plot by Preference and Intensity",
    x = "Preference",
    y = "PA",
    fill = "Strategy"
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 


#boxplot PA 
#Evenly&Scattered comparision
plot_r <- data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Scattered", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Evenly", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)

# boxplot PA with different col&shades
color_mapping <- c("ue25_uk25" = "darkred",
                   "ue25_uk50" = "red",
                   "ue25_uk75" = "lightcoral",
                   "ue50_uk25" = "darkblue",
                   "ue50_uk50" = "blue",
                   "ue50_uk75" = "lightblue",
                   "ue75_uk25" = "darkgreen",
                   "ue75_uk50" = "green",
                   "ue75_uk75" = "lightgreen")

plot_r <- data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_violin() +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Scattered", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = preference, y = PA, fill = strategy)) +
  geom_violin() +
  facet_wrap(~ intensity, ncol = 1) +
  geom_hline(yintercept = 0.5, color = "red", size = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Evenly", x = "Preference", y = "PA") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)

# boxplot D with different col&shades
color_mapping <- c("ue25_uk25" = "darkred",
                   "ue25_uk50" = "red",
                   "ue25_uk75" = "lightcoral",
                   "ue50_uk25" = "darkblue",
                   "ue50_uk50" = "blue",
                   "ue50_uk75" = "lightblue",
                   "ue75_uk25" = "darkgreen",
                   "ue75_uk50" = "green",
                   "ue75_uk75" = "lightgreen")

plot_r <- data %>%
  filter(disType == "R") %>%
  ggplot(aes(x = preference, y = D, fill = strategy)) +
  geom_boxplot(notch = TRUE) + #geom_violin() +
  facet_wrap(~ intensity, ncol = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Scattered", x = "Preference", y = "D") +
  theme(plot.title = element_text(hjust = 0.5))

plot_a <- data %>%
  filter(disType == "A") %>%
  ggplot(aes(x = preference, y = D, fill = strategy)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ intensity, ncol = 1) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Evenly", x = "Preference", y = "D") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot_r, plot_a, ncol=2)
