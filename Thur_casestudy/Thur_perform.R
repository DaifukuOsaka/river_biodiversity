##run model
rm(list=ls())

getwd()
setwd("/Users/fuqing/rhistory/Thur_data")






#read Thur Data
EPT <- read.table("EPT_Genera_Luca.txt", header = TRUE, sep = "\t")
Kicknet <- read.table("Kicknet_Luca.txt", header = TRUE, sep = "\t")
head(EPT, n = 5)
head(Kicknet, n = 5)



#start with one genus





