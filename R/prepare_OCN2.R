rm(list=ls())

#set working directory
getwd()
setwd("/Users/fuqing/rhistory/OCN_Task") 

library(Rcpp)
library(OCNet)#version 1.0.1
library(spam) #package dealing with sparse matrices. Maybe not needed?
library(rivnet)#version 0.1.0.9000

source("eDITH_functions_Rcpp.R")


if (!file.exists("OCN_40.rda")){
#create an OCN 
  set.seed(1)
  OCN <- create_OCN(40, 40, typeInitialState = "T", 
                    coolingRate = 0.5, initialNoCoolingPhase = 0.1, cellsize = 50)
  OCN <- landscape_OCN(OCN,zMin = 100,slope0=0.02)
  thr <- find_area_threshold_OCN(OCN)
  OCN <- aggregate_OCN(OCN, thrA = thr$thrValues[max(which(abs(thr$nNodesAG - 100) == min(abs(thr$nNodesAG - 100))))])
  #AG:102 nodes / RN:278 nodes
  OCN <- paths_OCN(OCN, includePaths =TRUE, includeDownstreamNode = TRUE)
  OCN <- rivergeometry_OCN(OCN)
  OCN <- path_velocities_river(OCN)
  
  save(OCN,file="OCN_40.rda")
} else {
  load(file="OCN_40.rda")
}

#generate uniform &random distribution on AG level,sum(randomAG)==sum(uniformAG)
uniformAG <- 1+numeric(OCN$AG$nNodes)
randomAG <- runif(OCN$AG$nNodes)
randomAG <- randomAG/sum(randomAG)*sum(uniformAG)

ss <- sort(OCN$AG$A,index.return=T); ss <- ss$ix
# "ss" is a sequence of nodes in upstream-to-downstream direction (sorted by drainage area)
# if we calculate concentration moving upstream to downstream, we do the calculation only once per node

Conc_4h_unif <- evalConc2_cpp(OCN,ss,3600*4,uniformAG,"AG")
Conc_4h_random <- evalConc2_cpp(OCN,ss,3600*4,randomAG,"AG")

par(mfrow = c(1, 2))
draw_thematic_OCN(Conc_4h_unif,OCN,colLevels=c(1000,6000))
draw_thematic_OCN(Conc_4h_random,OCN,colLevels=c(1000,6000)) 

