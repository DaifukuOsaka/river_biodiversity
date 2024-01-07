rm(list=ls())

#set working directory
getwd()
setwd("/Users/fuqing/rhistory/simulation_optim_test") 

library(OCNet)
library(spam)
library(Rcpp)
library(tictoc)
library(rivnet) 
library(optimParallel)
source("eDITH_functions_RcppT3.R")

tic()

load("OCN_40.rda")
OCN$RN$downstreamPath <- NULL # save memory

tau_sol <- 4*3600
alpha_sol <- 1.5

modelType <- c("fixed","unknown","me")
samplingIntensity <- c("s16","s48")
samplingPreference <- c("k25","k50","k75")
strategyType_eDNA <- c("ue25", "ue50", "ue75")
strategyType_kicknet <- c("uk25", "uk50", "uk75")
distributionType <- c("R1","R2","R3","R4","R5","A1","A2","A3","A4","A5")

#k25:25% kicknet,75% eDNA
#ue25:25% eDNA upstream,75% downstream

# define hotspots
#R:scattered,rare taxon. A:evenly,abundant.
availableNodes <- setdiff(1:OCN$AG$nNodes, OCN$AG$outlet)
hotspotList <- vector("list",0)
set.seed(2); hotspotList$R1 <- sample(availableNodes, 5)
set.seed(5); hotspotList$R2 <- sample(availableNodes, 5)
set.seed(8); hotspotList$R3 <- sample(availableNodes, 5)
set.seed(9); hotspotList$R4 <- sample(availableNodes, 5)
set.seed(15); hotspotList$R5 <- sample(availableNodes, 5)
set.seed(11); hotspotList$A1 <- sample(availableNodes, 50)
set.seed(12); hotspotList$A2 <- sample(availableNodes, 50)
set.seed(13); hotspotList$A3 <- sample(availableNodes, 50)
set.seed(14); hotspotList$A4 <- sample(availableNodes, 50)
set.seed(15); hotspotList$A5 <- sample(availableNodes, 50)

#define upstream&downstream
Nodes <- setdiff(1:OCN$AG$nNodes,OCN$AG$outlet)
NodesDownstream <- intersect(which(OCN$AG$A > median(OCN$AG$A[Nodes])), Nodes)
NodesUpstream <- intersect(which(OCN$AG$A <= median(OCN$AG$A[Nodes])), Nodes)

# pick sampling sites ####
samplingSitesList <- vector("list",length(samplingIntensity))
names(samplingSitesList) <- samplingIntensity
k <- 0
for (intensity in samplingIntensity){
  #  nSampling <- as.numeric(gsub("[^0-9.]", "",  sampling))
  nSampling <- floor(as.numeric(gsub("[^0-9.]", "",  intensity))/100*OCN$AG$nNodes)
  nSampling <- min(nSampling, length(NodesDownstream), length(NodesUpstream)) # it can't be higher than the no. of available nodes
  for (preference in samplingPreference){
    nSampling_kicknet <- as.numeric(gsub("[^0-9.]", "",  preference))/100*nSampling
    nSampling_eDNA <- nSampling- nSampling_kicknet
    for (strategyE in strategyType_eDNA){
      weightUpstream_eDNA <- as.numeric(gsub("[^0-9.]", "",  strategyE))/100
      weightDownstream_eDNA  <- 1 - weightUpstream_eDNA
      for (strategyK in strategyType_kicknet){
        samplingSitesList[[intensity]][[preference]][[strategyE]][[strategyK]] <- 
          list(eDNA=matrix(0,ncol=nSampling_eDNA,nrow=10), 
               kicknet=matrix(0,ncol=nSampling_kicknet,nrow=10))
        weightUpstream_kicknet <- as.numeric(gsub("[^0-9.]", "",  strategyK))/100
        weightDownstream_kicknet <- 1 - weightUpstream_kicknet
        for (indSim in 1:10){
          k <- k+1
          set.seed(k) # seed changing everytime
          samplingSitesList[[intensity]][[preference]][[strategyE]][[strategyK]]$eDNA[indSim,] <- 
            c( sample(NodesDownstream,round(weightDownstream_eDNA*nSampling_eDNA)), 
               sample(NodesUpstream,round(weightUpstream_eDNA*nSampling_eDNA)))
          samplingSitesList[[intensity]][[preference]][[strategyE]][[strategyK]]$kicknet[indSim,] <-     
            c(sample(NodesDownstream,round(weightDownstream_kicknet*nSampling_kicknet)),
              sample(NodesUpstream,round(weightUpstream_kicknet*nSampling_kicknet))) 
        }
      }
    }
  }
}

# perform simulations ####
results_df <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(results_df) <- c("model","preference","strategyE","strategyK",
                          "intensity","distrib","indSim","PA","D","value","convergence")

# open parallel clusters
cl <- makeCluster(spec=detectCores(), outfile="")
setDefaultCluster(cl=cl)
clusterEvalQ(cl, {library("OCNet"); library("Rcpp")})
clusterExport(cl, "evalConc2_cpp")

set.seed(200) # this seed controls initial ParamValues 
for (model in modelType){
  for (preference in samplingPreference){
    for (strategyE in strategyType_eDNA){
      for (strategyK in strategyType_kicknet){
        for (intensity in samplingIntensity){
          for (distrib in distributionType){
            for (indSim in 1:10){
              
              p_sol <- numeric(OCN$AG$nNodes)
              
              hotspots <- hotspotList[[distrib]]
              halfHotspots <- vector("numeric",0)
              for (ih in hotspots){
                halfHotspots <- c(halfHotspots, OCN$AG$downNode[ih]) 
                halfHotspots <- c(halfHotspots, which(OCN$AG$downNode==ih))
                p_sol[ih] <- p_sol[ih] + 1
              }
              for (ihh in halfHotspots){
                p_sol[ihh] <- p_sol[ihh] + 0.5
              }
              # normalize p_sol(p_sim)
              p_sol <- p_sol/sum(p_sol)*sum(OCN$AG$leng * OCN$AG$width) / (OCN$AG$leng * OCN$AG$width)
              p_sol[OCN$AG$outlet] <- 0
              
              ss <- sort(OCN$AG$A,index.return=T); ss <- ss$ix 
              C_sol <- evalConc2_cpp(OCN, ss, tau_sol, p_sol, "AG")
              maxC <- sum(OCN$AG$width*OCN$AG$leng)/max(OCN$AG$velocity*OCN$AG$width*OCN$AG$depth) 
              C_sol <- C_sol/maxC #normalization
              D_sol <- alpha_sol*p_sol
              
              # add measurement error
              if (model == "me"){
                C_obs <- C_sol
                C_obs[runif(length(C_sol)) < exp(-C_sol/1)] <- 0 # non detection probability 
                C_obs <- exp(rnorm(length(C_obs),log(C_obs),0.5)) 
                D_obs <- round(D_sol)
                D_obs[runif(length(D_sol)) < exp(-D_sol/1)] <- 0 # non detection probability 
                D_obs <- round(exp(rnorm(length(D_obs),log(D_obs),0.5)) )
              } else {
                C_obs <- C_sol
                D_obs <- round(D_sol)
              }
              
              # choose sampling sites
              samplingSites = list()  
              samplingSites_kicknet <- samplingSitesList[[intensity]][[preference]][[strategyE]][[strategyK]]$kicknet[indSim,] 
              samplingSites_eDNA <- samplingSitesList[[intensity]][[preference]][[strategyE]][[strategyK]]$eDNA[indSim,]
              
              file_name=paste0("results_bestseeds/results_",model,"_",intensity,"_",preference,"_",strategyE,"_",strategyK,"_",distrib,"_",as.character(indSim),".rda")
              cat("\n")
              cat(sprintf("indsim: %d    -   model: %s - intensity: %s   -  preference: %s   -   strategyType_eDNA: %s -   strategyType_kicknet: %s  -  distrib: %s   -   %11s \n",
                          indSim,model,intensity,preference,strategyE,strategyK,distrib,format(Sys.time(),"%b%d %H:%M")))
              
              if (file.exists(file_name)) {
                load(file_name)  
              } else {
                # fit model
                if (model == "fixed"){
                  out <- NULL
                  out <- simulation_fixed(OCN, tau_sol, alpha_sol, p_sol, samplingSites_kicknet, samplingSites_eDNA) 
                  #result_data, csv file 
                  result_data <- list(
                    optList = out$optList,
                    initialParamValues = out$initialParamValues,
                    C_sol = out$C_sol,
                    C_obs = out$C_obs,
                    D_sol = out$D_sol,
                    D_obs = out$D_obs,
                    p_sol = p_sol,
                    samplingSites_kicknet = out$samplingSites_kicknet,
                    samplingSites_eDNA = out$samplingSites_eDNA,
                    preference = preference,
                    strategyE = strategyE,
                    strategyK = strategyK,
                    intensity = intensity,
                    distrib = distrib,
                    indSim = indSim)
                  save(result_data,file=file_name)
                } else {
                  out <- NULL
                  out <- simulation_unknown(OCN, tau_sol, alpha_sol, p_sol, samplingSites_kicknet, samplingSites_eDNA) 
                  result_data <- list(
                    optList = out$optList,
                    initialParamValues = out$initialParamValues,
                    C_sol = out$C_sol,
                    C_obs = out$C_obs,
                    D_sol = out$D_sol,
                    D_obs = out$D_obs,
                    p_sol = p_sol,
                    samplingSites_kicknet = out$samplingSites_kicknet,
                    samplingSites_eDNA = out$samplingSites_eDNA,
                    preference = preference,
                    strategyE = strategyE,
                    strategyK = strategyK,
                    intensity = intensity,
                    distrib = distrib,
                    indSim = indSim)
                  save(result_data,file=file_name)  
                }
              }
              
                ##data prep
                value <- result_data$optList$value
                convergence <- result_data$optList$convergence
                
                param <- result_data$optList$par
                if (model != "fixed"){
                  p <- param[-c(1,2)] # remove tau & alpha
                } else {p <- param}
                #p[outlet] <- NaN # remove outlet
                p_sol <- result_data$p_sol
                #p_sol[outlet] <- NaN
                C_sol <- result_data$C_sol
                C <- result_data$C_obs
                D_sol <- result_data$D_sol
                D <- result_data$D_obs
                
                #calculate PA&D (p_sol:p_sim, p:p_mol)
                p_a <- sum( (p_sol <= 0.25 & p <= 0.25) | (p_sol > 0.2 & p > 0.2) , na.rm=T) / 102
                strict<- sum( (p_sol <= 0.25 & p <= 0.25) | (p >= 0.75*p_sol & p <= 4/3*p_sol) , na.rm=T) / 102
                
                #add new row to data frame
                new_row <- data.frame(model = model,
                                      preference = preference,
                                      strategyE = strategyE,
                                      strategyK = strategyK,
                                      intensity = intensity,
                                      distrib = distrib,
                                      indSim =indSim,
                                      PA = p_a,
                                      D = strict,
                                      value = value,
                                      convergence = convergence)
                results_df <- rbind(results_df, new_row)
            }
          }
        }
      }
    }
  }
}
write.csv(results_df, "results_variant_bestseeds.csv", row.names = FALSE)
setDefaultCluster(cl=NULL); stopCluster(cl)

toc()