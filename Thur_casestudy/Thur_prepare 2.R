#####prepare Thur
rm(list=ls())

setwd("/Users/fuqing/rhistory/Thur_data")

library(rivnet) #packageVersion‘0.3.1’
library(OCNet)
library(Rcpp)
library(terra)
library(readxl)
library(fields)
#library(eDITH) # from github.com/lucarraro/eDITH 
#here, we stick to the eDITH function
sourceCpp("evalConc2.cpp")

#simple plot
if (!file.exists("Thur.rda")){
  # extract Thur 
  Thur <- rivnet::extract_river(outlet = c(735010, 261530),
                                EPSG = 21781,
                                ext = c(700000, 770000, 220000, 270000),
                                z = 11,
                                displayUpdates=TRUE)
  # ground resolution about 54 meters; Latitude: 47°
  
  Thur <- rivnet::aggregate_river(Thur, thrA = 0.25e6, maxReachLength = 1000)
  
  # prepare hydrological data
  hydroData <- data.frame(x=c(723675, 727110, 718840, 737270),
                          y=c(252720, 247290, 248440, 251290),
                          w=c(    35,     20,    2.5,      8),
                          d=c(1.2034, 0.6784,     NA, 0.2500),
                          Q=c(52.8667,7.9083, 0.2595, 1.3392))
  
  # find coordinates of hydrological sites
  siteAG <- numeric(length(hydroData$x))
  for (i in 1:length(hydroData$x)){
    ss <- locate_site(hydroData$x[i], hydroData$y[i], Thur)
    siteAG[i] <- ss$AGnode
    #readline(prompt="Press [enter] to continue")
  }
  # build hydraulic model for Thur
  hyData <- data.frame(type=c(rep("w",4),rep("d",3),rep("Q",4)),
                       data=c(hydroData$w, hydroData$d[!is.na(hydroData$d)], hydroData$Q),
                       node=c(siteAG, siteAG[!is.na(hydroData$d)], siteAG))
  Thur <- hydro_river(hyData, Thur)
  
  # retrieve landcover raster
  landCovCH <- read.csv("/Users/fuqing/rhistory/Thur_data/ag-b-00.03-37-area-csv.csv", sep = ";")
  # landcover data can be downloaded at https://www.bfs.admin.ch/bfs/de/home/dienstleistungen/geostat/geodaten-bundesstatistik/boden-nutzung-bedeckung-eignung/arealstatistik-schweiz.assetdetail.20104753.html
  landCovCH <- terra::rast(data.frame(landCovCH$E, landCovCH$N, landCovCH$LU18_4),
                           type = "xyz", crs = "EPSG:2056") # convert into raster
  # Legend: 1-Urban; 2-Agriculture; 3-Forest; 4-Improductive
  landCovCH <- terra::project(landCovCH, crs("EPSG:21781"))
  
  # calculate landcover classes for the Thur
  Thur <- rivnet::covariate_river(landCovCH, Thur, overwrite=TRUE)
  rm(landCovCH)
  
  # retrieve geological classes
  geol_rast <- terra::rast("data/Clipped_Geology.asc")
  crs(geol_rast) <- "EPSG:21781"
  # simplfiy geological classes into fewer classes
  vv <- values(geol_rast)
  vv[vv==7 | vv==5 | vv==2] <- 1001 #"alluvial"
  vv[vv==30 | vv==57 | vv==108 | vv==48 | vv==45] <- 1002 # "alpine"
  vv[vv==14] <- 1003 # "loess"
  vv[vv==1 | vv==3 | vv==8 | vv==9 |  vv==20 |  vv==27] <- 1004 #"molasses"
  vv[vv==6] <- 1005 #"moraines"
  vv[vv==37 | vv==40 | vv==62 | vv==69 | vv==76 | vv==87] <- 1006 # "other"
  vv[vv==4] <- 1007 #"peat"
  vv[vv==35 | vv==43 | vv==97] <- 1008 #"scree"
  vv[vv==106] <- 1009 # "water"
  values(geol_rast) <- vv
  
  # calculate geological classes for the Thur
  Thur <- rivnet::covariate_river(geol_rast, Thur)
  
  # assign covariate names
  names(Thur$SC$locCov) <- names(Thur$SC$upsCov) <-c("urban", "agriculture", "forest", "improductive",
                                                     "alluvial","alpine","loess","molasses","moraines",
                                                     "other","peat","scree","water")
  
  save(Thur,file="Thur.rda",compress="xz")
} else {load("Thur.rda")}

# build geographical covariates (clusters of reaches that are near to each other)
siteClusters <- read_excel("data/siteClusters.xlsx")
geo_cluster <- numeric(Thur$AG$nNodes)
geoCov <- data.frame(matrix(0,Thur$AG$nNodes,dim(siteClusters)[1]))
for (i in 1:dim(siteClusters)[1]){
  geo_cluster[Thur$AG$upstream[[siteClusters$AG[i]]]] <- i
}
for (i in 1:dim(siteClusters)[1]){
  geoCov[geo_cluster==i, i] <- 1
}
names(geoCov) <- siteClusters$ID

covariates <- cbind(data.frame(logDrainageArea=log(Thur$AG$A), # streamOrder=Thur$AG$streamOrder,
                               elevation=Thur$AG$Z, slope=Thur$AG$slope),
                    Thur$SC$locCov[c("urban","agriculture","forest")],
                    Thur$SC$upsCov[c("alluvial","alpine","molasses","moraines","peat")], #"water", "loess", "scree"
                    geoCov)

# determine position of sampling sites
samplingSites <- read_excel("data/Coordinates_new.xlsx")
X <- samplingSites$X_new; Y <- samplingSites$Y_new
sampSiteAG <- length(X)
# manually edit coordinates of some sites to ensure that they are correctly snapped by locate_site
X[4] <- 726000;
X[7] <- 723700; Y[7] <- 254100
X[10] <- 733350
Y[15] <- 256300
X[21] <- 728400; Y[21] <- 251900
X[27] <- 737700; Y[27] <- 247200
X[39] <- 731300
Y[43] <- 239100
for (i in 1:length(X)){
  tmp <- locate_site(X[i], Y[i], Thur,showPlot=F)
  #title(sprintf("Site %d  -  X %d  -  Y %d",samplingSites$SiteID[i], X[i], Y[i]))
  sampSiteAG[i] <- tmp$AGnode
  #readline("Press Enter to continue:")
}
X_site <- X; Y_site <- Y;
site_key <- sort(samplingSites$SiteName, index.return=T); site_key <- site_key$ix
save(X_site,Y_site,site_key,sampSiteAG,file="samplingSitesInfo.rda")

#read Thur Data
eDNA <- read.table("EPT_Genera_Luca.txt", header = TRUE, sep = "\t")
#drop site no.22 for eDNA where kicknet not taken
eDNA <- eDNA[, !(colnames(eDNA) %in% c("M22.2", "M22.3", "M22.4"))]
Kicknet <- read.table("Kicknet_Luca.txt", header = TRUE, sep = "\t")

# ID of sampling sites
tmp <- as.numeric(substr(names(eDNA)[-1],2,3)) # the two digits after the 1st character of the site code are the site ID
                                               # e.g. site code "A07" -> site ID: 7
samplingSites_eDNA <- sampSiteAG[match(tmp,samplingSites$SiteID)]
# with match, we pass from SiteID 1--62 (with 38 missing), to indices 1--61. Then sampSiteAG gives the corresponding AG nodes
  
tmp <- Kicknet[,"site"]
samplingSites_kicknet <- sampSiteAG[match(tmp,samplingSites$SiteID)]

#define upstream&downstream
NodesDownstream <- samplingSites_kicknet[which(Thur$AG$A[samplingSites_kicknet] > median(Thur$AG$A[samplingSites_kicknet]))]
NodesUpstream <- setdiff(samplingSites_kicknet, NodesDownstream)

# 30 d 30 u
kicknet_downstream <- samplingSites_kicknet[samplingSites_kicknet %in% NodesDownstream]
kicknet_upstream <- samplingSites_kicknet[samplingSites_kicknet %in% NodesUpstream]
eDNA_downstream <- samplingSites_eDNA[samplingSites_eDNA %in% NodesDownstream]
eDNA_upstream <- samplingSites_eDNA[samplingSites_eDNA %in% NodesUpstream]

# data for Baetis; Caenis; Rhyacophila; Habroleptoides; Leuctra
genus <- c("Baetis","Caenis","Rhyacophila","Habroleptoides","Leuctra")#(most common speices, found at every kicknet site)
C_obs <- as.numeric(eDNA[which(eDNA[,1]==genus),-1])
D_obs <- Kicknet[,genus]

samplingIntensity <- "s30"
strategyType_eDNA <- c("ue10", "ue15", "ue20")#ue10:10 eDNA sites upstream
strategyType_kicknet <- c("uk10", "uk15", "uk20")

## randomly choose half sampling sites for calibration, half for validation
# pick sampling sites
samplingSitesList <- vector("list",length(samplingIntensity))
names(samplingSitesList) <- samplingIntensity
set.seed(200)
for (intensity in samplingIntensity){
  nSampling <- as.numeric(gsub("[^0-9.]", "",  intensity))
  nSampling_kicknet <- nSampling
  nSampling_eDNA <- nSampling
  for (strategyE in strategyType_eDNA){
    nUpstream_eDNA <- as.numeric(gsub("[^0-9.]", "",  strategyE))
    nDownstream_eDNA  <- nSampling_eDNA - nUpstream_eDNA
    for (strategyK in strategyType_kicknet){
      nUpstream_kicknet <- as.numeric(gsub("[^0-9.]", "",  strategyK))
      nDownstream_kicknet  <- nSampling_kicknet - nUpstream_kicknet
      samplingSitesList[[intensity]][[strategyE]][[strategyK]] <- 
        list(eDNA=matrix(0,ncol=nSampling_eDNA,nrow=10), 
             kicknet=matrix(0,ncol=nSampling_kicknet,nrow=10))
      samplingSitesList[[intensity]][[strategyE]][[strategyK]]$eDNA<- 
          c( sample(eDNA_downstream,nDownstream_eDNA), 
             sample(eDNA_upstream,nUpstream_eDNA))
      samplingSitesList[[intensity]][[strategyE]][[strategyK]]$kicknet<-     
          c(sample(kicknet_downstream,nDownstream_kicknet),
            sample(kicknet_upstream,nUpstream_kicknet))
    }
  }
}

# likelihood function ##river <- Thur
eval_loglik <- function(params, samplingSites_eDNA, samplingSites_kicknet, 
                        ss=ss, C_obs, D_obs, Thur, covariates){
  tau <- params[1]*3600 
  alpha <- params[2]    
  logp0 <- params[3]
  beta <- params[-c(1:3)] #remove tau, alpha & logp0
  p <- as.numeric(exp(logp0 + as.matrix(covariates) %*% beta)) 
  ss <- sort(Thur$AG$A,index.return=T); ss <- ss$ix
  conc <- evalConc2_cpp(Thur, ss, AS=Thur$AG$width*Thur$AG$leng,
                        tau, p, "AG")
  densi <- alpha*p

  # define flat prior to force ranges for parameter values
  logprior_values <- log(c(dunif(tau,1*3600,24*3600), # tau is bound between 1 h and 24 h
                   dunif(alpha,0,100), # alpha is bound between 0 and 100
                   dunif(logp0,-20,20), # logp0 is bound between -20 and 20
                    dnorm(beta,0,1)))  # gaussian prior for beta coefficients
  logprior <- sum(logprior_values)
  # geometric distribution for both read numbers and kicknet data
  # if prob = 1/(1+M), then the mean of the distribution is M
  # let's use geometric instead of Poisson distribution for kicknet because it's wider (observed kicknet data have high variance & zeros)
  # another advantage of the geometric distribution is that it doesn't require a parameter to define its variance
  loglik_values <- c(log(dgeom(C_obs, prob=1/(1+conc[samplingSites_eDNA]))), 
                  log(dgeom(D_obs, prob=1/(1+densi[samplingSites_kicknet])))) 
  loglik_values[loglik_values==-Inf] = -1e4 # set -Inf loglik values to very small value
  loglik <- sum(loglik_values)
  return(loglik+logprior) # we actually maximize the logposterior (=logprior + loglikelihood)
}

simulation_unknown <- function(Thur, C_obs, D_obs, samplingSites_kicknet, 
                               samplingSites_eDNA, covariates) {
  
  # Z-normalize covariates: each covariate will have mean=0 and sd=1
  for (i in 1:length(covariates)){
    covariates[,i] <- (covariates[,i]-mean(covariates[,i]))/sd(covariates[,i])}
  
  ss <- sort(Thur$AG$A,index.return=T); ss <- ss$ix
  # "ss" is a sequence of nodes in upstream-to-downstream direction (sorted by drainage area)
  # if we calculate concentration moving upstream to downstream, we do the calculation only once per node
  
  set.seed(1) 
  n.attempts <- 10 # try to optimize 10 times starting from different initial param values
                   # it might be increased to 100, if time allows
  ll_final_vec  <- numeric(n.attempts)
  
  for (ind in 1:n.attempts){
    initialParamValues <- c(runif(1,1,24), runif(1,0,100), runif(1,-20,20), # tau, alpha, logp0 (randomly chosen within the allowed ranges)
                            runif(length(covariates),-3,3)) # beta
    
    # we use the default (Nelder-Mead) optimization method
    # this method doesn't allow lower/upper bounds for the parameter
    # but we "cheat" it by using a prior distribution
    out <- optim(initialParamValues, eval_loglik,
                 samplingSites_kicknet=samplingSites_kicknet, 
                 samplingSites_eDNA=samplingSites_eDNA, ss=ss,
                 C_obs=C_obs, D_obs=D_obs, Thur=Thur, covariates=covariates,
                 control=list(fnscale=-1, maxit=1e6, trace=1)) # set trace=0 if you don't want to see output on console
   
    ll_final_vec[ind] <- out$value
    if (ind > 1){
      if (ll_final_vec[ind] > max(ll_final_vec[1:(ind-1)])){ 
        # only store the simulation with highest likelihood (actually, it's log posterior=likelihood+logprior)
        optList <- out
        iPV <- initialParamValues}
    } else {
      optList <- out
    iPV <- initialParamValues}
  }
  
  # re-run calibrated model
  tau <- optList$par[1]*3600
  alpha <- optList$par[2]
  logp0 <- optList$par[3]
  beta <- optList$par[-(1:3)]
  p <- as.numeric(exp(logp0 + as.matrix(covariates) %*% beta)) 
  ss <- sort(Thur$AG$A,index.return=T); ss <- ss$ix
  conc <- evalConc2_cpp(Thur, ss, AS=Thur$AG$width*Thur$AG$leng,
                        tau, p, "AG")
  densi <- alpha*p
  
  # prepare list for export
  out <- list(optList=optList,
              samplingSites_kicknet = samplingSites_kicknet,
              samplingSites_eDNA = samplingSites_eDNA,
              C_obs = C_obs,
              D_obs = D_obs,
              initialParamValues=iPV,
              p_mod = p, C_mod = conc, D_mod = densi) # export model output
  return(out)
}

# RUN THE SIMULATION
set.seed(200) 
for (species in genus){
  for(intensity in samplingIntensity){
    for (strategyE in strategyType_eDNA){
      for (strategyK in strategyType_kicknet){
        for (indSim in 1:10){
          
          C_obs <- as.numeric(eDNA[which(eDNA[,1]==genus),-1])
          D_obs <- Kicknet[,genus]
          
          ss <- sort(Thur$AG$A,index.return=T); ss <- ss$ix
          
          # choose sampling sites
          samplingSites = list()  
          samplingSites_kicknet <- samplingSitesList[[intensity]][[strategyE]][[strategyK]]$kicknet
          samplingSites_eDNA <- samplingSitesList[[intensity]][[strategyE]][[strategyK]]$eDNA
          
          file_name=paste0("results/results_",intensity,"_",strategyE,"_",strategyK,"_",as.character(indSim),".rda")
          cat("\n")
          cat(sprintf("indsim: %d   - intensity: %s   -  strategyType_eDNA: %s -   strategyType_kicknet: %s",
                      indSim,intensity,strategyE,strategyK,format(Sys.time(),"%b%d %H:%M")))
          
          if (file.exists(file_name)) {
            load(file_name)  
          } else {
              out <- NULL
              out <- simulation_unknown(Thur, C_obs, D_obs, samplingSites_kicknet, 
                                        samplingSites_eDNA, covariates)
              
              #out2 <- simulation_unknown(Thur, C_obs[1:10], D_obs, samplingSites_kicknet, 
                                         #samplingSites_eDNA[1:10], covariates)
              
              #export results data list
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
                strategyE = strategyE,
                strategyK = strategyK,
                intensity = intensity,
                indSim = indSim)
              save(result_data,file=file_name)}
          
             #access prediction skill
            
            
        } 
      }
    }
  } 
}



