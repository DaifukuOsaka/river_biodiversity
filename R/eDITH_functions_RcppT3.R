# 2b) advanced algorithm which use C++
####faster version of eval_conc + eDITH
sourceCpp("evalConc2.cpp")

eval_loglik <- function(params, samplingSites_eDNA, samplingSites_kicknet, ConcAG, DensAG, OCN){
  tau <- (params[1] + 1)*3600 # tau can't be lower than 1 h
  alpha <- params[2] + 0.5    # alpha can't be lower than 0.5    
  p <- params[-c(1,2)] #remove tau & alpha
  ss <- sort(OCN$AG$A,index.return=T); ss <- ss$ix
  conc <- evalConc2_cpp(OCN, ss, tau, p, "AG")
  maxC <- sum(OCN$AG$width*OCN$AG$leng)/max(OCN$AG$velocity*OCN$AG$width*OCN$AG$depth) 
  conc <- conc/maxC#nomalization
  densi <- alpha*p
  loglik <- sum(log(dnorm(conc[samplingSites_eDNA], ConcAG[samplingSites_eDNA], 2)) +
                  log(dpois(DensAG[samplingSites_kicknet], densi[samplingSites_kicknet]))) # modelled density ("densi") is the mean; this imposes that densi[i] can't be 0 when DensAG[i]>0
  if (!isFALSE(loglik < -1e100 )){ # if loglik is NA or -Inf, set random very low value
    loglik <- -(1+runif(1))*1e100
  }
  return(loglik)
}

simulation_unknown <- function(OCN, tau_sol, alpha_sol, p_sol, samplingSites_kicknet, samplingSites_eDNA) {
  ss <- sort(OCN$AG$A,index.return=T); ss <- ss$ix
  # "ss" is a sequence of nodes in upstream-to-downstream direction (sorted by drainage area)
  # if we calculate concentration moving upstream to downstream, we do the calculation only once per node
  C_sol <- evalConc2_cpp(OCN, ss, tau_sol, p_sol, "AG")
  maxC <- sum(OCN$AG$width*OCN$AG$leng)/max(OCN$AG$velocity*OCN$AG$width*OCN$AG$depth) 
  C_sol <- C_sol/maxC #nomalization
  C_obs <- C_sol#assume no measurement error
  D_sol <- alpha_sol*p_sol
  D_obs <- round(D_sol)
  initialParamValues <- c(runif(1,0,7), runif(1,0,5), runif(OCN$AG$nNodes,0,5))
  optList <- optim(initialParamValues, eval_loglik,
                           method = "L-BFGS-B",
                           lower=0, upper=250,
                           samplingSites_kicknet=samplingSites_kicknet, samplingSites_eDNA=samplingSites_eDNA,
                           ConcAG=C_obs, DensAG=D_obs, OCN=OCN,
                           control=list(fnscale=-1, maxit=200,trace=1, REPORT=10))
  out <- list(optList=optList,
              samplingSites_kicknet = samplingSites_kicknet,
              samplingSites_eDNA = samplingSites_eDNA,
              C_sol = C_sol,
              C_obs = C_obs,
              D_sol = D_sol,
              D_obs = D_obs,
              initialParamValues=initialParamValues
              )
  return(out)
}

eval_loglik_taufixed_alphafixed <- function(params, tau_sol, alpha_sol, samplingSites_eDNA, samplingSites_kicknet, ConcAG, DensAG, OCN){
  p <- params
  ss <- sort(OCN$AG$A,index.return=T); ss <- ss$ix
  conc <- evalConc2_cpp(OCN, ss, tau_sol, p, "AG")
  maxC <- sum(OCN$AG$width*OCN$AG$leng)/max(OCN$AG$velocity*OCN$AG$width*OCN$AG$depth) 
  conc <- conc/maxC#nomalization
  densi <- alpha_sol*p
  loglik <- sum(log(dnorm(conc[samplingSites_eDNA], ConcAG[samplingSites_eDNA], 2)) +
                  log(dpois(DensAG[samplingSites_kicknet], densi[samplingSites_kicknet]))) # modelled density ("densi") is the mean; this imposes that densi[i] can't be 0 when DensAG[i]>0
  if (!isFALSE(loglik < -1e100 )){ # if loglik is NA or -Inf, set random very low value
    loglik <- -(1+runif(1))*1e100
  }
  return(loglik)
}

simulation_fixed <- function(OCN, tau_sol, alpha_sol, p_sol, samplingSites_kicknet, samplingSites_eDNA){
  ss <- sort(OCN$AG$A,index.return=T); ss <- ss$ix
  # "ss" is a sequence of nodes in upstream-to-downstream direction (sorted by drainage area)
  # if we calculate concentration moving upstream to downstream, we do the calculation only once per node
  C_sol <- evalConc2_cpp(OCN, ss, tau_sol, p_sol, "AG")
  maxC <- sum(OCN$AG$width*OCN$AG$leng)/max(OCN$AG$velocity*OCN$AG$width*OCN$AG$depth) 
  C_sol <- C_sol/maxC #nomalization
  C_obs <- C_sol#assume no measurement error
  D_sol <- alpha_sol*p_sol
  D_obs <- round(D_sol)
  initialParamValues <- runif(OCN$AG$nNodes,0,5)
  optList <- optim(initialParamValues, eval_loglik_taufixed_alphafixed,
                           method = "L-BFGS-B",
                           lower=0, upper=250,
                           tau_sol=tau_sol, alpha_sol=alpha_sol, 
                           samplingSites_kicknet=samplingSites_kicknet, samplingSites_eDNA=samplingSites_eDNA,
                           ConcAG=C_obs, DensAG=D_obs, OCN=OCN,
                           control=list(fnscale=-1, maxit=200,trace=1, REPORT=10))
  out <- list(optList=optList,
              samplingSites_kicknet = samplingSites_kicknet,
              samplingSites_eDNA = samplingSites_eDNA,
              C_sol = C_sol,
              C_obs = C_obs,
              D_sol = D_sol,
              D_obs = D_obs,
              initialParamValues=initialParamValues
             )
  return(out)
}