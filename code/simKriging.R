# Functions for estimating bias and SDs of Kriging predictions based on 
# the sampling designs from simulatePP.R

unitSquareKriging = function(nsims = 500, seed=1) {
  set.seed(seed)
  
  unifPreds = 1:nsims
  clustPreds = 1:nsims
  prefPreds = 1:nsims
  trueVals = 1:nsims
  
  #preallocate data vectors
  unifPreds1 = numeric(nsims)
  clustPreds1 = numeric(nsims)
  prefPreds1 = numeric(nsims)
  trueVals1 = numeric(nsims)
  unifPreds2 = numeric(nsims)
  clustPreds2 = numeric(nsims)
  prefPreds2 = numeric(nsims)
  trueVals2 = numeric(nsims)
  
  #generate Kriging predictor statistics
  unifPPs1 = genUnifPP(numSamples = nsims)
  clustPPs1 = genClusterPP2(numSamples = nsims, method="cutoff")
  prefPPs1 = genPreferentialPP2(numSamples = nsims, method="cutoff")
  GPs1 = prefPPs1$GPs
  unifPPs2 = genUnifPP(numSamples = nsims)
  clustPPs2 = genClusterPP2(numSamples = nsims,mu=1.515, sigmasq=.138, phi=.313, kappa=.5, beta=-2.198, tausq=0, method="cutoff")
  prefPPs2 = genPreferentialPP2(numSamples = nsims,mu=1.515, sigmasq=.138, phi=.313, kappa=.5, beta=-2.198, tausq=0, method="cutoff")
  GPs2 = prefPPs2$GPs
  lambda = .059/.138
  
  #add random noise to model 2 GPs:
  addNoise = function(x) {
    x$v = x$v + rnorm(length(x$v), sd=sqrt(.059))
    return(x)
  }
  GPs2 = lapply(GPs2, addNoise)
  
  for(i in 1:nsims) {
    print(paste0("generating ", i, "th set of simulations"))
    
    #get PP i of each type
    unifPP1 = unifPPs1[[i]]
    clustPP1 = clustPPs1$sims[[i]]
    prefPP1 = prefPPs1$sims[[i]]
    unifPP2 = unifPPs2[[i]]
    clustPP2 = clustPPs2$sims[[i]]
    prefPP2 = prefPPs2$sims[[i]]
    
    #get corresponding GP
    GP1 = GPs1[[i]]
    GP2 = GPs2[[i]]
    
    #compute Kriging estimators and true values of GP (model 1)
    unifStats1 = getKrigingPreds(unifPP1, GP1)
    clustStats1 = getKrigingPreds(clustPP1, GP1)
    prefStats1 = getKrigingPreds(prefPP1, GP1)
    unifPreds1[i] = unifStats1$preds
    clustPreds1[i] = clustStats1$preds
    prefPreds1[i] = prefStats1$preds
    trueVals1[i] = unifStats1$trueVals
    
    #compute Kriging estimators and true values of GP (model 2)
    unifStats2 = getKrigingPreds(unifPP2, GP2, lambda=lambda)
    clustStats2 = getKrigingPreds(clustPP2, GP2, lambda=lambda)
    prefStats2 = getKrigingPreds(prefPP2, GP2, lambda=lambda)
    unifPreds2[i] = unifStats2$preds
    clustPreds2[i] = clustStats2$preds
    prefPreds2[i] = prefStats2$preds
    trueVals2[i] = unifStats2$trueVals
  }
  
  zThresh = qnorm(.975)
  
  #get bias confidence intervals (model 1)
  unifBiasMean1 = mean(unifPreds1 - trueVals1)
  unifBiasSE1 = sd(unifPreds1 - trueVals1)/sqrt(nsims)
  unifBiasLB1 = unifBiasMean1 - zThresh*unifBiasSE1
  unifBiasUB1 = unifBiasMean1 + zThresh*unifBiasSE1
  unifBiasConf1 = c(unifBiasLB1, unifBiasUB1)
  clustBiasMean1 = mean(clustPreds1 - trueVals1)
  clustBiasSE1 = sd(clustPreds1 - trueVals1)/sqrt(nsims)
  clustBiasLB1 = clustBiasMean1 - zThresh*clustBiasSE1
  clustBiasUB1 = clustBiasMean1 + zThresh*clustBiasSE1
  clustBiasConf1 = c(clustBiasLB1, clustBiasUB1)
  prefBiasMean1 = mean(prefPreds1 - trueVals1)
  prefBiasSE1 = sd(prefPreds1 - trueVals1)/sqrt(nsims)
  prefBiasLB1 = prefBiasMean1 - zThresh*prefBiasSE1
  prefBiasUB1 = prefBiasMean1 + zThresh*prefBiasSE1
  prefBiasConf1 = c(prefBiasLB1, prefBiasUB1)
  
  #get bias confidence intervals (model 2)
  unifBiasMean2 = mean(unifPreds2 - trueVals2)
  unifBiasSE2 = sd(unifPreds2 - trueVals2)/sqrt(nsims)
  unifBiasLB2 = unifBiasMean2 - zThresh*unifBiasSE2
  unifBiasUB2 = unifBiasMean2 + zThresh*unifBiasSE2
  unifBiasConf2 = c(unifBiasLB2, unifBiasUB2)
  clustBiasMean2 = mean(clustPreds2 - trueVals2)
  clustBiasSE2 = sd(clustPreds2 - trueVals2)/sqrt(nsims)
  clustBiasLB2 = clustBiasMean2 - zThresh*clustBiasSE2
  clustBiasUB2 = clustBiasMean2 + zThresh*clustBiasSE2
  clustBiasConf2 = c(clustBiasLB2, clustBiasUB2)
  prefBiasMean2 = mean(prefPreds2 - trueVals2)
  prefBiasSE2 = sd(prefPreds2 - trueVals2)/sqrt(nsims)
  prefBiasLB2 = prefBiasMean2 - zThresh*prefBiasSE2
  prefBiasUB2 = prefBiasMean2 + zThresh*prefBiasSE2
  prefBiasConf2 = c(prefBiasLB2, prefBiasUB2)
  
  # get RMS and confidence intervals (model 1)
  unifRMSMean1 = sqrt(mean((unifPreds1 - trueVals1)^2))
  clustRMSMean1 = sqrt(mean((clustPreds1 - trueVals1)^2))
  prefRMSMean1 = sqrt(mean((prefPreds1 - trueVals1)^2))
  unifRMSE1 = sd((unifPreds1 - trueVals1)^2)/(2*sqrt(nsims)*sqrt(unifRMSMean1))
  clustRMSE1 = sd((clustPreds1 - trueVals1)^2)/(2*sqrt(nsims)*sqrt(clustRMSMean1))
  prefRMSE1 = sd((prefPreds1 - trueVals1)^2)/(2*sqrt(nsims)*sqrt(prefRMSMean1))
  unifRMSConf1 = c(unifRMSMean1 - zThresh*unifRMSE1, unifRMSMean1 + zThresh*unifRMSE1)
  clustRMSConf1 = c(clustRMSMean1 - zThresh*clustRMSE1, clustRMSMean1 + zThresh*clustRMSE1)
  prefRMSConf1 = c(prefRMSMean1 - zThresh*prefRMSE1, prefRMSMean1 + zThresh*prefRMSE1)
  
  # get RMS and confidence intervals (model 2)
  unifRMSMean2 = sqrt(mean((unifPreds2 - trueVals2)^2))
  clustRMSMean2 = sqrt(mean((clustPreds2 - trueVals2)^2))
  prefRMSMean2 = sqrt(mean((prefPreds2 - trueVals2)^2))
  unifRMSE2 = sd((unifPreds2 - trueVals2)^2)/(2*sqrt(nsims)*sqrt(unifRMSMean2))
  clustRMSE2 = sd((clustPreds2 - trueVals2)^2)/(2*sqrt(nsims)*sqrt(clustRMSMean2))
  prefRMSE2 = sd((prefPreds2 - trueVals2)^2)/(2*sqrt(nsims)*sqrt(prefRMSMean2))
  unifRMSConf2 = c(unifRMSMean2 - zThresh*unifRMSE2, unifRMSMean2 + zThresh*unifRMSE2)
  clustRMSConf2 = c(clustRMSMean2 - zThresh*clustRMSE2, clustRMSMean2 + zThresh*clustRMSE2)
  prefRMSConf2 = c(prefRMSMean2 - zThresh*prefRMSE2, prefRMSMean2 + zThresh*prefRMSE2)
  
  return(list(unifBiasConf1=unifBiasConf1, unifBiasConf2=unifBiasConf2, unifRMSConf1=unifRMSConf1, unifRMSConf2=unifRMSConf2, 
         clustBiasConf1=clustBiasConf1, clustBiasConf2=clustBiasConf2, clustRMSConf1=clustRMSConf1, clustRMSConf2=clustRMSConf2, 
         prefBiasConf1=prefBiasConf1, prefBiasConf2=prefBiasConf2, prefRMSConf1=prefRMSConf1, prefRMSConf2=prefRMSConf2))
}

# NOTE: if a window object is passed, it must have a ``bdry'' field with a 
# polygon of the boundary.  Further note that the ``res'' argument is ignored 
# unless the domainWin is passed as an argument.
getKrigingPreds = function(PP, GP, predCoords = rbind(c(.49,.49)), 
                               domainWin = NULL, res=101, lambda=0) {
  # get Kriging predictions at x0=(.49, .49)
  
    PPcoords = cbind(PP$x, PP$y)
    if(is.null(domainWin))
      GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
    else {
      GPCoords = attr(GP, "coords")
      GPCoords = GPCoords[in.poly(GPCoords, domainWin$bdry[[1]]),]
    }
    
    #####get the associated marks for the PP at the GP
    
    # round PP and prediction coordinates in coords to be on GP coordinate grid
    roundToRange = function(coordVec, coordRange) {
      inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
      inds = round(inds*(length(coordRange)-1)) + 1
      roundedCoordVec = coordRange[inds]
      return(roundedCoordVec)
    }
    if(is.null(domainWin)) {
      roundXPP = roundToRange(PPcoords[,1], GP$xcol)
      roundYPP = roundToRange(PPcoords[,2], GP$yrow)
      roundXPred = roundToRange(predCoords[,1], GP$xcol)
      roundYPred = roundToRange(predCoords[,2], GP$yrow)
    }
    else{
      roundXPP = roundToRange(PPcoords[,1], seq(domainWin$xrange[1], domainWin$xrange[2], l=res))
      roundYPP = roundToRange(PPcoords[,2], seq(domainWin$yrange[1], domainWin$yrange[2], l=res))
      roundXPred = roundToRange(predCoords[,1], seq(domainWin$xrange[1], domainWin$xrange[2], l=res))
      roundYPred = roundToRange(predCoords[,2], seq(domainWin$yrange[1], domainWin$yrange[2], l=res))
    }
    roundPPCoords = cbind(roundXPP, roundYPP)
    roundPredCoords = cbind(roundXPred, roundYPred)
    
    # find index of GP coords that are the same as the PP and prediction coords
    #GPCoords = attr(GP, "coords")
    findIndex = function(rCoords, gCoords) {
      return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
    }
    inds = findIndex(roundPPCoords, GPCoords)
    indsPred = findIndex(roundPredCoords, GPCoords)
    
    #get PP and prediction marks from GP using the indices found above
    #PPmarks = t(matrix(GP$v, ncol=length(GP$xcol)))[inds]
    PPmarks = GP$v[inds]
    predMarks = GP$v[indsPred]
    
    ##### Compute MLE
    # fixed parameters:
    kappa=1
    covArgs=list(Covariance="Matern", nu=kappa)
    
    # initial covariance guesses:
    phi=.2
    #sigma = sd(PPmarks)
    initCovParams = list(theta=phi)
    
    #do joint optimization using fields function
#     optArgs = list(method = "BFGS", 
#                       control=list(fnscale = -1, 
#                                    ndeps = rep(log(1.1), length(initCovParams)), 
#                                    reltol=1e-05, maxit=10))
    theta.grid=seq(.05, .6, l=20)
    PPcoords = jitter(PPcoords, amount=.01)
    opt = mKrig.MLE(x=PPcoords, y=PPmarks, lambda.profile=FALSE, 
                             lambda = lambda, par.grid=list(theta=theta.grid), 
                             cov.fun="stationary.cov", cov.args=covArgs, m=1)
    phi.MLE = opt$cov.args.MLE$theta
    MLE.ind = which(opt$par.grid == phi.MLE)
    sigma.MLE = sqrt(opt$summary[MLE.ind, 5])
    
    #make final mKrig fit at optimal parameters
    mKrigObj = mKrig(x=PPcoords, y=PPmarks, lambda = lambda, m=1, 
                      cov.args=list(Covariance="Matern", theta=phi.MLE, nu=kappa))
    
    #make predictions at prediction locations:
    preds = predict(mKrigObj, xnew=roundPredCoords)
    
    return(list(preds=preds, trueVals=predMarks, phi.MLE=phi.MLE))
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
# functions for optimizing GP likelihood





#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    

# Modification of mKrig.MLE.joint that assumes tausq (nugget) is negligible.  In 
# this case, assume lambda is very close to 0 and do not optimize it
My.mKrig.MLE.joint <- function(x, y, weights = rep(1, nrow(x)), 
                            lambda = 10^(-4), cov.params.guess=NULL, 
                            cov.fun="stationary.cov", cov.args=NULL, 
                            Z = NULL, optim.args=NULL, find.trA.MLE = FALSE, 
                            ..., verbose = FALSE) {
  
  #set default optim.args if necessary
  if(is.null(optim.args))
    optim.args = list(method = "BFGS", 
                      control=list(fnscale = -1, 
                                   ndeps = rep(log(1.1), length(cov.params.guess)), 
                                   reltol=1e-04, maxit=10))
  
  #check which optimization options the covariance function supports
  supportsDistMat = supportsArg(cov.fun, "distMat")
  
  #precompute distance matrix if possible so it only needs to be computed once
  if(supportsDistMat) {
    
    #Get distance function and arguments if available
    #
    Dist.fun= c(cov.args, list(...))$Distance
    Dist.args=c(cov.args, list(...))$Dist.args
    
    #If user left all distance settings NULL, use rdist with compact option.
    #Use rdist function by default in general.
    #
    if(is.null(Dist.fun)) {
      Dist.fun = "rdist"
      if(is.null(Dist.args))
        Dist.args = list(compact=TRUE)
    }
    
    distMat = do.call(Dist.fun, c(list(x), Dist.args))
  }
  
  #set cov.args for optimal performance if possible
  if(supportsDistMat)
    cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
  
  # these are all the arguments needed to call mKrig except lambda and cov.args
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, cov.fun=cov.fun), 
                  list(...))
  mKrig.args$find.trA = find.trA.MLE
  
  # output matrix to summarize results
  ncolSummary = 8 + length(cov.params.guess)
  summary <- matrix(NA, nrow = 1, ncol = ncolSummary)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", "GCV", "sigma.MLE", 
                                    "rho.MLE", "llambda.MLE", names(cov.params.guess), 
                                    "counts eval","counts grad"))
  
  # Define the objective function as a tricksy call to mKrig
  # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
  lnProfileLike.max <- -Inf
  temp.fn <- function(parameters) {
    # Seperate lambda from covariance parameters.
    # Optimization is over log-scale so exponentiate log-parameters.
    if(length(parameters) > 1) {
      otherParams = as.list(exp(parameters))
      names(otherParams) = names(cov.params.guess)
    }
    else
      otherParams = NULL
    
    #get all this eval's covariance arguments using the input parameters
    cov.args.temp = c(cov.args, otherParams)
    
    # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
    # assign to hold the last mKrig object
    hold <- do.call("mKrig", c(mKrig.args, list(lambda = lambda),
                               cov.args.temp))
    
    #save best mKrig object to global environment
    if(hold$lnProfileLike.FULL > lnProfileLike.max) {
      out <<- hold
      lnProfileLike.max = hold$lnProfileLike.FULL
    }
    hold = hold[c("rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
    
    # add this evalution to an object (i.e. here a matrix) in the calling frame
    temp.eval <- get("capture.evaluations")
    assign("capture.evaluations", rbind(temp.eval, c(lambda, parameters, unlist(hold))), envir = capture.env)
    return(hold$lnProfileLike.FULL)
  }
  
  #
  # optimize over covariance parameters
  
  # list of covariance arguments from par.grid with right names (some R arcania!)
  # note that this only works because 1) temp.fn will search in this frame for this object
  # par.grid has been coerced to a data frame so one has a concept of a row subscript.
  
  # set up matrix to store evaluations from within optim
  capture.evaluations <- matrix(NA, ncol = 4+length(cov.params.guess), nrow = 1,
                                dimnames = list(NULL, c("lambda", names(cov.params.guess), "rho.MLE",
                                                        "sigma.MLE", "lnProfileLike.FULL")))
  capture.env <- environment()
  
  # call to optim with initial guess (on log-scale)
  init.guess = log(unlist(cov.params.guess))
  look <- do.call(optim, c(list(par=init.guess), list(temp.fn), optim.args))
  
  #get optim results
  optim.counts <- look$counts
  llambda.opt <- log(lambda)
  lambda.opt <- lambda
  if(length(look$par) > 0) {
    params.opt <- exp(look$par)
    params.opt <- as.list(params.opt)
    names(params.opt) <- names(cov.params.guess)
  }
  else
    params.opt=NULL
  
  # call to 1-d search
  #            opt.summary     <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
  #            llambda.opt <- opt.summary$maximum
  #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
  # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
  lnLik.eval <- capture.evaluations[-1,]
  
  #exponentiate lambda and covariance parameters in lnLik.eval
  lnLik.eval[, 2:(1+length(look$par))] = exp(lnLik.eval[, 2:(1+length(look$par))])
  
  # calculate trace of best mKrig object if necessary
  #
  find.trA = list(...)$find.trA
  if(is.null(find.trA) || find.trA) {
    
    #get arguments for mKrig.trace
    iseed = list(...)$iseed
    NtrA = list(...)$NtrA
    
    #set iseed and NtrA to default values of mKrig if NULL
    if(is.null(iseed))
      iseed = 123
    if(is.null(NtrA))
      NtrA = 20
    
    #update best mKrig object with trace results
    out2 <- mKrig.trace(out, iseed, NtrA)
    out$eff.df <- out2$eff.df
    out$trA.info <- out2$trA.info
    np <- nrow(x)
    out$GCV <- (sum(out$residuals^2)/np)/(1 - out2$eff.df/np)^2
    if (NtrA < np)
      out$GCV.info <- (sum(out$residuals^2)/np)/(1 - out2$trA.info/np)^2
    else
      out$GCV.info <- NA
  }
  
  # save results of the best covariance model evaluation in a neat table
  summary[1, 1:ncolSummary] <- unlist(c(out$eff.df, out$lnProfileLike.FULL, 
                                        out$GCV, out$sigma.MLE.FULL, out$rho.MLE.FULL, llambda.opt, 
                                        params.opt, optim.counts))
  if (verbose) {
    cat("Summary: ", 1, summary[1, 1:ncolSummary], fill = TRUE)
  }
  
  #add summary table to output mKrig object and ensure it is still 
  #of class mKrig
  out = c(out, list(summary=summary, lnLik.eval=lnLik.eval))
  class(out) = "mKrig"
  
  return(out)
}





# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
mKrig.MLE <- function(x, y, weights = rep(1, nrow(x)), cov.fun="stationary.cov", cov.args = NULL, 
                      Z = NULL, par.grid = NULL, lambda = NULL, lambda.profile = TRUE, 
                      verbose = FALSE, relative.tolerance = 1e-04, ...) {
  
  #check which optimization options the covariance function supports
  supportsDistMat = supportsArg(cov.fun, "distMat")
  
  #precompute distance matrix if possible so it only needs to be computed once
  if(supportsDistMat) {
    
    #Get distance function and arguments if available.  Otherwise use 'dist' function
    #to compute upper triangle of distance matrix
    #
    Dist.fun= c(cov.args, list(...))$Distance
    Dist.args=c(cov.args, list(...))$Dist.args
    
    if(is.null(Dist.fun))
      Dist.fun = "dist"
    
    distMat = do.call(Dist.fun, c(list(x), Dist.args))
  }
  
  # mKrig.args has all the arguments needed to call mKrig except lambda and cov.args
  if(supportsDistMat)
    cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, cov.fun=cov.fun), 
                  list(...))
  
  lnProfileLike.max <- -1e+20
  
  # find NG --  number of parameters to try
  par.grid <- data.frame(par.grid)
  if (nrow(par.grid) == 0) {
    if (is.null(lambda)) {
      NG <- 1
    }
    else {
      NG <- length(lambda)
    }
  }
  else {
    NG <- nrow(par.grid)
  }
  
  # output matrix to summarize results
  summary <- matrix(NA, nrow = NG, ncol = 8)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", 
                                    "GCV", "sigma.MLE", "rho.MLE", "llambda.MLE", "counts eval", 
                                    "counts grad"))
  lambda.best <- NA
  
  # default for lambda is 1.0 for first value and exp(llambda.opt) for subsequent ones
  # this is controlled by NAs for lambda starting values.
  if (is.null(lambda)) {
    lambda <- rep(NA, NG)
  }
  
  # default starting value for lambda is 1 or log lambda is 0
  llambda.opt <- 0
  optim.counts <- c(NA, NA)
  lnLike.eval <- list()
  
  # Define the objective function as a tricksy call to mKrig
  # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
  temp.fn <- function(x) {
    # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
    # assign to hold only a few components returned by mKrig
    hold <- do.call("mKrig", c(mKrig.args, list(find.trA = FALSE, lambda = exp(x), 
                                                cov.args=c(cov.args.temp, cov.args)))
    )[c("lambda.fixed", "rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
    
    # add this evalution to an  object (i.e. here a matrix) in the calling frame
    temp.eval <- get("capture.evaluations")
    assign("capture.evaluations", rbind(temp.eval, unlist(hold)), 
           envir = capture.env)
    return(hold$lnProfileLike.FULL)
  }
  #
  # begin loop over covariance arguments
  for (k in 1:NG) {
    llambda.start <- ifelse(is.na(lambda[k]), llambda.opt, log(lambda[k]))
    
    # list of covariance arguments from par.grid with right names (some R arcania!)
    # note that this only works because 1) temp.fn will search in this frame for this object
    # par.grid has been coerced to a data frame so one has a concept of a row subscript.
    cov.args.temp <- as.list(par.grid[k, ])
    names(cov.args.temp) <- names(par.grid)
    
    #optimize over lambda if lambda.profile is TRUE
    if (lambda.profile) {
      # set up matrix to store evaluations from within optim
      capture.evaluations <- matrix(NA, ncol = 4, nrow = 1, 
                                    dimnames = list(NULL, c("lambda", "rho.MLE", 
                                                            "sigma.MLE", "lnProfileLike.FULL")))
      capture.env <- environment()
      
      # call to optim
      look <- optim(llambda.start, temp.fn, method = "BFGS", 
                    control = list(fnscale = -1, parscale = 0.1, 
                                   ndeps = 0.05, reltol = relative.tolerance))
      llambda.opt <- look$par
      optim.counts <- look$counts
      
      # call to 1-d search
      #            opt.summary <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
      #            llambda.opt <- opt.summary$maximum
      #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
      
      # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
      lnLike.eval <- c(lnLike.eval, list(capture.evaluations[-1, ]))
    }
    else {
      # no refinement for lambda so just save the the 'start' value as final one.
      llambda.opt <- llambda.start
    }
    
    # final fit at optimal value (or starting value if not refinement/maximization for lambda)
    obj <- do.call("mKrig", c(mKrig.args, list(lambda = exp(llambda.opt), 
                                               cov.args=c(cov.args.temp, cov.args))))
    if (obj$lnProfileLike.FULL > lnProfileLike.max) {
      lnProfileLike.max <- obj$lnProfileLike.FULL
      cov.args.MLE <- cov.args.temp
      lambda.best <- exp(llambda.opt)
    }
    
    # save results of the kth covariance model evaluation
    summary[k, 1:8] <- c(obj$eff.df, obj$lnProfileLike.FULL, 
                         obj$GCV, obj$sigma.MLE.FULL, obj$rho.MLE.FULL, llambda.opt, 
                         optim.counts)
    if (verbose) {
      cat("Summary: ", k, summary[k, 1:8], fill = TRUE)
    }
  }
  return(list(summary = summary, par.grid = par.grid, cov.args.MLE = cov.args.MLE, 
              mKrig.args = list(...), lambda.best = lambda.best, lambda.MLE = lambda.best, 
              call = match.call(), lnLike.eval = lnLike.eval))
}










