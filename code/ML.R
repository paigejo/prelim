# maximize likelihood under standard geostatistical paradigm.  Use Nelder-Mead, 
# although the exact maximization method was not given in the paper.  Choose the 
# MLE parameters as the start value as a first test, but SDs not Vars and on log 
# scale.  Fit the parameters separately for the 1997 and 2000 log-transformed data.
# If not using MLEs as initial value, use:
# initParams = c(mean(log97), mean(log00), var(c(log00, log97)), .3, 0.001)
# But since we are on log scale, use instead
# initParams = c(mean(log97), mean(log00), log(sd(c(log00, log97))), log(.3), log(sqrt(10^-3)))
MLStandard = function(coords97, log97, coords00, log00, 
                      initParams = c(mu97MLE, mu00MLE, log(sqrt(sillMLE)), 
                                     log(scaleMLE), log(sqrt(nuggetMLE)))) {
  ##### decide on initial parameters if they aren't provided
  if(is.null(initParams)) {
    initParamsJoint = c(mean(log97), mean(log00), log(sd(c(log00, log97))), log(.3), log(sqrt(10^-3)))
    initParams00 = c(mean(log00), log(sd(log00)), log(.3), log(sqrt(10^-3)))
    initParams97 = c(mean(log97), log(sd(log97)), log(.3), log(sqrt(10^-3)))
  }
  else {
    initParamsJoint = initParams
    initParams00 = initParams[-1]
    initParams97 = initParams[-2]
  }
  
  ##### do precomputations:
  
  #precompute distance matrices
  distMat97 = rdist(coords97)
  distMat00 = rdist(coords00)
  
  ##### prepare log likelihood functions for maximization
  
  # make a summary tables showing parameter evaluations and likelihoods throughout
  # optimization
  summTable97 = matrix(nrow=0, ncol=5)
  colnames(summTable97) = c("mu97", "log(sigma)", "log(phi)", 
                          "log(tau)", "loglik97")
  summTable00 = summTable97
  colnames(summTable97) = c("mu00", "log(sigma)", "log(phi)", 
                            "log(tau)", "loglik00")
  
  getLik = function(params, is97 = TRUE) {
    #get parameters:
    mu = params[1]
    sigma = exp(params[2])
    phi = exp(params[3])
    tau = exp(params[4])
    
    #return very zero likelihood if parameters are negative
    if(any(c(sigma, phi, tau) < 0)) {
      if(is97)
        summTable97 <<- rbind(summTable97, c(params, -Inf))
      else
        summTable00 <<- rbind(summTable00, c(params, -Inf))
      return(log(10^-150))
    }
    
    if(is97) {
      distMat = distMat97
      logDat = log97
      coords = coords97
    }
    else {
      distMat = distMat00
      logDat = log00
      coords = coords00
    }
    
    # get Cholesky decomposition of Covariance matrices
    Sigma = Exp.cov(coords, distMat=distMat, theta=phi, onlyUpper=TRUE)*sigma^2
    invisible(.Call("addToDiagC", Sigma, as.double(rep(tau^2, nrow(Sigma))), nrow(Sigma)))
    SigmaU = chol(Sigma)
    
    # get log-likelihood assuming mean-centered data
    loglik = likGP(logDat - mu, SigmaU)
    
    # update summary table
    if(is97)
      summTable97 <<- rbind(summTable97, c(params, loglik))
    else
      summTable00 <<- rbind(summTable00, c(params, loglik))
    
    if(!is.finite(loglik)) {
      return(log(10^-150))
    }
    
    return(loglik)
  }
  
  summTableJoint = matrix(nrow=0, ncol=8)
  colnames(summTableJoint) = c("mu97", "mu00", "log(sigma)", "log(phi)", 
                            "log(tau)", "loglik97", "loglik00", "loglik")
  
  getLikJoint = function(params) {
    #get parameters:
    mu97 = params[1]
    mu00 = params[2]
    sigma = exp(params[3])
    phi = exp(params[4])
    tau = exp(params[5])
    
    #return very zero likelihood if parameters are negative
    if(any(c(sigma, phi, tau) < 0)) {
      if(any(c(sigma, phi, tau) < 0)) {
        summTableJoint <<- rbind(summTableJoint, c(params, -Inf))
        return(log(10^-150))
      }
    }
    
    # get Cholesky decomposition of Covariance matrices.  Use internal fields C function to
    # add nugget effect to diagonal of covariance matrices
    Sigma97 = Exp.cov(coords97, distMat=distMat97, theta=phi, onlyUpper=TRUE)*sigma^2
    invisible(.Call("addToDiagC", Sigma97, as.double(rep(tau^2, nrow(Sigma97))), nrow(Sigma97)))
    Sigma97U = chol(Sigma97)
    Sigma00 = Exp.cov(coords00, distMat=distMat00, theta=phi, onlyUpper=TRUE)*sigma^2
    invisible(.Call("addToDiagC", Sigma00, as.double(rep(tau^2, nrow(Sigma00))), nrow(Sigma00)))
    Sigma00U = chol(Sigma00)
    
    # get log-likelihood assuming mean-centered data
    loglik97 = likGP(log97 - mu97, Sigma97U)
    loglik00 = likGP(log00 - mu00, Sigma00U)
    loglik = loglik97 + loglik00
    
    if(!is.finite(loglik)) {
      summTableJoint <<- rbind(summTableJoint, c(params, c(loglik97, loglik00, -Inf)))
      return(log(10^-150))
    }
    
    # update summary table
    summTableJoint <<- rbind(summTableJoint, c(params, c(loglik97, loglik00, loglik)))
    
    return(loglik)
  }
  
  ##### perform optimization
  controls = list(fnscale=-1, parscale=.1)
  MLEs97 = optim(initParams97, getLik, control=list(fnscale=-1, parscale=rep(.1, 4)), 
                 hessian=TRUE, is97=TRUE)
  MLEs00 = optim(initParams00, getLik, control=list(fnscale=-1, parscale=rep(.1, 4)), 
                 hessian=TRUE, is97=FALSE)
  MLEsJoint = optim(initParamsJoint, getLikJoint, control=list(fnscale=-1, parscale=rep(.1, 5)), 
                    hessian=TRUE)
  
  #convert units back to standard (non-log) and Var rather than SD scales
  MLEs97$par[2:4] = exp(MLEs97$par[2:4])
  MLEs97$par[c(2,4)] = MLEs97$par[c(2,4)]^2
  names(MLEs97$par) = c("mu97", "sigmasq", "phi", "tausq")
  MLEs00$par[2:4] = exp(MLEs00$par[2:4])
  MLEs00$par[c(2,4)] = MLEs00$par[c(2,4)]^2
  names(MLEs00$par) = c("mu00", "sigmasq", "phi", "tausq")
  MLEsJoint$par[3:5] = exp(MLEsJoint$par[3:5])
  MLEsJoint$par[c(3,5)] = MLEsJoint$par[c(3,5)]^2
  names(MLEsJoint$par) = c("mu97", "mu00", "sigmasq", "phi", "tausq")
  
  initParams97[2:4] = exp(initParams97[2:4])
  initParams97[c(2,4)] = initParams97[c(2,4)]^2
  names(initParams97) = c("mu97", "sigmasq", "phi", "tausq")
  initParams00[2:4] = exp(initParams00[2:4])
  initParams00[c(2,4)] = initParams00[c(2,4)]^2
  names(initParams00) = c("mu00", "sigmasq", "phi", "tausq")
  initParamsJoint[3:5] = exp(initParamsJoint[3:5])
  initParamsJoint[c(3,5)] = initParamsJoint[c(3,5)]^2
  names(initParamsJoint) = c("mu97", "mu00", "sigmasq", "phi", "tausq")
  
  summTable97[,2:4] = exp(summTable97[,2:4])
  summTable97[,c(2,4)] = summTable97[,c(2,4)]^2
  colnames(summTable97) = c("mu97", "sigmasq", "phi", "tausq", "lik97")
  summTable00[,2:4] = exp(summTable00[,2:4])
  summTable00[,c(2,4)] = summTable00[,c(2,4)]^2
  colnames(summTable00) = c("mu00", "sigmasq", "phi", "tausq", "lik00")
  summTableJoint[,3:5] = exp(summTableJoint[,3:5])
  summTableJoint[,c(3,5)] = summTableJoint[,c(3,5)]^2
  colnames(summTableJoint) = c("mu97", "mu00", "sigmasq", "phi", "tausq", "loglik97", "loglik00", "loglik")
  
  return(list(summTable97=summTable97,summTable00=summTable00, 
              summTableJoint=summTableJoint, MLEs97 = MLEs97, 
              MLEs00 = MLEs00, MLEsJoint = MLEsJoint))
}

# test if sigma phi and tau are the same
likRatioTest = function(loglik97, loglik00, loglikJoint) {
  testStat = 2*(loglik97 + loglik00 - loglikJoint)
  dof = 3 # this is the number of extra parameters in alternative model
  pVal = 1 - pchisq(testStat, dof)
  return(pVal)
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


# Now do preferential sampling likelihood maximization


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

# maximize the parameters and using Monte Carlo likelihood calculations on a grid
# the paper uses 10,000 MC samples of S given Y per likelihood evaluation.  The 1997
# data is treated as preferential, and the 2000 data is not in the joint model.  It 
# looks like the authors assume conditional independence of the two realizations 
# given the parameters???  Or maybe they assume they take place at the same time???
MLPref = function(coords97, log97, coords00, log00, 
                  initParams = c(mu97MLE, mu00MLE, log(sqrt(sillMLE)), 
                                 log(scaleMLE), log(sqrt(nuggetMLE)), beta=betaMLE), 
                  res=60, nMCSamples = 10000, doPar=FALSE, nProc=4) {
  
  ##### decide on initial parameters if they aren't provided.  Note that the year 
  ##### 2000 model doesn't require a beta parameter since not preferential
  if(is.null(initParams)) {
    initParamsJoint = c(mean(log97), mean(log00), log(sd(c(log00, log97))), 
                        log(.3), log(sqrt(10^-3)), 0)
    initParams00 = c(mean(log00), log(sd(log00)), log(.3), log(sqrt(10^-3)))
    initParams97 = c(mean(log97), log(sd(log97)), log(.3), log(sqrt(10^-3)), 0)
    # or should initial guess for beta be related to difference in means?
  }
  else {
    initParamsJoint = initParams
    initParams00 = initParams[-c(1, 6)]
    initParams97 = initParams[-2]
  }
  
  ##### do precomputations:
  
  # precompute C matrix:
  # get range of points
  xrange = win$xrange
  yrange = win$yrange
  xs = seq(xrange[1], xrange[2], l=res)
  ys = seq(yrange[1], yrange[2], l=res)
  
  # get lattice coordinates, the coordinates for the grid approximation to the GP
  win.poly = matrix(unlist(win$bdry[[1]]), ncol=2)
  latticeCoords = make.surface.grid(list(x=xs, y=ys))
  
  # only take points within domain
  latticeCoords = latticeCoords[in.poly(latticeCoords, win.poly),]
  
  # round true coordinates to the lattice coordinates
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords97[,1], xs)
  roundY = roundToRange(coords97[,2], ys)
  roundCoords97 = cbind(roundX, roundY)
  
  # find indices of lattice coords corresponding to data coords
#   findIndex = function(coords, len=res) {
#     #convert coords to x and y grid indices
#     coords = coords*(len-1) + 1
#     xInd = coords[1]
#     yInd = coords[2]
#     return((yInd-1)*len + xInd)
#   }
#   inds = apply(roundCoords97, 1, findIndex, len=length(xrange))
  findIndex = function(coordPair) {
    coords.eq = function(latticePair) {
      return(all(coordPair == latticePair))
    }
    eq = apply(latticeCoords, 1, coords.eq)
    return(which(eq))
  }
  inds = unlist(apply(roundCoords97, 1, findIndex))
  
  #generate C matrix giving indices of X in X* (indices of rounded 
  #observation locations in grid).  Also make C transpose
  n = length(log97)
  N = nrow(latticeCoords)
  C = Matrix(0, nrow=n, ncol=N, sparse=TRUE)
  C[matrix(c(1:n, inds), ncol=2)] = 1
  Ct = t(C)
  
  # precompute distance matrices
  distMat97Lattice = rdist(latticeCoords)
  distMat00 = rdist(coords00)
  
  ##### prepare log likelihood functions for maximization
  
  # make a summary tables showing parameter evaluations and likelihoods throughout
  # optimization
  summTable97 = matrix(nrow=0, ncol=6)
  colnames(summTable97) = c("mu97", "log(sigma)", "log(phi)", 
                            "log(tau)", "beta", "loglik97")
  summTable00 = matrix(nrow=0, ncol=5)
  colnames(summTable00) = c("mu00", "log(sigma)", "log(phi)", 
                            "log(tau)", "loglik00")
  
  getLik00 = function(params) {
    #get parameters:
    mu = params[1]
    sigma = exp(params[2])
    phi = exp(params[3])
    tau = exp(params[4])
    
    #return very zero likelihood if parameters are negative
    if(any(c(sigma, phi, tau) < 0)) {
      summTable00 <<- rbind(summTable00, c(params, -Inf))
      return(log(10^-150))
    }
    
    distMat = distMat00
    logDat = log00
    coords = coords00
    
    # get Cholesky decomposition of Covariance matrices
    Sigma = Exp.cov(coords, distMat=distMat, theta=phi, onlyUpper=TRUE)*sigma^2
    invisible(.Call("addToDiagC", Sigma, as.double(rep(tau^2, nrow(Sigma))), nrow(Sigma)))
    SigmaU = chol(Sigma)
    
    # get log-likelihood assuming mean-centered data
    loglik = likGP(logDat - mu, SigmaU)
    
    # update summary table
    summTable00 <<- rbind(summTable00, c(params, loglik))
    
    if(!is.finite(loglik)) {
      return(log(10^-150))
    }
    
    return(loglik)
  }
  
  getLik97 = function(params) {
    # Evaluate log likelihood
    loglik = likPreferentialMC(params, latticeCoords, log97, distMat97Lattice, C, Ct, inds, 
                             res=res, nsims=nMCSamples, doPar, nProc)
    
    # update summary table
    summTable97 <<- rbind(summTable97, c(params, loglik))
    
    if(!is.finite(loglik))
      return(log(10^-150))
    else
      return(loglik)
  }
  
  summTableJoint = matrix(nrow=0, ncol=9)
  colnames(summTableJoint) = c("mu97", "mu00", "log(sigma)", "log(phi)", "log(tau)", "beta", 
                               "loglik97", "loglik00", "loglik")
  
  getLikJoint = function(params) {
    # Evaluate log likelihood
    loglik97 = getLik97(params[-2])
    loglik00 = getLik00(params[-c(1, 6)])
    loglik = loglik97 + loglik00
    
    # update summary table
    summTableJoint <<- rbind(summTableJoint, c(params, loglik97, loglik00, loglik), 
                             doPar=doPar, nProc=nProc)
    
    if(!is.finite(loglik))
      return(log(10^-150))
    else
      return(loglik)
  }
  
  ##### maximize likelihood
  controls = list(fnscale=-1, parscale=.1)
  MLEs97 = optim(initParams97, getLik97, control=list(fnscale=-1, parscale=rep(.1, 5)), 
                 hessian=TRUE)
  MLEs00 = optim(initParams00, getLik00, control=list(fnscale=-1, parscale=rep(.1, 4)), 
                 hessian=TRUE)
  MLEsJoint = optim(initParamsJoint, getLikJoint, control=list(fnscale=-1, parscale=rep(.1, 6)), 
                    hessian=TRUE)
  
  #convert units back to standard (non-log) and Var rather than SD scales
  MLEs97$par[2:4] = exp(MLEs97$par[2:4])
  MLEs97$par[c(2,4)] = MLEs97$par[c(2,4)]^2
  names(MLEs97$par) = c("mu97", "sigmasq", "phi", "tausq", "beta")
  MLEs00$par[2:4] = exp(MLEs00$par[2:4])
  MLEs00$par[c(2,4)] = MLEs00$par[c(2,4)]^2
  names(MLEs00$par) = c("mu00", "sigmasq", "phi", "tausq")
  MLEsJoint$par[3:5] = exp(MLEsJoint$par[3:5])
  MLEsJoint$par[c(3,5)] = MLEsJoint$par[c(3,5)]^2
  names(MLEsJoint$par) = c("mu97", "mu00", "sigmasq", "phi", "tausq", "beta")
  
  initParams97[2:4] = exp(initParams97[2:4])
  initParams97[c(2,4)] = initParams97[c(2,4)]^2
  names(initParams97) = c("mu97", "sigmasq", "phi", "tausq", "beta")
  initParams00[2:4] = exp(initParams00[2:4])
  initParams00[c(2,4)] = initParams00[c(2,4)]^2
  names(initParams00) = c("mu00", "sigmasq", "phi", "tausq")
  initParamsJoint[3:5] = exp(initParamsJoint[3:5])
  initParamsJoint[c(3,5)] = initParamsJoint[c(3,5)]^2
  names(initParamsJoint) = c("mu97", "mu00", "sigmasq", "phi", "tausq", "beta")
  
  summTable97[,2:4] = exp(summTable97[,2:4])
  summTable97[,c(2,4)] = summTable97[,c(2,4)]^2
  colnames(summTable97) = c("mu97", "sigmasq", "phi", "tausq", "beta", "lik97")
  summTable00[,2:4] = exp(summTable00[,2:4])
  summTable00[,c(2,4)] = summTable00[,c(2,4)]^2
  colnames(summTable00) = c("mu00", "sigmasq", "phi", "tausq", "lik00")
  summTableJoint[,3:5] = exp(summTableJoint[,3:5])
  summTableJoint[,c(3,5)] = summTableJoint[,c(3,5)]^2
  colnames(summTableJoint) = c("mu97", "mu00", "sigmasq", "phi", "tausq", "beta",
                               "loglik97", "loglik00", "loglik")
  
  return(list(summTable97=summTable97,summTable00=summTable00, 
              summTableJoint=summTableJoint, MLEs97 = MLEs97, 
              MLEs00 = MLEs00, MLEsJoint = MLEsJoint))
}

# make a function to calculate log-likelihood under preferential model.  Note that 
# this calculates the likelihood for any SINGLE dataset, not both 97 and 00 data.
# based on Eq. (9) and (10) from diggle (2010).  The input distMat should be for the
# gridCoords
likPreferentialMC = function(params, gridCoords, dat, distMat, C, Ct, inds, 
                           res=60, nsims=10000, doPar=FALSE, nProc=4) {
  # n and N as defined in paper: n is number of observations, N is number of pts in grid
  n = length(dat)
  N = length(gridCoords)
  
  #get parameters:
  mu97 = params[1]
  #mu00 = params[2]
  sigma = exp(params[2])
  phi = exp(params[3])
  tau = exp(params[4])
  beta = params[5] #(note: we don't care about beta until we calculate likelihood of X given S)
  
  #generate marginal covariance matrix of S, Sigma, for the given parameters
  Sigma = Exp.cov(gridCoords, theta=phi, distMat=distMat)*sigma^2
  
  #Now compute Sigma0, the covariance matrix of Y
  #Sigma0 = as.matrix(C %*% Sigma %*% Ct)
  Sigma0 = Sigma[inds, inds]
  
  # add nugget error to Sigma with fast internal function of fields in C
  invisible(.Call("addToDiagC", Sigma0, as.double(rep(tau^2, nrow(Sigma0))), nrow(Sigma0)))
  
  #Cholesky decomposition of Sigma0 required for likelihood
  Sigma0U = chol(Sigma0)
  
  # compute muc, expectation of Sj based on this formula:
  # Sigma %*% Ct %*% Sigma0^-1 %*% (dat - mu)
  vec = dat - mu97
  z = backsolve(Sigma0U, vec, transpose=TRUE)
  x = cbind(backsolve(Sigma0U, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
  #muc = as.numeric(Sigma %*% Ct %*% x)
  muc = Sigma[,inds] %*% x
  
  #calculate conditional mean of S given Y (just get sample mean)
  # muc = sigma^2*(dat - mu)/((sigma + tau)*sqrt(sigma^2 + tau^2))
  
  #estimate alpha using MOM estimator (I assume this is what Diggle does?)
  alpha = log(n/Area) - beta^2*sigma^2/2
  
  # calculate log-likelihood for given simulation, S, under preferential model
  getSimLik = function(S) {
    #simulate S given Y (authors denote this by Sj).  We want:
    # Sj = S + Sigma %*% Ct %*% Sigma0Inv %*% (dat - mu + Zs - C %*% S)
    # let x = Sigms0Inv %*% (dat - mu + Zs - C %*% S)
    # calculate x using Cholesky decomp of Sigma0
    Zs = rnorm(n)*tau
    #vec = dat - mu97 + Zs - C %*% S
    vec = dat - mu97 + Zs - S[inds]
    z = backsolve(Sigma0U, vec, transpose=TRUE)
    x = cbind(backsolve(Sigma0U, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
    #Sj = as.numeric(S + Sigma %*% Ct %*% x)
    Sj = S + Sigma[,inds] %*% x
    S0j = Sj[inds]
    
    # calculate lambda(x) given Y since needed in likelihood
    lambdas = exp(alpha + beta*Sj)
    
    # calculate likelihood of X given S|Y simulation.  Likelihood given by Eq. (8) 
    # in Diggle (2013). There might be faster ways to calculate it though.
    lambda = mean(lambdas) # N*lambda*dA should be about n (where dA = Area/N)
    # lik.XGivenS = prod(lambdas[inds])*lambda^(-n)
    lik.XGivenS = sum(log(lambdas[inds])) - n*log(Area*lambda) #use log-likelihood, not likelihood
    
    # calculate P(Y | S_{0j}) / P(S_{0j} | Y)
    if(tau == 0) {
      lik.frac = 1
    }
    else {
      # Calculate P(Y | S_{0j}):
      # Use iid normals since Y's conditionally indep given S
      lik.YGivenS0j = sum(log(dnorm((S0j + mu97 - dat)/tau))) #use log-likelihood not likelihood
      
      # Get P(S_{0j} | Y):
      # First calculate SigmaS0jGivenY =? tau^2*(I - tau^2*Sigma0Inv)
      UInv = backsolve(r=Sigma0U, x=diag(nrow=nrow(Sigma0U)))
      Sigma0Inv = UInv %*% t(UInv)
      SigmaS0jGivenY = tau^2*(diag(nrow=n) - tau^2*Sigma0Inv)
      
      # second calculate expectation of S0j|Y = Sigma %*% t(C) %*% Sigma0Inv %*% cbind(y)
      # we've already calculated this over all lattice points, just take right values
      muc0 = muc[inds]
      
      # then do Cholesky decomposition and solve using lik.GP
      SigmaS0jGivenYU = chol(SigmaS0jGivenY)
      lik.S0jGivenY = likGP(S0j - muc0, SigmaS0jGivenYU)
      
      # lik.frac = lik.YGivenS0j/lik.S0jGivenY
      lik.frac = lik.YGivenS0j - lik.S0jGivenY
    }
    
    # Calcualte log-likelihood of S0j marginally:
    # first calculate covariance matrix (same as Sigma0 but without nugget)
    Sigma0Tilde = Sigma0 - diag(tau^2, nrow=nrow(Sigma0))
    Sigma0TildeU = chol(Sigma0Tilde)
    lik.GP = likGP(S0j, Sigma0TildeU)
    
    lik = lik.XGivenS + lik.frac + lik.GP
    return(lik)
  }
  
  # compute likelihood for simulated GP S and its antithetic pair: 2mu_c - S
  getLikPair = function(i) {
    S = genS(params, gridCoords) ## TODO: S is list with 1 elmnt
    S = as.numeric(S$variable1)
    
    # muc = sigma^2*(dat - mu)/((sigma + tau)*sqrt(sigma^2 + tau^2))
    SAnti = 2*muc - S
    return(c(getSimLik(S), getSimLik(SAnti)))
  }
  
  #now we can compute all log-likelihoods using apply
  nPairs = ceil(nsims/2)
  
  if(!doPar) {
    allLogLiks = c(sapply(1:nPairs, getLikPair))
  }
  else {
    clust = makeCluster(nProc)
    clusterEvalQ(clust, setwd("~/git/prelim/code/"))
    clusterEvalQ(clust, source("loadAll.R"))
    allLogLiks = c(parSapply(clust, 1:nPairs, getLikPair))
    stopCluster(clust)
  }
  
  # final MC likelihood estimate is average of likelihoods, but return log-likelihood:
  return(log(mean(exp(allLogLiks))))
}

#calculate Gaussian Process log-likelihood given data covariance matrix Cholesky 
# decomposition assuming dat (data) is mean zero
likGP = function(dat, Sigma0U) {
  n = length(dat)
  
  #calcualte GP log likelihood (log likelihood of S0j)
  # (-1/2) t(y) %*% Sigma0^-1 y - (1/2) log det Sigma0 - (n/2) log(2pi)
  # define z = U^-1 %*% y, and x = L^-1 %*% U^-1 %*% y
  # then:
  z = backsolve(Sigma0U, dat, transpose=TRUE)
  x = cbind(backsolve(Sigma0U, z)) #make sure x is a column
  log.likGP = -(1/2) * rbind(dat) %*% x - sum(log(diag(Sigma0U))) - (n/2)*log(2*pi)
  lik.GP = log.likGP
  
  return(lik.GP)
}

#generate S from marginal distribution at data locations
#params: a vector of elements in order (mu97, mu00, logSill, logPhi, logTau, beta)
genS = function(params, gridCoords) {
  #get parameters:
  mu97 = params[1]
  #mu00 = params[2] #this is only called from getLikPair, which is only called from getSimLik
  sigmasq= exp(params[2])^2
  phi = exp(params[3])
  tausq = exp(params[4])^2
  beta = params[5]
  
  #use RMexp and RMsimulate to perform circulant embedding to quickly simulate S
  #TODO: generate multiple simulations to improve code efficiency
  RMobj = RMexp(var=sigmasq, scale=phi)
  sim = RFsimulate(RMobj, x=gridCoords[,1], y=gridCoords[,2]) 
  sim$variable1 = sim$variable1 #add mean? No!
  
  return(sim)
}

#based on http://onlinelibrary.wiley.com/doi/10.1002/9780470824566.app1/pdf and 
#http://blogs.sas.com/content/iml/2010/12/10/converting-between-correlation-and-covariance-matrices.html
hessianToCorrMat = function(hessMat) {
  covMat = solve(-hessMat)
  D = diag(1/sqrt(diag(covMat)))
  corrMat = D %*% covMat %*% D
  return(corrMat)
}








