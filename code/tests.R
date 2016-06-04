##### functions for testing likelihood evaluation

##### Their likelihood: L_(theta) = p(x, y) ~ (1/m) sum([X | Sj] ([Y|S0j] / [S0j|Y]) [S0j])

##### mcmcse: estimates mcmc and mc error

##### simulate data under 2000 model and fit joint likelihood of simulated data
##### look into additional variance of simulation study kriging predictions
##### email Peter about predictive distribution
##### check naive likelihood with simulated data
##### compare multivariate normal likelihood with that of some package

##### test 8) make VG bias plots under model with same parameters as fitted MLEs

##### test 9) simulate data on Galicia domain (preferentially and non-preferentially) and try to find MLEs

##### test 1): p(y, x, s) = p(y | X, S) p(X | S) p(S)
# this is the naive way of calculating the likelihood.  Doesn't work when tau is close to 0 since 
# the number of samples required explodes due to the likelihood being so peaked.

naivePrefLik = function(coords97, log97, coords00, log00, 
                        initParams = c(mu97MLE, mu00MLE, log(sqrt(sillMLE)), 
                                       log(scaleMLE), log(sqrt(nuggetMLE)), beta=betaMLE), 
                        res=60, nMCSamples = 10000) {
  
}



##### test 2) p(y,s) = int_x p(x, y | S) p(S) = int_x p(S, x | y) p(y) (use S0's?)
# this tests to see if my p(S | y) agrees with my P(y | S)
# not really sure how to implment this.  bootstrap for marginal of y?  what is p(S, x | y)?

#####  other tests:
#
# 3) do predictive simulations given MLEs to see if the Sj simulations are right
# also generate standard predictions from fields and and Sj, see if they agree
#
testPredsStandard = function(nsim=500, mu=1.55157286, sigma=sqrt(0.13635714), 
                             phi=0.30503468, tau=sqrt(0.05239249), nu=.5) {
  
  #get prediction coordinates (win and predRes are globally defined variables)
  xrange = win$xrange
  yrange = win$yrange
  xs = seq(xrange[1], xrange[2], l=predRes)
  ys = seq(yrange[1], yrange[2], l=predRes)
  
  # make sure they're in spatial domain
  win.poly = matrix(unlist(win$bdry[[1]]), ncol=2)
  latticeCoords = make.surface.grid(list(x=xs, y=ys))
  latticeCoords = latticeCoords[in.poly(latticeCoords, win.poly),]
  
  # get Kriging object 
  lambda = tau^2/(tau^2 + sigma^2)
  krigObj = mKrig(coords97, log97, cov.args=list(range=phi, Covariance="Exponential"), m=1, lambda=lambda)
  
  #generate simulations from predictive distribution
  preds = exp(sim.mKrig.approx(krigObj, predictionPoints=latticeCoords, M=nsim)$Ensemble)
  
  # get 5%, 50%, and 95% quantiles
  quant5 = apply(preds, 1, quantile, probs=.05)
  quant50 = apply(preds, 1, quantile, probs=.5)
  quant95 = apply(preds, 1, quantile, probs=.95)
  
  #get distributions for proportion larger than 3, 5, and 7 micrograms of lead
  propBigger = function(dat, thresh = 3) {
    mean(dat > thresh)
  }
  above3Prop = apply(preds, 2, propBigger, thresh=3)
  above5Prop = apply(preds, 2, propBigger, thresh=5)
  above7Prop = apply(preds, 2, propBigger, thresh=7)
  return(list(quant5=quant5, quant50=quant50, quant95=quant95, above3Prop=above3Prop, above5Prop=above5Prop, 
              above7Prop=above7Prop, simCoords=latticeCoords))
}

# fieldsPredsMLE = testPredsStandard(nsim=10000, mu=mu97MLE, sigma=sqrt(varianceMLE), phi=scaleMLE, tau=sqrt(nuggetMLE))
# makePlots(fieldsPredsMLE)

# generate predictions directly from multivariate normal formulation and see if it matches Sj's
genPredictionsTest = function(nsim=10000, nx=predRes, ny=predRes, mu=mu97MLE, sigma=sqrt(varianceMLE), 
                              phi=scaleMLE, tau=sqrt(nuggetMLE), nu=.5) {
  
  # generate the simulation coordinates
  xrange = seq(win$xrange[1], win$xrange[2], l=nx)
  yrange = seq(win$yrange[1], win$yrange[2], l=ny)
  simCoords = make.surface.grid(list(x=xrange, y=yrange))
  
  # simulate S from marginal distribution (nu is assumed to be 1/2 in Matern covariance)
  # use intrinsic method instead of cutoff embedding due to rectangular grid (faster)
  # to save space, will overwrite Ss with Sjs then predictions
  params = c(NA, log(sigma), log(phi), log(tau))
  preds = genS(params, gridCoords=simCoords, nsim=nsim)
  
  # filter simulations and coordinates based on data domain (better to simulate before this because it's efficient to 
  # simulate on square grid)
  win.poly = matrix(unlist(win$bdry[[1]]), ncol=2)
  winInds = in.poly(simCoords, win.poly)
  simCoords = simCoords[winInds,]
  preds = preds[winInds,]
  
  # n and N as defined in paper: n is number of observations, N is number of pts in grid
  n = length(log97)
  N = nrow(simCoords)
  
  # round true coordinates to the lattice coordinates
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords97[,1], xrange)
  roundY = roundToRange(coords97[,2], yrange)
  roundCoords97 = cbind(roundX, roundY)
  
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords97, simCoords)
  
  #generate marginal covariance matrix of S, Sigma, for the given parameters
  distMat = rdist(simCoords)
  Sigma = Exp.cov(simCoords, theta=phi, distMat=distMat)*sigma^2
  
  #Now compute Sigma0, the covariance matrix of Y
  #Sigma0 = as.matrix(C %*% Sigma %*% Ct)
  Sigma0 = Sigma[inds, inds]
  Sigma10 = Sigma[,inds]
  Sigma01 = Sigma[inds,]
  
  # add nugget error to Sigma with fast internal function of fields in C
  invisible(.Call("addToDiagC", Sigma0, as.double(rep(tau^2, nrow(Sigma0))), nrow(Sigma0)))
  
  #Cholesky decomposition of Sigma0 required for likelihood
  Sigma0U = chol(Sigma0)
  
  # compute muc, expectation of Sj based on this formula:
  # Sigma %*% Ct %*% Sigma0^-1 %*% (dat - mu)
  vec = log97 - mu
  z = backsolve(Sigma0U, vec, transpose=TRUE)
  x = cbind(backsolve(Sigma0U, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
  #muc = as.numeric(Sigma %*% Ct %*% x)
  muc = Sigma10 %*% x
  
  #simulate S given Y (authors denote this by Sj).  We want:
  # Sj = S + Sigma %*% Ct %*% Sigma0Inv %*% (dat - mu + Zs - C %*% S)
  # let x = Sigms0Inv %*% (dat - mu + Zs - C %*% S)
  # calculate x using Cholesky decomp of Sigma0
  genSj = function(S) {
    Zs = rnorm(n)*tau
    #vec = dat - mu97 + Zs - C %*% S
    vec = log97 - mu + Zs - S[inds]
    z = backsolve(Sigma0U, vec, transpose=TRUE)
    x = cbind(backsolve(Sigma0U, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
    #Sj = as.numeric(S + Sigma %*% Ct %*% x)
    Sj = S + Sigma10 %*% x
    return(Sj)
  }
  
  #compute conditional covariance matrix and Cholesky decomp
  UInv = backsolve(r=Sigma0U, x=diag(nrow=nrow(Sigma0U)))
  Sigma0Inv = UInv %*% t(UInv)
  SigmaCond = Sigma - Sigma10 %*% Sigma0Inv %*% Sigma01
  SigmaCondL = t(chol(SigmaCond))
  
  # for each simulation of S, generate the conditional simulation, Sj
  for(i in 1:nsim) {
    #get predictions on untransformed (non-log) scale
    preds[,i] = exp(genSj(preds[,i]) + mu)
  }
  
  # now generate test Sjs
  predsTest = matrix(rnorm(nsim*N), nrow=N)
  predsTest = SigmaCondL %*% predsTest # simulate zero mean with correct covariance struct
  predsTest = sweep(predsTest, 1, muc, "+") + mu # add muc + mu
  predsTest = exp(predsTest) # undo log transform
  
  # get 5%, 50%, and 95% quantiles
  quant5 = apply(preds, 1, quantile, probs=.05)
  quant50 = apply(preds, 1, quantile, probs=.5)
  quant95 = apply(preds, 1, quantile, probs=.95)
  
  quant5Test = apply(predsTest, 1, quantile, probs=.05)
  quant50Test = apply(predsTest, 1, quantile, probs=.5)
  quant95Test = apply(predsTest, 1, quantile, probs=.95)
  
  #get distributions for proportion larger than 3, 5, and 7 micrograms of lead
  propBigger = function(dat, thresh = 3) {
    mean(dat > thresh)
  }
  above3Prop = apply(preds, 2, propBigger, thresh=3)
  above5Prop = apply(preds, 2, propBigger, thresh=5)
  above7Prop = apply(preds, 2, propBigger, thresh=7)
  
  above3PropTest = apply(preds, 2, propBigger, thresh=3)
  above5PropTest = apply(preds, 2, propBigger, thresh=5)
  above7PropTest = apply(preds, 2, propBigger, thresh=7)
  
  return(list(simCoords=simCoords, 
              quant5=quant5, quant50=quant50, quant95=quant95, above3Prop=above3Prop, above5Prop=above5Prop, 
              above7Prop=above7Prop, 
              quant5All=quant5Test, quant50All=quant50Test, quant95All=quant95Test, above3PropAll=above3PropTest, above5PropAll=above5PropTest, 
              above7PropAll=above7PropTest))
}

# SjTest = genPredictionsTest(nsim=10000)
# makePlotsJoint(SjTest, zlim=c(0,11.5))

#
# 4) start optimizing at the exploded parameters and check likelihood of X given S
#
# 5) compute average Lambda(x) and plot each iterationx
# 
# 6) try more samples and check for standard errors of parts of likelihood.  Maybe there aren't enough samples
# DONE: for 10k iterations, didn't work :\ tau and the means still exploded and the range shrunk
#
# start:
# 
#      mu97         sigma           phi           tau          beta           lik   lik.XGivenS lik.YGivenS0j lik.S0jGivenY        lik.GP
# 1.5150000     0.3714835     0.3130000     0.2428992    -2.1980000  -204.6975649   -76.1604483   -74.1892978    34.3076407    -0.1800151 
#
# end:
#      mu97         sigma           phi           tau          beta           lik   lik.XGivenS lik.YGivenS0j lik.S0jGivenY        lik.GP 
# 2.0444096     0.2897138     0.2560072     0.7174064    -2.8024180  -137.0789725   -56.0486294   -67.1620509    22.0036604    15.4244144 
# 

# Check to see if preferential likelihood works when using all data (including 2000 data) in pref. lik.

# maximize the parameters and using Monte Carlo likelihood calculations on a grid
# the paper uses 10,000 MC samples of S given Y per likelihood evaluation.  The 1997
# data is treated as preferential, and the 2000 data is not in the joint model.  It 
# looks like the authors assume conditional independence of the two realizations 
# given the parameters???  Or maybe they assume they take place at the same time???
MLPrefJoint = function(coords97, log97, coords00, log00, 
                  initParams = c(mu97MLE, mu00MLE, log(sqrt(sillMLE)), 
                                 log(scaleMLE), log(sqrt(nuggetMLE)), beta=betaMLE), 
                  res=60, nMCSamples = 10000, doPar=FALSE, nProc=4) {
  debug=FALSE
  allCoords = rbind(coords97, coords00)
  logAll = c(log97, log00)
  
  ##### decide on initial parameters if they aren't provided.  Note that the year 
  ##### 2000 model doesn't require a beta parameter since not preferential
  if(is.null(initParams)) {
    initParamsJoint = c(mean(log97), mean(log00), log(sd(c(log00, log97))), 
                        log(.3), log(sqrt(10^-3)), 0)
    # or should initial guess for beta be related to difference in means?
  }
  else {
    initParamsJoint = initParams
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
  roundX = roundToRange(allCoords[,1], xs)
  roundY = roundToRange(allCoords[,2], ys)
  roundCoordsAll = cbind(roundX, roundY)
  
  # find indices of lattice coords corresponding to data coords
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoordsAll, latticeCoords)
  
  #generate C matrix giving indices of X in X* (indices of rounded 
  #observation locations in grid).  Also make C transpose
  n = length(logAll)
  N = nrow(latticeCoords)
  C = Matrix(0, nrow=n, ncol=N, sparse=TRUE)
  C[matrix(c(1:n, inds), ncol=2)] = 1
  Ct = t(C)
  
  # precompute distance matrices
  distMatLattice = rdist(latticeCoords)
  distMatAll = rdist(allCoords)
  
  ##### prepare log likelihood functions for maximization
  
  summTableJoint = matrix(nrow=0, ncol=12)
  colnames(summTableJoint) = c("mu97", "mu00", "sigma", "phi", "tau", "beta", "loglik", 
                               "lik.XGivenS", "lik.YGivenS0j", "lik.S0jGivenY", "lik.S0j", "MCSE")
  
  getLikJoint = function(params) {
    # Evaluate log likelihood
    out = likPreferentialMCJoint(params, latticeCoords, log97, distMatLattice, C, Ct, inds, 
                                 res=res, nsims=nMCSamples, npDat=log00, distMatData=distMatAll, doPar, nProc)
    logliks = out$allLogLiks
    loglik = logliks[1]
    MCSE = out$MCSE[1]
    
    # update summary table
    params[3:5] = exp(params[3:5])
    names(params) = colnames(summTableJoint)[1:6]
    summTableJoint <<- rbind(summTableJoint, c(params, loglik, logliks[-1], MCSE))
    tmp = c(params, loglik, logliks[-1], MCSE)
    names(tmp) = colnames(summTableJoint)
    print(tmp)
    
    if(!is.finite(loglik))
      return(log(10^-150))
    else
      return(loglik)
  }
  
  ##### maximize likelihood
  controls = list(fnscale=-1, parscale=.1)
  MLEsJoint = optim(initParamsJoint, getLikJoint, control=list(fnscale=-1, parscale=rep(.1, 6), maxit=50, 
                                                               reltol=getRelTol(v=-150, nMC=nMCSamples)), 
                    hessian=TRUE)
  MLEsJoint$par[3:5] = exp(MLEsJoint$par[3:5])
  MLEsJoint$par[c(3,5)] = MLEsJoint$par[c(3,5)]^2
  names(MLEsJoint$par) = c("mu97", "mu00", "sigmasq", "phi", "tausq", "beta")
  
  initParamsJoint[3:5] = exp(initParamsJoint[3:5])
  initParamsJoint[c(3,5)] = initParamsJoint[c(3,5)]^2
  names(initParamsJoint) = c("mu97", "mu00", "sigmasq", "phi", "tausq", "beta")
  
  return(list(summTableJoint=summTableJoint, MLEsJoint=MLEsJoint))
}

# make a function to calculate log-likelihood under preferential model.  Note that 
# this calculates the likelihood for any SINGLE dataset, not both 97 and 00 data.
# based on Eq. (9) and (10) from diggle (2010).  The input distMat should be for the
# gridCoords

#modifications for joint likelihood calculations:
# params should have mu00 at second index
# npDat should be set to log00
# C, Ct, and inds should account for increased data size
# res might need to be higher
likPreferentialMCJoint = function(params, gridCoords, dat, distMat, C, Ct, inds, 
                             res=60, nsims=10000, npDat=NULL, distMatData=NULL, doPar=FALSE, nProc=4) {
  # n and N as defined in paper: n is number of observations, N is number of pts in grid
  n97 = length(dat)
  n = n97 + length(npDat)
  N = nrow(gridCoords)
  
  #get parameters:
  mu97 = params[1]
  mu00 = params[2]
  sigma = exp(params[3])
  phi = exp(params[4])
  tau = exp(params[5])
  beta = params[6] #(note: we don't care about beta until we calculate likelihood of X given S)
  
  #generate marginal covariance matrix of S, Sigma, for the given parameters
  Sigma = Exp.cov(gridCoords, theta=phi, distMat=distMat)*sigma^2
  Sigma0 = Exp.cov(NA, theta=phi, distMat=distMatData)*sigma^2
  
  # add nugget error to Sigma0 with fast internal function of fields in C
  invisible(.Call("addToDiagC", Sigma0, as.double(rep(tau^2, nrow(Sigma0))), nrow(Sigma0)))
  
  #Now compute Sigma0, the covariance matrix of Y
  #Sigma0 = as.matrix(C %*% Sigma %*% Ct)
  Sigma10 = Sigma[,inds]
  
  #Cholesky decomposition of Sigma0 required for likelihood
  Sigma0U = chol(Sigma0)
  
  # compute muc, expectation of Sj based on this formula:
  # Sigma %*% Ct %*% Sigma0^-1 %*% (dat - mu97)
  vec = c(dat - mu97, npDat - mu00)
  z = backsolve(Sigma0U, vec, transpose=TRUE)
  x = cbind(backsolve(Sigma0U, z))
  #muc = as.numeric(Sigma %*% Ct %*% x)
  muc = Sigma10 %*% x # checked, muc correct given above formula
  
  #calculate conditional mean of S given Y (just get sample mean)
  # muc = sigma^2*(dat - mu)/((sigma + tau)*sqrt(sigma^2 + tau^2))
  
  #estimate alpha using MOM estimator (I assume this is what Diggle does?)
  alpha = log(n97/Area) - beta*mu97 - beta^2*sigma^2/2
  lambdas <<- matrix(nrow=nrow(gridCoords), ncol=nsims)
  Sjs <<- lambdas
  
  # calculate log-likelihood for given simulation, S, under preferential model
  getSimLikJoint = function(S) {
    #simulate S given Y (authors denote this by Sj).  We want:
    # Sj = S + Sigma %*% Ct %*% Sigma0Inv %*% (dat - mu + Zs - C %*% S)
    # let x = Sigms0Inv %*% (dat - mu + Zs - C %*% S)
    # calculate x using Cholesky decomp of Sigma0
    Zs = rnorm(n)*tau
    #vec = dat - mu97 + Zs - C %*% S
    vec = c(dat - mu97, npDat - mu00) + Zs - S[inds]
    z = backsolve(Sigma0U, vec, transpose=TRUE)
    x = cbind(backsolve(Sigma0U, z)) # checked, it works given above formula
    #Sj = as.numeric(S + Sigma %*% Ct %*% x)
    Sjs[,count] <<- S + Sigma10 %*% x
    S0j = Sjs[inds, count]
    
    # calculate lambda(x) given Y since needed in likelihood
    lambdas[,count] <<- exp(alpha + beta*(mu97 + Sjs[,count]))
    
    # calculate likelihood of X given S|Y simulation.  Likelihood given by Eq. (8) 
    # in Diggle (2013). There might be faster ways to calculate it though.
    lambda = mean(lambdas[,count]) # N*lambda*dA should be about n (where dA = Area/N)
    # lik.XGivenS = prod(lambdas[inds])*(lambda*Area)^(-n)
    lik.XGivenS = sum(log(lambdas[inds[1:n97],count])) - n97*log(Area*lambda) #use log-likelihood, not likelihood
    
    # calculate P(Y | S_{0j}) / P(S_{0j} | Y)
    if(tau == 0) {
      lik.frac = 0
      lik.YGivenS0j=1
      lik.S0jGivenY=1
    }
    else {
      # Calculate P(Y | S_{0j}):
      # Use iid normals since Y's conditionally indep given S
      EY = S0j
      EY[1:n97] = EY[1:n97] + mu97
      EY[(n97+1):n] = EY[(n97+1):n] + mu00
      lik.YGivenS0j = sum(log(dnorm((c(dat, npDat) - EY)/tau))) #use log-likelihood not likelihood
      
      # Get P(S_{0j} | Y):
      # First calculate SigmaS0jGivenY =? tau^2*(I - tau^2*Sigma0Inv)
      UInv = backsolve(r=Sigma0U, x=diag(nrow=nrow(Sigma0U)))
      Sigma0Inv = UInv %*% t(UInv)
      SigmaS0jGivenY = tau^2*(diag(nrow=n) - tau^2*Sigma0Inv)
      #was compared with:
      # Sigma0Tilde - Sigma0Tilde %*% (Sigma0Tilde + tau^2*I)^-1 %*% Sigma0Tilde
      
      # second calculate expectation of S0j|Y = Sigma %*% t(C) %*% Sigma0Inv %*% cbind(y)
      # we've already calculated this over all lattice points, just take right values
      muc0 = muc[inds]
      
      # then do Cholesky decomposition and solve using lik.GP
      SigmaS0jGivenYU = chol(SigmaS0jGivenY)
      lik.S0jGivenY = likGP(S0j - muc0, SigmaS0jGivenYU)
      
      # lik.frac = lik.YGivenS0j/lik.S0jGivenY
      # the ratio of likelihoods should be 1 when tausq=0, i.e. the difference in log liks should be 0
      lik.frac = lik.YGivenS0j - lik.S0jGivenY
    }
    
    # Calcualte log-likelihood of S0j marginally:
    # first calculate covariance matrix (same as Sigma0 but without nugget)
    Sigma0Tilde = Sigma0 - diag(tau^2, nrow=n) 
    Sigma0TildeU = chol(Sigma0Tilde)
    lik.GP = likGP(S0j, Sigma0TildeU)
    
    lik = lik.XGivenS + lik.frac + lik.GP
    count <<- count + 1
    
    #     if(tau > .45){
    #       print("exploded...")
    #       print("badly...")
    #       print("very badly...")
    #     }
    
    return(c(lik=lik, lik.XGivenS=lik.XGivenS, lik.YGivenS0j=lik.YGivenS0j, 
             lik.S0jGivenY=lik.S0jGivenY, lik.S0j=lik.GP))
  }
  
  # compute likelihood for simulated GP S and its antithetic pair: 2mu_c - S
  getLikPair = function(i, nPairs=1) {
    Ss = genS(params, gridCoords, nsim=nPairs)
    
    # muc = sigma^2*(dat - mu)/((sigma + tau)*sqrt(sigma^2 + tau^2))
    #     
    #     return(c(getSimLik(S), getSimLik(SAnti)))
    
    #get likelihoods given original simulations of S
    count<<- 1
    liks = apply(Ss, 2, getSimLikJoint)
    # testLiks = apply(Ss, 2, naiveSimLik)
    
    # Now generate the ``antithetic'' pairs and calculate their likelihoods: SAnti = 2*muc - S
    Ss = sweep(-Ss, 1, 2*muc, FUN="+")
    liks = cbind(liks, apply(Ss, 2, getSimLikJoint))
    # testLiks = cbind(testLiks, apply(-Ss, 2, naiveSimLik))
    return(list(liks=liks))
    # return(list(liks=liks, testLiks=testLiks))
  }
  
  #now we can compute all log-likelihoods using apply
  nPairs = ceil(nsims/2)
  
  if(!doPar) {
    # allLogLiks = c(sapply(1:nPairs, getLikPair))
    allLogLiks = getLikPair(1, nPairs=nPairs)
  }
  else {
    clust = makeCluster(nProc)
    clusterEvalQ(clust, setwd("~/git/prelim/code/"))
    clusterEvalQ(clust, source("loadAll.R"))
    allLogLiks = c(parSapply(clust, 1:nPairs, getLikPair))
    stopCluster(clust)
  }
  
  # function for converting log likelihoods to likelihoods, averaging, and logging again
  # that is numerically stable
  shiftExpFun = function(logs, fun="mean") {
    M = max(logs)
    expThresh = 500 #make sure all logs < expThresh before exponentiating
    C = M - expThresh
    
    return(C + log(do.call(fun, list(exp(logs - C)))))
  }
  shiftExpMean = function(logs) {
    shiftExpFun(logs)
  }
  
  # final MC likelihood estimate is average of likelihoods, but return log-likelihood:
  MCSD = apply(allLogLiks$liks, 1, sd)
  MCSE = MCSD/sqrt(nsims)
  #   MCSDTest = apply(allLogLiks$testLiks, 2, sd)
  #   MCSETest = MCSDTest/sqrt(nsims)
#   return(list(logLik = log(mean(exp(allLogLiks$liks[1,]))), MCSD=MCSD, MCSE=MCSE, 
#               allLogLiks=log(apply(exp(allLogLiks$liks), 1, mean))))
    return(list(logLik = shiftExpMean(allLogLiks$liks[1,]), MCSD=MCSD, MCSE=MCSE, 
                allLogLiks=apply(allLogLiks$liks, 1, shiftExpMean)))
}


# 7) assume tau=0, see if convergence is correct

# results for joint likelihood:
# start:
# mu97          mu00         sigma           phi          beta      loglik97      loglik00        loglik   lik.XGivenS lik.YGivenS0j lik.S0jGivenY 
# 1.5150000     0.7620000     0.3714835     0.3130000     -2.198000  -165.2977595   -43.2264694  -208.5242289   -82.6206642     1.0000000     1.0000000 
# lik.S0j 
# -82.6770953 
# end: 
# mu97          mu00         sigma           phi          beta      loglik97      loglik00        loglik   lik.XGivenS lik.YGivenS0j lik.S0jGivenY 
# 1.9764589     0.7130755     0.4219840     0.1306503     -1.807657  -116.2168128   -43.8593381  -160.0761509   -64.4326351     1.0000000     1.0000000 
# lik.S0j 
# -51.7841777 
# 
# notes: I can't tell if the above result shows a mistake in the likelihood calculations.  If tau = 0, then we must do perfect 
# interpolation, but sometimes the data changes on short scales.  This could be causing the large sigma and low phi.  Not sure 
# what's causing the high mu97, although it could just be drift.  mu97 changes the simulations of Sj, but otherwise doesn't 
# affect the likelihood directly in the no nugget case.  Also, lik.XGivenS should work to increase phi, but not seen here...

MLPrefNoNugget = function(coords97, log97, coords00, log00, 
                  initParams = c(mu97MLE, mu00MLE, log(sqrt(sillMLE)), 
                                 log(scaleMLE), beta=betaMLE), 
                  res=60, nMCSamples = 10000, doPar=FALSE, nProc=4) {
  
  ##### decide on initial parameters if they aren't provided.  Note that the year 
  ##### 2000 model doesn't require a beta parameter since not preferential
  if(is.null(initParams)) {
    initParamsJoint = c(mean(log97), mean(log00), log(sd(c(log00, log97))), 
                        log(.3), 0)
    initParams00 = c(mean(log00), log(sd(log00)), log(.3))
    initParams97 = c(mean(log97), log(sd(log97)), log(.3), 0)
    # or should initial guess for beta be related to difference in means?
  }
  else {
    initParamsJoint = initParams
    initParams00 = initParams[-c(1, 5)]
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
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords97, latticeCoords)
  
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
  summTable97 = matrix(nrow=0, ncol=9)
  colnames(summTable97) = c("mu97", "sigma", "phi", 
                            "beta", 
                            "loglik97", "lik.XGivenS", "lik.YGivenS0j", "lik.S0jGivenY", "lik.S0j")
  summTable00 = matrix(nrow=0, ncol=4)
  colnames(summTable00) = c("mu00", "sigma", "phi", 
                            "loglik00")
  
  getLik00 = function(params, catTable=TRUE) {
    #get parameters:
    mu = params[1]
    sigma = exp(params[2])
    phi = exp(params[3])
    
    #return very zero likelihood if parameters are negative
    if(any(c(sigma, phi) < 0)) {
      if(catTable)
        summTable00 <<- rbind(summTable00, c(params, -Inf))
      return(log(10^-150))
    }
    
    distMat = distMat00
    logDat = log00
    coords = coords00
    
    # get Cholesky decomposition of Covariance matrices
    Sigma = Exp.cov(coords, distMat=distMat, theta=phi, onlyUpper=TRUE)*sigma^2
    #invisible(.Call("addToDiagC", Sigma, as.double(rep(tau^2, nrow(Sigma))), nrow(Sigma)))
    SigmaU = chol(Sigma)
    
    # get log-likelihood assuming mean-centered data
    loglik = likGP(logDat - mu, SigmaU)
    
    # update summary table
    params[2:3] = exp(params[2:3])
    if(catTable) {
      summTable00 <<- rbind(summTable00, c(params, loglik))
      c(params, loglik)
    }
    
    if(!is.finite(loglik)) {
      return(log(10^-150))
    }
    
    return(loglik)
  }
  
  getLiks97 = function(params, catTable=TRUE) {
    # Evaluate log likelihood
    lt = params[1:3]
    rt = params[4:length(params)]
    withNuggetParams = c(lt, -Inf, rt) #exp(-Inf) = 0.  Do this since nugget is 0
    out = likPreferentialMC(withNuggetParams, latticeCoords, log97, distMat97Lattice, C, Ct, inds, 
                            res=res, nsims=nMCSamples, doPar, nProc)
    loglik = out$logLik
    # update summary table
    params[2:3] = exp(params[2:3])
    
    if(!is.finite(loglik))
      out$allLogLiks[1] = log(10^-150)
    
    if(catTable) {
      summTable97 <<- rbind(summTable97, c(params, out$allLogLiks))
      print(c(params, out$allLogLiks))
    }
    
    return(out$allLogLiks)
  }
  
  getLik97 = function(params, catTable=TRUE) {
    out = getLiks97(params, catTable=catTable)
    out[1]
  }
  
  summTableJoint = matrix(nrow=0, ncol=12)
  colnames(summTableJoint) = c("mu97", "mu00", "sigma", "phi", "beta", 
                               "loglik97", "loglik00", "loglik", 
                               "lik.XGivenS", "lik.YGivenS0j", "lik.S0jGivenY", "lik.S0j")
  
  getLikJoint = function(params) {
    # Evaluate log likelihood
    logliks97 = getLiks97(params[-2], catTable=FALSE)
    loglik00 = getLik00(params[-c(1, 5)], catTable=FALSE)
    loglik = logliks97[1] + loglik00
    
    # update summary table
    params[3:4] = exp(params[3:4])
    summTableJoint <<- rbind(summTableJoint, c(params, logliks97[1], loglik00, loglik, logliks97[-1]))
    tmp = c(params, logliks97[1], loglik00, loglik, logliks97[-1])
    names(tmp) = colnames(summTableJoint)
    print(tmp)
    
    if(!is.finite(loglik))
      return(log(10^-150))
    else
      return(loglik)
  }
  
  ##### maximize likelihood
  controls = list(fnscale=-1, parscale=.1)
  MLEs97 = optim(initParams97, getLik97, control=list(fnscale=-1, parscale=rep(.1, 4), 
                                                      reltol=getRelTol(nMC=nMCSamples), maxit=50), 
                 hessian=TRUE)
  MLEs97$par[2:3] = exp(MLEs97$par[2:3])
  MLEs97$par[2] = MLEs97$par[2]^2
  names(MLEs97$par) = c("mu97", "sigmasq", "phi", "beta")
  
  MLEs00 = optim(initParams00, getLik00, control=list(fnscale=-1, parscale=rep(.1, 3), maxit=50, 
                                                      reltol=getRelTol(nMC=nMCSamples)), 
                 hessian=TRUE)
  MLEs00$par[2:3] = exp(MLEs00$par[2:3])
  MLEs00$par[2] = MLEs00$par[2]^2
  names(MLEs00$par) = c("mu00", "sigmasq", "phi")
  
  MLEsJoint = optim(initParamsJoint, getLikJoint, control=list(fnscale=-1, parscale=rep(.1, 5), maxit=50, 
                                                               reltol=getRelTol(v=-150, nMC=nMCSamples)), 
                    hessian=TRUE)
  MLEsJoint$par[3:4] = exp(MLEsJoint$par[3:4])
  MLEsJoint$par[3] = MLEsJoint$par[3]^2
  names(MLEsJoint$par) = c("mu97", "mu00", "sigmasq", "phi", "beta")
  
  initParams97[2:3] = exp(initParams97[2:3])
  initParams97[2] = initParams97[2]^2
  names(initParams97) = c("mu97", "sigmasq", "phi", "beta")
  initParams00[2:3] = exp(initParams00[2:3])
  initParams00[2] = initParams00[2]^2
  names(initParams00) = c("mu00", "sigmasq", "phi")
  initParamsJoint[3:4] = exp(initParamsJoint[3:4])
  initParamsJoint[4] = initParamsJoint[3]^2
  names(initParamsJoint) = c("mu97", "mu00", "sigmasq", "phi", "beta")
  
  return(list(summTable97=summTable97,summTable00=summTable00, 
              summTableJoint=summTableJoint, MLEs97 = MLEs97, 
              MLEs00 = MLEs00, MLEsJoint = MLEsJoint))
}

##### Simulate data on Galicia spatial domain
sim97Data = function(npts=63, res=60, 
                     mu=mu97MLE, sigmasq=varianceMLE, phi=scaleMLE, tausq=nuggetMLE, beta=betaMLE) {
  ##### First simulate data
  #get grid over spatial domain
  xrange = win$xrange
  yrange = win$yrange
  xs = seq(xrange[1], xrange[2], l=res)
  ys = seq(yrange[1], yrange[2], l=res)
  
  # get lattice coordinates, the coordinates for the grid approximation to the GP
  win.poly = matrix(unlist(win$bdry[[1]]), ncol=2)
  latticeCoords = make.surface.grid(list(x=xs, y=ys))
  
  # only take points within domain
  latticeCoords = latticeCoords[in.poly(latticeCoords, win.poly),]
  
  # generate LGCP
  LGCP = genPreferentialPP2(numSamples=1, xyRes = res, 
                               mu=mu, sigmasq=sigmasq, phi=phi, kappa=0.5, beta=beta, 
                               npts=npts, GPcoords=latticeCoords, method="circulant", replace=FALSE)
  PP = LGCP$sims[[1]]
  GP = LGCP$GPs[[1]]$v
  coords = cbind(PP$x, PP$y)
  
  # round true coordinates to the lattice coordinates
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords[,1], xs)
  roundY = roundToRange(coords[,2], ys)
  roundCoords = cbind(roundX, roundY)
  
  # find indices of lattice coords corresponding to data coords
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords, latticeCoords)
  
  #generate observations with noise
  PPmarks = GP[inds]
  PPmarks = PPmarks + rnorm(length(PPmarks), sd = sqrt(tausq))
  
  return(list(PPcoords=coords, PPmarks=PPmarks, GPcoords=latticeCoords, GPmarks=GP))
}

sim00Data = function(GPcoords, GPmarks, res=60, 
                     mu=mu00MLE, sigmasq=varianceMLE, phi=scaleMLE, tausq=nuggetMLE) {
  
  #get grid over spatial domain
  xrange = win$xrange
  yrange = win$yrange
  xs = seq(xrange[1], xrange[2], l=res)
  ys = seq(yrange[1], yrange[2], l=res)
  
  # round true coordinates to the lattice coordinates
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords00[,1], xs)
  roundY = roundToRange(coords00[,2], ys)
  roundCoords = cbind(roundX, roundY)
  
  # find indices of lattice coords corresponding to data coords
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords, GPcoords)
  
  #get log-lead data at observation locations and add measurement noise
  PPmarks = GPmarks[inds] + rnorm(nrow(coords00), sd=sqrt(tausq))
  
  return(list(PPcoords=coords00, PPmarks=PPmarks, GPcoords=GPcoords, GPmarks=GPmarks))
}

# sim97 = sim97Data()
# quilt.plot(sim97$GPcoords, sim97$GPmarks)
# quilt.plot(sim97$PPcoords, sim97$PPmarks)
# 
# sim00 = sim00Data(sim97$GPcoords, sim97$GPmarks)
# quilt.plot(sim00$GPcoords, sim00$GPmarks)
# quilt.plot(sim00$PPcoords, sim00$PPmarks)
# 
# xs = seq(0, 1.2, l=100)
# trueVG = function(x) {
#   sqrt(varianceMLE*(1 - Matern(x, nu=.5, range=scaleMLE)) + nuggetMLE)
# }
# noNuggetVG = function(x) {
#   sqrt(varianceMLE*(1 - Matern(x, nu=.5, range=scaleMLE)))
# }
# plot(vgram(sim97$PPcoords, sim97$PPmarks), pch=19, cex=.6, xlim=c(0, 1.2), N=40)
# lines(xs, trueVG(xs))
# plot(vgram(sim00$PPcoords, sim00$PPmarks), pch=19, cex=.6, xlim=c(0, 1.2), N=40)
# lines(xs, trueVG(xs))
# 
# # simML = MLPref(sim97$PPcoords, sim97$PPmarks, sim00$PPcoords, sim00$PPmarks)
# 
# > simML = MLPref(sim97$PPcoords, sim97$PPmarks, sim00$PPcoords, sim00$PPmarks)
# Called from: eval(expr, envir, enclos)
# Browse[1]> n
# debug at ~/git/prelim/code/ML.R#409: MLEs97 = optim(initParams97, getLik97, control = list(fnscale = -1, 
# parscale = rep(0.1, 5), reltol = getRelTol(nMC = nMCSamples), 
# maxit = 50), hessian = TRUE)
# Browse[2]> MLEsJoint = optim(initParamsJoint, getLikJoint, control=list(fnscale=-1, parscale=rep(.1, 6), maxit=50, 
#                                                                         +                                                                reltol=getRelTol(v=-150, nMC=nMCSamples)), 
#                              +                     hessian=TRUE)

#start:
# mu97          mu00         sigma           phi           tau          beta      loglik97      loglik00        loglik   lik.XGivenS lik.YGivenS0j 
# 1.51500000    0.76200000    0.37148351    0.31300000    0.24289916   -2.19800000 -203.94838195  -91.79134564 -295.73972759  -82.57053178  -73.84205381 
# lik.S0jGivenY       lik.S0j          MCSE 
# 31.71047597    2.59396786    0.06148979 
#end:
# mu97          mu00         sigma           phi           tau          beta      loglik97      loglik00        loglik   lik.XGivenS lik.YGivenS0j 
# 1.54472689    1.35544492    0.24856076    0.49861147    0.52790627   -2.88850253 -159.36076023  -82.15911362 -241.51987385  -79.91149927  -70.20372666 
# lik.S0jGivenY       lik.S0j          MCSE 
# 44.14377630   39.06908877    0.07618589 
##### joint likelihood optimization:
#start:
#           mu97      mu00     sigma       phi       tau      beta  loglik97  loglik00    loglik lik.XGivenS lik.YGivenS0j lik.S0jGivenY      lik.S0j       MCSE
# [1,]  1.515000 0.7620000 0.3714835 0.3130000 0.2428992 -2.198000 -204.8851 -41.06680 -245.9519   -76.34801     -74.05363      33.47853   0.06400469 0.07983608
#end:
# [49,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434425 -2.446210 -156.0276 -59.86810 -215.8957   -61.34278     -73.69954      43.91348  34.42287030 0.06574451
#full table:
#          mu97      mu00     sigma       phi       tau      beta  loglik97  loglik00    loglik lik.XGivenS lik.YGivenS0j lik.S0jGivenY      lik.S0j       MCSE
# [1,] 1.515000 0.7620000 0.3714835 0.3130000 0.2428992 -2.198000 -204.8851 -41.06680 -245.9519   -76.34801     -74.05363      33.47853   0.06400469 0.07983608
# [2,] 1.734800 0.7620000 0.3714835 0.3130000 0.2428992 -2.198000 -201.8599 -41.06680 -242.9267   -72.42533     -73.35933      32.99966  -0.47594544 0.07791612
# [3,] 1.515000 0.9818000 0.3714835 0.3130000 0.2428992 -2.198000 -204.4281 -43.37749 -247.8055   -75.89094     -74.16559      33.55725   2.49497475 0.08060571
# [4,] 1.515000 0.7620000 0.4628044 0.3130000 0.2428992 -2.198000 -203.3898 -46.38312 -249.7730   -76.25404     -74.28351      28.20852 -10.68936449 0.09889938
# [5,] 1.515000 0.7620000 0.3714835 0.3899440 0.2428992 -2.198000 -204.8422 -39.91091 -244.7531   -75.03704     -74.70989      36.55905   5.93640440 0.08536491
# [6,] 1.515000 0.7620000 0.3714835 0.3130000 0.3026105 -2.198000 -189.3936 -46.96986 -236.3635   -76.51982     -73.18931      27.11693   3.82093252 0.08414683
# [7,] 1.515000 0.7620000 0.3714835 0.3130000 0.2428992 -1.978200 -204.2323 -41.06680 -245.2991   -75.69523     -74.24931      32.49255   2.03018287 0.07018769
# [8,] 1.588267 0.8352667 0.2981821 0.3367935 0.2613637 -2.124733 -199.5943 -41.08233 -240.6767   -73.51100     -75.58686      40.46430  14.09747664 0.06721510
# [9,] 1.569950 0.8169500 0.3328209 0.3306807 0.2566200 -2.143050 -201.1427 -40.95444 -242.0972   -75.28222     -74.50137      36.18744   7.69277315 0.07294371
# [10,] 1.612689 0.5666222 0.3452393 0.3451200 0.2678254 -2.100311 -194.5236 -42.66760 -237.1912   -72.10282     -74.35599      34.65137   8.98734741 0.07428281
# [11,] 1.588267 0.6704167 0.3516212 0.3367935 0.2613637 -2.124733 -197.9653 -41.56660 -239.5319   -74.05530     -74.23454      34.69717   6.00959843 0.07529398
# [12,] 1.645252 0.7212963 0.3369099 0.3565431 0.2766901 -2.067748 -192.6984 -42.08884 -234.7872   -72.21725     -74.33009      36.11770  10.39500401 0.07255784
# [13,] 1.710378 0.7009444 0.3208492 0.3805361 0.2953096 -2.002622 -188.5506 -43.21529 -231.7659   -71.77533     -74.89494      36.81965  14.79466120 0.06835736
# [14,] 1.710378 0.7009444 0.3208492 0.3805361 0.2953096 -2.295689 -187.7488 -43.21529 -230.9640   -70.97349     -74.97899      36.36889  15.27887815 0.08193308
# [15,] 1.808067 0.6704167 0.2981821 0.4195866 0.3256141 -2.454433 -179.2426 -45.68031 -224.9229   -67.53212     -75.29258      39.57395  22.00349831 0.08596410
# [16,] 1.808067 0.6704167 0.2981821 0.3130000 0.3256141 -2.161367 -178.7243 -46.80558 -225.5299   -67.13810     -74.11717      34.42362  15.38021632 0.06681529
# [17,] 1.734800 0.6933125 0.3150258 0.3306807 0.3026105 -2.170525 -184.9161 -44.41119 -229.3273   -69.83813     -74.32102      36.12382  13.50868828 0.07162158
# [18,] 1.612689 0.6398889 0.2771165 0.3899440 0.3590286 -2.149156 -174.3805 -49.59262 -223.9731   -70.88230     -74.74303      39.39767  27.41175778 0.06616412
# [19,] 1.551633 0.5788333 0.2393446 0.4352424 0.4364962 -2.124733 -161.9663 -59.52264 -221.4890   -69.86894     -74.35344      47.23824  41.20670553 0.05794142
# [20,] 1.747011 0.4811444 0.3208492 0.3947349 0.3958720 -2.222422 -166.3308 -57.18126 -223.5120   -68.25724     -72.76407      32.31810  19.50747929 0.08074197
# [21,] 1.707325 0.5696750 0.3150258 0.3793761 0.3568432 -2.198000 -172.7740 -50.92755 -223.7015   -68.88188     -73.60311      34.09822  19.73434347 0.07756835
# [22,] 1.767363 0.7212963 0.2704306 0.4028508 0.4400641 -2.287548 -159.4301 -59.87315 -219.3033   -66.29000     -73.02895      39.65820  30.21311824 0.06803171
# [23,] 1.844700 0.7986333 0.2393446 0.4352424 0.5640896 -2.381167 -142.7435 -77.41520 -220.1587   -59.31122     -71.05695      45.74547  39.43426711 0.06298820
# [24,] 1.949173 0.5123506 0.2260869 0.4831736 0.4412599 -2.219709 -160.1955 -60.92915 -221.1246   -62.56828     -75.61541      52.88993  42.16924213 0.05256172
# [25,] 1.840630 0.5747630 0.2559714 0.4334744 0.4015520 -2.214281 -165.6183 -55.04255 -220.6608   -65.88081     -74.75773      46.11490  33.59781461 0.06205202
# [26,] 1.797212 0.5313457 0.2427788 0.4150569 0.5009438 -2.485640 -151.4531 -69.41331 -220.8664   -63.80968     -72.41235      44.88154  37.68742367 0.06421307
# [27,] 1.775504 0.5737454 0.2603055 0.4061436 0.4389460 -2.364885 -158.7663 -60.36963 -219.1359   -65.22276     -73.28728      42.49309  33.25931552 0.06620584
# [28,] 1.688669 0.5296497 0.2496686 0.5503904 0.5018792 -2.394735 -151.9831 -68.90120 -220.8843   -66.04909     -73.03140      49.52195  43.95176613 0.07119370
# [29,] 1.718519 0.5648414 0.2610017 0.4779574 0.4504278 -2.336393 -156.7557 -61.61218 -218.3679   -65.36654     -73.51046      45.56953  37.41965370 0.06921270
# [30,] 1.658820 0.4944579 0.2388276 0.4287948 0.5592078 -2.062321 -148.8569 -78.18424 -227.0412   -67.58469     -71.05579      45.74697  40.99681142 0.05763922
# [31,] 1.770755 0.6264270 0.2820865 0.4218699 0.3727526 -2.356405 -170.2831 -51.22675 -221.5098   -67.74974     -74.53456      40.80393  29.91660834 0.07613818
# [32,] 1.727790 0.7321577 0.2126306 0.4660043 0.4508275 -2.338993 -158.3835 -60.33698 -218.7204   -66.15638     -75.44413      56.14526  49.29284287 0.05606650
# [33,] 1.732595 0.6694044 0.2356647 0.4470631 0.4364118 -2.309850 -159.5218 -58.63632 -218.1581   -65.89257     -74.69027      48.47723  40.80625924 0.06052363
# [34,] 1.691326 0.6012009 0.2277797 0.4445429 0.5046260 -2.189492 -151.9748 -68.69001 -220.6648   -66.05476     -72.96389      49.85524  44.27311947 0.05734812
# [35,] 1.711183 0.6075074 0.2402879 0.4387629 0.4678239 -2.231220 -156.1405 -63.55913 -219.6996   -66.50743     -73.41413      47.56190  40.28169472 0.05874952
# [36,] 1.963631 0.6583526 0.2688086 0.4320432 0.4410053 -2.456659 -157.6207 -59.90988 -217.5306   -61.03685     -73.39709      42.02160  29.71452871 0.07086510
# [37,] 2.169630 0.6981123 0.2848740 0.4304525 0.4432773 -2.622622 -160.4116 -60.47761 -220.8892   -59.19383     -73.04899      38.72192  22.57038789 0.08134312
# [38,] 1.715635 0.6902862 0.2554814 0.4333274 0.4945951 -2.447904 -151.3189 -67.23124 -218.5502   -64.41842     -72.22240      42.42093  36.69683422 0.06929808
# [39,] 1.746884 0.6614054 0.2556038 0.4333642 0.4694862 -2.389498 -154.4768 -63.69293 -218.1697   -64.74947     -72.84985      43.42330  35.62700012 0.06733956
# [40,] 1.856982 0.6755077 0.2778203 0.4263414 0.4250331 -2.483724 -160.5678 -57.82229 -218.3901   -64.37493     -73.36905      39.91952  28.80274269 0.07761087
# [41,] 1.820532 0.6585077 0.2679205 0.4294134 0.4353491 -2.420598 -159.8663 -59.11148 -218.9778   -65.38405     -73.50222      41.44586  32.03222346 0.07219035
# [42,] 1.830675 0.5464560 0.2490822 0.4732032 0.4466405 -2.492789 -157.0859 -61.21611 -218.3020   -63.52538     -74.13436      47.19095  39.55187818 0.06867304
# [43,] 1.814847 0.5901661 0.2542559 0.4545396 0.4449873 -2.441478 -156.3798 -60.61150 -216.9913   -62.95054     -73.77566      45.30484  37.53908036 0.06851838
# [44,] 1.835649 0.6994805 0.2567504 0.4873147 0.4498231 -2.440983 -156.8968 -60.51151 -217.4083   -63.84268     -73.90383      46.73002  38.20193361 0.07052395
# [45,] 1.820613 0.6680467 0.2576346 0.4656151 0.4470788 -2.421958 -156.5481 -60.33139 -216.8795   -63.37383     -73.70720      44.44694  36.78105926 0.06934912
# [46,] 1.742048 0.5952312 0.2345906 0.4780664 0.4724443 -2.301555 -155.3111 -64.09080 -219.4019   -65.74514     -73.85608      51.23848  43.68225792 0.05964125
# [47,] 1.828248 0.6554386 0.2663181 0.4387228 0.4364201 -2.438182 -158.8882 -59.18759 -218.0757   -64.40778     -73.53399      43.24471  33.38021290 0.07267698
# [48,] 1.917088 0.7360965 0.2513913 0.4144353 0.4411382 -2.482816 -157.6907 -59.59414 -217.2848   -61.44631     -73.76490      43.12847  32.89290084 0.06507756
# [49,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434425 -2.446210 -156.0276 -59.86810 -215.8957   -61.34278     -73.69954      43.91348  34.42287030 0.06574451
# [50,] 1.928909 0.6501583 0.2560816 0.4557101 0.4152553 -2.448615 -161.8213 -56.03439 -217.8557   -62.34612     -74.69477      46.63222  35.55746992 0.06924155
# [51,] 1.883403 0.6529701 0.2559620 0.4500178 0.4281955 -2.433836 -158.7679 -57.81652 -216.5844   -62.00596     -74.29933      44.53007  35.00488421 0.06902763
# [52,] 1.868545 0.6932827 0.2537602 0.4294770 0.4434425 -2.446210 -158.0561 -59.86810 -217.9242   -63.34897     -73.71780      44.53135  36.24319830 0.06631801
# [53,] 1.868345 0.6932827 0.2537602 0.4294770 0.4434425 -2.446210 -158.5313 -59.86810 -218.3994   -63.82826     -73.65549      44.00111  33.62113769 0.06677081
# [54,] 1.868445 0.6933827 0.2537602 0.4294770 0.4434425 -2.446210 -157.5276 -59.86789 -217.3955   -62.82249     -73.70324      44.05141  34.01190374 0.06645268
# [55,] 1.868445 0.6931827 0.2537602 0.4294770 0.4434425 -2.446210 -158.1430 -59.86831 -218.0113   -63.43789     -73.63313      44.45348  35.47190575 0.06674804
# [56,] 1.868445 0.6932827 0.2537856 0.4294770 0.4434425 -2.446210 -157.5067 -59.86852 -217.3753   -62.80238     -73.70516      44.32777  34.22838516 0.06653859
# [57,] 1.868445 0.6932827 0.2537348 0.4294770 0.4434425 -2.446210 -158.3987 -59.86768 -218.2664   -63.69294     -73.68763      44.35462  34.79111415 0.06642062
# [58,] 1.868445 0.6932827 0.2537602 0.4295199 0.4434425 -2.446210 -158.0206 -59.86784 -217.8884   -63.31573     -73.64371      43.91732  35.44780119 0.06581615
# [59,] 1.868445 0.6932827 0.2537602 0.4294340 0.4434425 -2.446210 -157.6868 -59.86835 -217.5551   -62.98146     -73.61410      44.52946  35.23951907 0.06650822
# [60,] 1.868445 0.6932827 0.2537602 0.4294770 0.4434868 -2.446210 -158.5847 -59.87435 -218.4590   -63.88494     -73.71223      43.71167  35.29456049 0.06590105
# [61,] 1.868445 0.6932827 0.2537602 0.4294770 0.4433981 -2.446210 -157.4910 -59.86184 -217.3528   -62.78052     -73.74904      43.80429  34.78797052 0.06687304
# [62,] 1.868445 0.6932827 0.2537602 0.4294770 0.4434425 -2.446110 -157.6928 -59.86810 -217.5609   -62.98773     -73.71277      43.69479  34.83342745 0.06699726
# [63,] 1.868445 0.6932827 0.2537602 0.4294770 0.4434425 -2.446310 -157.3358 -59.86810 -217.2039   -62.63069     -73.73985      44.33648  34.06371571 0.06631029
# [64,] 1.866545 0.6932827 0.2537602 0.4294770 0.4434425 -2.446210 -157.7859 -59.86810 -217.6540   -63.11926     -73.73507      44.02767  34.14553138 0.06700384
# [65,] 1.866345 0.6932827 0.2537602 0.4294770 0.4434425 -2.446210 -154.0823 -59.86810 -213.9504   -59.41973     -73.80035      43.86410  34.26258928 0.06607507
# [66,] 1.866445 0.6933827 0.2537602 0.4294770 0.4434425 -2.446210 -156.6521 -59.86789 -216.5199   -61.98744     -73.72083      44.74403  34.52743640 0.06677966
# [67,] 1.866445 0.6931827 0.2537602 0.4294770 0.4434425 -2.446210 -157.6028 -59.86831 -217.4711   -62.93814     -73.71931      44.25083  34.90380297 0.06619985
# [68,] 1.866445 0.6932827 0.2537856 0.4294770 0.4434425 -2.446210 -157.8959 -59.86852 -217.7645   -63.23203     -73.69163      44.49873  34.51009614 0.06534026
# [69,] 1.866445 0.6932827 0.2537348 0.4294770 0.4434425 -2.446210 -158.0929 -59.86768 -217.9606   -63.42755     -73.69947      43.60310  34.26297016 0.06633872
# [70,] 1.866445 0.6932827 0.2537602 0.4295199 0.4434425 -2.446210 -157.6099 -59.86784 -217.4777   -62.94548     -73.59460      43.77353  34.05649158 0.06728636
# [71,] 1.866445 0.6932827 0.2537602 0.4294340 0.4434425 -2.446210 -157.9481 -59.86835 -217.8165   -63.28330     -73.71509      43.61715  34.97051681 0.06719502
# [72,] 1.866445 0.6932827 0.2537602 0.4294770 0.4434868 -2.446210 -158.2294 -59.87435 -218.1037   -63.57011     -73.72156      43.92306  34.68460321 0.06662024
# [73,] 1.866445 0.6932827 0.2537602 0.4294770 0.4433981 -2.446210 -158.1017 -59.86184 -217.9635   -63.43169     -73.75837      43.83911  34.68277423 0.06665998
# [74,] 1.866445 0.6932827 0.2537602 0.4294770 0.4434425 -2.446110 -158.3232 -59.86810 -218.1913   -63.65853     -73.72126      44.47055  34.06383429 0.06617335
# [75,] 1.866445 0.6932827 0.2537602 0.4294770 0.4434425 -2.446310 -157.9582 -59.86810 -217.8263   -63.29357     -73.71584      43.62588  33.95452769 0.06707308
# [76,] 1.867545 0.6942827 0.2537602 0.4294770 0.4434425 -2.446210 -158.3385 -59.86605 -218.2045   -63.65163     -73.71144      44.72850  34.43862922 0.06670777
# [77,] 1.867345 0.6942827 0.2537602 0.4294770 0.4434425 -2.446210 -157.6922 -59.86605 -217.5582   -63.00935     -73.66917      44.10883  34.44993236 0.06646819
# [78,] 1.867445 0.6943827 0.2537602 0.4294770 0.4434425 -2.446210 -158.4976 -59.86585 -218.3635   -63.81282     -73.58268      43.88155  35.57958879 0.06747276
# [79,] 1.867445 0.6941827 0.2537602 0.4294770 0.4434425 -2.446210 -156.6444 -59.86625 -216.5107   -61.95961     -73.71626      43.72071  34.67286449 0.06690413
# [80,] 1.867445 0.6942827 0.2537856 0.4294770 0.4434425 -2.446210 -158.1578 -59.86647 -218.0243   -63.47374     -73.68067      43.98109  35.65355113 0.06645165
# [81,] 1.867445 0.6942827 0.2537348 0.4294770 0.4434425 -2.446210 -158.0539 -59.86563 -217.9196   -63.36841     -73.83432      45.09155  34.59775275 0.06663258
# [82,] 1.867445 0.6942827 0.2537602 0.4295199 0.4434425 -2.446210 -157.7689 -59.86579 -217.6347   -63.08431     -73.67971      45.03534  34.88751644 0.06610597
# [83,] 1.867445 0.6942827 0.2537602 0.4294340 0.4434425 -2.446210 -157.5574 -59.86630 -217.4237   -62.87239     -73.70548      44.20156  34.16768495 0.06603597
# [84,] 1.867445 0.6942827 0.2537602 0.4294770 0.4434868 -2.446210 -158.3364 -59.87230 -218.2087   -63.65691     -73.68205      43.46013  33.73409904 0.06575809
# [85,] 1.867445 0.6942827 0.2537602 0.4294770 0.4433981 -2.446210 -157.6093 -59.85979 -217.4691   -62.91916     -73.70751      44.05130  35.60110395 0.06605132
# [86,] 1.867445 0.6942827 0.2537602 0.4294770 0.4434425 -2.446110 -158.7659 -59.86605 -218.6320   -64.08113     -73.73159      44.27832  35.59215689 0.06566936
# [87,] 1.867445 0.6942827 0.2537602 0.4294770 0.4434425 -2.446310 -158.8629 -59.86605 -218.7289   -64.17803     -73.78261      43.84669  34.50589908 0.06532855
# [88,] 1.867545 0.6922827 0.2537602 0.4294770 0.4434425 -2.446210 -158.3248 -59.87023 -218.1950   -63.63795     -73.66101      43.40896  34.46826176 0.06676917
# [89,] 1.867345 0.6922827 0.2537602 0.4294770 0.4434425 -2.446210 -158.6702 -59.87023 -218.5405   -63.98745     -73.67571      44.16682  34.36829045 0.06643337
# [90,] 1.867445 0.6923827 0.2537602 0.4294770 0.4434425 -2.446210 -157.2118 -59.87001 -217.0818   -62.52697     -73.67153      44.32198  34.17257298 0.06588179
# [91,] 1.867445 0.6921827 0.2537602 0.4294770 0.4434425 -2.446210 -158.6739 -59.87045 -218.5444   -63.98911     -73.65043      44.37620  34.82400706 0.06645223
# [92,] 1.867445 0.6922827 0.2537856 0.4294770 0.4434425 -2.446210 -158.0765 -59.87065 -217.9471   -63.39235     -73.64516      45.03101  33.99650176 0.06642953
# [93,] 1.867445 0.6922827 0.2537348 0.4294770 0.4434425 -2.446210 -156.0732 -59.86981 -215.9430   -61.38768     -73.73553      43.93530  33.97792224 0.06659996
# [94,] 1.867445 0.6922827 0.2537602 0.4295199 0.4434425 -2.446210 -158.1264 -59.86998 -217.9964   -63.44185     -73.66055      44.24248  33.59436752 0.06549599
# [95,] 1.867445 0.6922827 0.2537602 0.4294340 0.4434425 -2.446210 -157.9198 -59.87049 -217.7903   -63.23479     -73.68364      43.54072  34.52104150 0.06643387
# [96,] 1.867445 0.6922827 0.2537602 0.4294770 0.4434868 -2.446210 -157.5154 -59.87649 -217.3919   -62.83591     -73.71728      43.95032  34.63867968 0.06664731
# [97,] 1.867445 0.6922827 0.2537602 0.4294770 0.4433981 -2.446210 -158.2010 -59.86398 -218.0650   -63.51085     -73.74480      44.05519  34.98832907 0.06635119
# [98,] 1.867445 0.6922827 0.2537602 0.4294770 0.4434425 -2.446110 -158.5208 -59.87023 -218.3911   -63.83600     -73.64679      45.02549  34.21059753 0.06654401
# [99,] 1.867445 0.6922827 0.2537602 0.4294770 0.4434425 -2.446310 -157.7031 -59.87023 -217.5734   -63.01832     -73.73331      44.05136  34.36069632 0.06684181
# [100,] 1.867545 0.6932827 0.2540141 0.4294770 0.4434425 -2.446210 -158.2447 -59.87231 -218.1170   -63.56507     -73.73253      44.89975  34.62119376 0.06707432
# [101,] 1.867345 0.6932827 0.2540141 0.4294770 0.4434425 -2.446210 -157.8272 -59.87231 -217.6995   -63.15161     -73.71047      43.66048  34.09490749 0.06760757
# [102,] 1.867445 0.6933827 0.2540141 0.4294770 0.4434425 -2.446210 -157.7414 -59.87210 -217.6135   -63.06373     -73.71924      43.89587  33.64393837 0.06756735
# [103,] 1.867445 0.6931827 0.2540141 0.4294770 0.4434425 -2.446210 -158.5020 -59.87252 -218.3745   -63.82440     -73.65550      44.21724  34.52746937 0.06573907
# [104,] 1.867445 0.6932827 0.2540395 0.4294770 0.4434425 -2.446210 -157.1108 -59.87273 -216.9835   -62.43387     -73.66864      44.29809  34.29145442 0.06720218
# [105,] 1.867445 0.6932827 0.2539887 0.4294770 0.4434425 -2.446210 -158.3858 -59.87189 -218.2577   -63.70748     -73.72889      44.34282  34.07409470 0.06600161
# [106,] 1.867445 0.6932827 0.2540141 0.4295199 0.4434425 -2.446210 -156.3285 -59.87205 -216.2005   -61.65108     -73.67801      45.52203  34.87236460 0.06769299
# [107,] 1.867445 0.6932827 0.2540141 0.4294340 0.4434425 -2.446210 -158.5241 -59.87256 -218.3967   -63.84629     -73.62402      44.05072  35.10208568 0.06643690
# [108,] 1.867445 0.6932827 0.2540141 0.4294770 0.4434868 -2.446210 -158.3815 -59.87856 -218.2600   -63.70922     -73.63530      44.82079  34.14262254 0.06655004
# [109,] 1.867445 0.6932827 0.2540141 0.4294770 0.4433981 -2.446210 -157.3745 -59.86606 -217.2405   -62.69151     -73.63640      44.79349  34.92406938 0.06861583
# [110,] 1.867445 0.6932827 0.2540141 0.4294770 0.4434425 -2.446110 -157.6773 -59.87231 -217.5497   -62.99973     -73.62167      43.55415  34.38907109 0.06602559
# [111,] 1.867445 0.6932827 0.2540141 0.4294770 0.4434425 -2.446310 -157.9326 -59.87231 -217.8049   -63.25502     -73.66642      43.89637  34.43884945 0.06723679
# [112,] 1.867545 0.6932827 0.2535066 0.4294770 0.4434425 -2.446210 -158.4824 -59.86391 -218.3464   -63.78838     -73.72272      43.69333  34.21881946 0.06621123
# [113,] 1.867345 0.6932827 0.2535066 0.4294770 0.4434425 -2.446210 -157.0738 -59.86391 -216.9377   -62.38378     -73.69004      44.60883  34.50688908 0.06671306
# [114,] 1.867445 0.6933827 0.2535066 0.4294770 0.4434425 -2.446210 -157.6214 -59.86370 -217.4850   -62.92931     -73.66962      43.44496  35.01701734 0.06722375
# [115,] 1.867445 0.6931827 0.2535066 0.4294770 0.4434425 -2.446210 -158.5179 -59.86412 -218.3820   -63.82588     -73.69868      43.46724  34.54567489 0.06679918
# [116,] 1.867445 0.6932827 0.2535319 0.4294770 0.4434425 -2.446210 -158.3132 -59.86432 -218.1776   -63.62191     -73.68262      43.88193  34.96326633 0.06758483
# [117,] 1.867445 0.6932827 0.2534812 0.4294770 0.4434425 -2.446210 -157.4850 -59.86349 -217.3485   -62.79221     -73.69276      43.85328  34.52420946 0.06629184
# [118,] 1.867445 0.6932827 0.2535066 0.4295199 0.4434425 -2.446210 -158.6768 -59.86365 -218.5405   -63.98499     -73.62910      44.07962  34.37541650 0.06606114
# [119,] 1.867445 0.6932827 0.2535066 0.4294340 0.4434425 -2.446210 -156.7835 -59.86416 -216.6476   -62.09122     -73.73646      44.13358  33.95946795 0.06690936
# [120,] 1.867445 0.6932827 0.2535066 0.4294770 0.4434868 -2.446210 -158.8001 -59.87016 -218.6702   -64.11339     -73.73888      44.48413  35.18871187 0.06583606
# [121,] 1.867445 0.6932827 0.2535066 0.4294770 0.4433981 -2.446210 -157.8012 -59.85765 -217.6589   -63.10384     -73.70806      44.22691  33.87579596 0.06661096
# [122,] 1.867445 0.6932827 0.2535066 0.4294770 0.4434425 -2.446110 -156.5196 -59.86391 -216.3835   -61.82760     -73.71847      44.22597  34.95132271 0.06704417
# [123,] 1.867445 0.6932827 0.2535066 0.4294770 0.4434425 -2.446310 -157.4181 -59.86391 -217.2820   -62.72610     -73.77755      43.43546  34.82078225 0.06748092
# [124,] 1.867545 0.6932827 0.2537602 0.4299067 0.4434425 -2.446210 -157.5374 -59.86556 -217.4029   -62.85277     -73.68830      44.45031  33.69791201 0.06667469
# [125,] 1.867345 0.6932827 0.2537602 0.4299067 0.4434425 -2.446210 -158.7334 -59.86556 -218.5990   -64.05287     -73.65416      44.04840  34.90737702 0.06704746
# [126,] 1.867445 0.6933827 0.2537602 0.4299067 0.4434425 -2.446210 -156.9932 -59.86535 -216.8586   -62.31064     -73.75988      44.26960  34.31427635 0.06621644
# [127,] 1.867445 0.6931827 0.2537602 0.4299067 0.4434425 -2.446210 -157.8182 -59.86577 -217.6839   -63.13559     -73.73863      43.91419  34.41775157 0.06628491
# [128,] 1.867445 0.6932827 0.2537856 0.4299067 0.4434425 -2.446210 -158.2288 -59.86598 -218.0947   -63.54690     -73.70542      44.07129  34.64096357 0.06670596
# [129,] 1.867445 0.6932827 0.2537348 0.4299067 0.4434425 -2.446210 -156.5529 -59.86514 -216.4180   -61.86960     -73.68776      44.42739  34.70869247 0.06647710
# [130,] 1.867445 0.6932827 0.2537602 0.4299497 0.4434425 -2.446210 -158.2558 -59.86530 -218.1211   -63.57339     -73.73021      44.10342  34.32194806 0.06688042
# [131,] 1.867445 0.6932827 0.2537602 0.4298637 0.4434425 -2.446210 -157.1076 -59.86581 -216.9734   -62.42482     -73.69931      43.85295  34.12161875 0.06658289
# [132,] 1.867445 0.6932827 0.2537602 0.4299067 0.4434868 -2.446210 -158.6356 -59.87181 -218.5074   -63.95835     -73.68130      43.36199  34.58943181 0.06616541
# [133,] 1.867445 0.6932827 0.2537602 0.4299067 0.4433981 -2.446210 -157.4786 -59.85931 -217.3379   -62.79064     -73.68922      44.00757  33.82536369 0.06736148
# [134,] 1.867445 0.6932827 0.2537602 0.4299067 0.4434425 -2.446110 -157.7416 -59.86556 -217.6072   -63.05906     -73.76096      44.69924  36.01139048 0.06643377
# [135,] 1.867445 0.6932827 0.2537602 0.4299067 0.4434425 -2.446310 -157.9893 -59.86556 -217.8548   -63.30668     -73.71944      43.80502  33.69591423 0.06747832
# [136,] 1.867545 0.6932827 0.2537602 0.4290477 0.4434425 -2.446210 -154.9152 -59.87064 -214.7858   -60.22612     -73.65793      44.11705  34.64860917 0.06584598
# [137,] 1.867345 0.6932827 0.2537602 0.4290477 0.4434425 -2.446210 -158.0359 -59.87064 -217.9065   -63.35083     -73.69229      44.28832  34.83094593 0.06653386
# [138,] 1.867445 0.6933827 0.2537602 0.4290477 0.4434425 -2.446210 -158.2093 -59.87043 -218.0797   -63.52222     -73.69385      44.44415  33.89943166 0.06599744
# [139,] 1.867445 0.6931827 0.2537602 0.4290477 0.4434425 -2.446210 -157.6491 -59.87085 -217.5199   -62.96200     -73.74908      44.07338  34.47742940 0.06674973
# [140,] 1.867445 0.6932827 0.2537856 0.4290477 0.4434425 -2.446210 -157.8062 -59.87106 -217.6773   -63.11986     -73.68713      44.13086  34.55179883 0.06642702
# [141,] 1.867445 0.6932827 0.2537348 0.4290477 0.4434425 -2.446210 -157.5232 -59.87022 -217.3934   -62.83538     -73.68327      43.91469  33.98075210 0.06652844
# [142,] 1.867445 0.6932827 0.2537602 0.4290906 0.4434425 -2.446210 -157.7039 -59.87038 -217.5743   -63.01703     -73.69321      44.27444  34.92524170 0.06740522
# [143,] 1.867445 0.6932827 0.2537602 0.4290048 0.4434425 -2.446210 -155.2573 -59.87089 -215.1282   -60.56998     -73.71613      43.49838  34.84370637 0.06676765
# [144,] 1.867445 0.6932827 0.2537602 0.4290477 0.4434868 -2.446210 -157.7032 -59.87689 -217.5801   -63.02151     -73.70957      43.93166  34.86142731 0.06638687
# [145,] 1.867445 0.6932827 0.2537602 0.4290477 0.4433981 -2.446210 -158.4621 -59.86438 -218.3265   -63.76967     -73.71280      44.76223  34.39806387 0.06690519
# [146,] 1.867445 0.6932827 0.2537602 0.4290477 0.4434425 -2.446110 -158.2453 -59.87064 -218.1160   -63.55828     -73.70494      43.86699  33.80588505 0.06748308
# [147,] 1.867445 0.6932827 0.2537602 0.4290477 0.4434425 -2.446310 -158.3175 -59.87064 -218.1881   -63.63045     -73.56577      44.25341  34.32335910 0.06887623
# [148,] 1.867545 0.6932827 0.2537602 0.4294770 0.4438861 -2.446210 -158.3449 -59.93067 -218.2756   -63.71170     -73.70010      44.04185  33.97422094 0.06594813
# [149,] 1.867345 0.6932827 0.2537602 0.4294770 0.4438861 -2.446210 -157.4598 -59.93067 -217.3905   -62.83059     -73.73168      44.82791  34.32277284 0.06779386
# [150,] 1.867445 0.6933827 0.2537602 0.4294770 0.4438861 -2.446210 -158.5390 -59.93047 -218.4694   -63.90775     -73.70564      44.70258  34.02315379 0.06747224
# [151,] 1.867445 0.6931827 0.2537602 0.4294770 0.4438861 -2.446210 -156.2958 -59.93088 -216.2267   -61.66460     -73.75140      44.22825  33.75302045 0.06623901
# [152,] 1.867445 0.6932827 0.2537856 0.4294770 0.4438861 -2.446210 -154.1228 -59.93109 -214.0539   -59.49234     -73.57420      44.81331  36.81788150 0.06760822
# [153,] 1.867445 0.6932827 0.2537348 0.4294770 0.4438861 -2.446210 -156.2537 -59.93025 -216.1839   -61.62175     -73.70494      44.66828  34.36252594 0.06822252
# [154,] 1.867445 0.6932827 0.2537602 0.4295199 0.4438861 -2.446210 -158.3302 -59.93042 -218.2606   -63.69917     -73.65840      43.18987  34.24611824 0.06558533
# [155,] 1.867445 0.6932827 0.2537602 0.4294340 0.4438861 -2.446210 -158.1713 -59.93093 -218.1022   -63.53986     -73.62148      44.18534  34.46701100 0.06590023
# [156,] 1.867445 0.6932827 0.2537602 0.4294770 0.4439305 -2.446210 -157.1772 -59.93694 -217.1142   -62.55136     -73.60646      45.26966  34.74580179 0.06653021
# [157,] 1.867445 0.6932827 0.2537602 0.4294770 0.4438417 -2.446210 -157.9395 -59.92441 -217.8640   -63.30297     -73.64241      44.64164  35.38833737 0.06579952
# [158,] 1.867445 0.6932827 0.2537602 0.4294770 0.4438861 -2.446110 -156.0294 -59.93067 -215.9601   -61.39816     -73.63305      44.46897  34.82265516 0.06774509
# [159,] 1.867445 0.6932827 0.2537602 0.4294770 0.4438861 -2.446310 -158.3535 -59.93067 -218.2841   -63.72224     -73.66215      44.46215  33.68710417 0.06725194
# [160,] 1.867545 0.6932827 0.2537602 0.4294770 0.4429992 -2.446210 -156.6152 -59.80561 -216.4208   -61.87468     -73.72167      44.48552  35.36882296 0.06703061
# [161,] 1.867345 0.6932827 0.2537602 0.4294770 0.4429992 -2.446210 -154.9089 -59.80561 -214.7145   -60.17242     -73.70833      44.15009  34.25395494 0.06640306
# [162,] 1.867445 0.6933827 0.2537602 0.4294770 0.4429992 -2.446210 -157.6005 -59.80540 -217.4059   -62.86204     -73.69873      44.19348  34.96030484 0.06610378
# [163,] 1.867445 0.6931827 0.2537602 0.4294770 0.4429992 -2.446210 -157.4258 -59.80582 -217.2316   -62.68729     -73.73404      45.04133  34.10137513 0.06689759
# [164,] 1.867445 0.6932827 0.2537856 0.4294770 0.4429992 -2.446210 -158.5591 -59.80603 -218.3652   -63.82135     -73.76840      44.88260  34.92792314 0.06691754
# [165,] 1.867445 0.6932827 0.2537348 0.4294770 0.4429992 -2.446210 -157.9960 -59.80519 -217.8011   -63.25674     -73.64895      43.37765  33.55001755 0.06610959
# [166,] 1.867445 0.6932827 0.2537602 0.4295199 0.4429992 -2.446210 -157.1410 -59.80536 -216.9464   -62.40276     -73.64420      44.20673  34.03189745 0.06656346
# [167,] 1.867445 0.6932827 0.2537602 0.4294340 0.4429992 -2.446210 -158.3663 -59.80586 -218.1722   -63.62760     -73.68471      44.99432  33.73874880 0.06576721
# [168,] 1.867445 0.6932827 0.2537602 0.4294770 0.4430435 -2.446210 -158.1178 -59.81185 -217.9296   -63.38467     -73.72972      44.90027  34.95106161 0.06670671
# [169,] 1.867445 0.6932827 0.2537602 0.4294770 0.4429549 -2.446210 -157.0394 -59.79937 -216.8388   -62.29552     -73.66473      44.62973  35.20033739 0.06668026
# [170,] 1.867445 0.6932827 0.2537602 0.4294770 0.4429992 -2.446110 -157.2647 -59.80561 -217.0703   -62.52623     -73.77914      43.89965  34.29468136 0.06626930
# [171,] 1.867445 0.6932827 0.2537602 0.4294770 0.4429992 -2.446310 -158.3005 -59.80561 -218.1061   -63.56197     -73.70210      44.23639  33.91628662 0.06559807
# [172,] 1.867545 0.6932827 0.2537602 0.4294770 0.4434425 -2.445210 -156.5605 -59.86810 -216.4286   -61.87369     -73.71022      44.25239  34.76812296 0.06607482
# [173,] 1.867345 0.6932827 0.2537602 0.4294770 0.4434425 -2.445210 -157.1275 -59.86810 -216.9956   -62.44469     -73.68171      44.51308  35.60671751 0.06629825
# [174,] 1.867445 0.6933827 0.2537602 0.4294770 0.4434425 -2.445210 -158.1273 -59.86789 -217.9952   -63.44248     -73.65974      44.94456  34.42342466 0.06522092
# [175,] 1.867445 0.6931827 0.2537602 0.4294770 0.4434425 -2.445210 -158.3247 -59.86831 -218.1930   -63.63991     -73.72415      44.17607  33.90458562 0.06722778
# [176,] 1.867445 0.6932827 0.2537856 0.4294770 0.4434425 -2.445210 -157.0708 -59.86852 -216.9393   -62.38671     -73.68379      44.25958  34.10477303 0.06716771
# [177,] 1.867445 0.6932827 0.2537348 0.4294770 0.4434425 -2.445210 -158.4370 -59.86768 -218.3047   -63.75149     -73.69197      43.52892  33.70532745 0.06583431
# [178,] 1.867445 0.6932827 0.2537602 0.4295199 0.4434425 -2.445210 -158.4500 -59.86784 -218.3178   -63.76541     -73.68889      44.01015  34.34415485 0.06677887
# [179,] 1.867445 0.6932827 0.2537602 0.4294340 0.4434425 -2.445210 -157.9127 -59.86835 -217.7810   -63.22765     -73.68150      44.67663  34.59028406 0.06686218
# [180,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434868 -2.445210 -157.3846 -59.87435 -217.2590   -62.70515     -73.64318      44.40916  33.96385566 0.06722510
# [181,] 1.867445 0.6932827 0.2537602 0.4294770 0.4433981 -2.445210 -157.2465 -59.86184 -217.1084   -62.55634     -73.62597      44.22546  34.65127687 0.06547805
# [182,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434425 -2.445110 -158.8594 -59.86810 -218.7275   -64.17462     -73.72808      45.31239  34.61985194 0.06648294
# [183,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434425 -2.445310 -157.7924 -59.86810 -217.6605   -63.10759     -73.67348      43.95102  34.17163540 0.06574631
# [184,] 1.867545 0.6932827 0.2537602 0.4294770 0.4434425 -2.447210 -158.1453 -59.86810 -218.0134   -63.45843     -73.72301      44.95870  34.16468634 0.06666234
# [185,] 1.867345 0.6932827 0.2537602 0.4294770 0.4434425 -2.447210 -158.1508 -59.86810 -218.0189   -63.46802     -73.62597      44.64274  35.15677749 0.06667064
# [186,] 1.867445 0.6933827 0.2537602 0.4294770 0.4434425 -2.447210 -156.0598 -59.86789 -215.9277   -61.37499     -73.62469      43.63342  34.18551054 0.06726101
# [187,] 1.867445 0.6931827 0.2537602 0.4294770 0.4434425 -2.447210 -157.7500 -59.86831 -217.6183   -63.06514     -73.75473      43.91294  34.95987762 0.06739422
# [188,] 1.867445 0.6932827 0.2537856 0.4294770 0.4434425 -2.447210 -154.4216 -59.86852 -214.2901   -59.73747     -73.65077      44.60005  34.39195186 0.06676988
# [189,] 1.867445 0.6932827 0.2537348 0.4294770 0.4434425 -2.447210 -158.1678 -59.86768 -218.0354   -63.48222     -73.69195      44.15259  34.28383649 0.06690796
# [190,] 1.867445 0.6932827 0.2537602 0.4295199 0.4434425 -2.447210 -157.1835 -59.86784 -217.0514   -62.49893     -73.66867      44.22460  35.68901361 0.06729737
# [191,] 1.867445 0.6932827 0.2537602 0.4294340 0.4434425 -2.447210 -157.1399 -59.86835 -217.0083   -62.45490     -73.65881      44.36594  34.39369862 0.06650070
# [192,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434868 -2.447210 -157.8681 -59.87435 -217.7425   -63.18867     -73.69508      44.06829  37.35156508 0.06705644
# [193,] 1.867445 0.6932827 0.2537602 0.4294770 0.4433981 -2.447210 -157.4695 -59.86184 -217.3313   -62.77929     -73.71583      44.61704  34.42429851 0.06678614
# [194,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434425 -2.447110 -157.2669 -59.86810 -217.1350   -62.58206     -73.77160      43.68123  34.23486175 0.06708201
# [195,] 1.867445 0.6932827 0.2537602 0.4294770 0.4434425 -2.447310 -155.8695 -59.86810 -215.7376   -61.18467     -73.73585      44.33113  33.93908346 0.06701758