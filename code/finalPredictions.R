# Functions for generating the final predictions of the paper: Fig. 7-9
# NOTE: the figures are only for 1997
predRes <<- 64

genPredictions = function(nsim=10000, nx=predRes, ny=predRes, mu=mu97MLE, sigma=sqrt(varianceMLE), 
                          phi=scaleMLE, tau=sqrt(nuggetMLE), nu=.5) {
  # generate the simulation coordinates
  xrange = seq(win$xrange[1], win$xrange[2], l=nx)
  yrange = seq(win$yrange[1], win$yrange[2], l=ny)
  simCoords = make.surface.grid(list(x=xrange, y=yrange))
  
  # simulate S from marginal distribution (nu is assumed to be 1/2 in Matern covariance)
  # use intrinsic method instead of cutoff embedding due to rectangular grid (faster)
  # to save space, will overwrite Ss with Sjs then predictions
  params = c(mu, sigma, phi, tau)
  preds = genS(params, gridCoords=simCoords, nsim=nsim, intrinsic = TRUE)
  
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
  
  #generate C matrix giving indices of X in X* (indices of rounded 
  #observation locations in grid).  Also make C transpose
  C = Matrix(0, nrow=n, ncol=N, sparse=TRUE)
  C[matrix(c(1:n, inds), ncol=2)] = 1
  Ct = t(C)
  
  #generate marginal covariance matrix of S, Sigma, for the given parameters
  distMat = rdist(simCoords)
  Sigma = Exp.cov(simCoords, theta=phi, distMat=distMat)*sigma^2
  
  #Now compute Sigma0, the covariance matrix of Y
  #Sigma0 = as.matrix(C %*% Sigma %*% Ct)
  Sigma0 = Sigma[inds, inds]
  Sigma10 = Sigma[,inds]
  
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
  
  # for each simulation of S, generate the conditional simulation, Sj
  for(i in 1:nsim) {
    #get predictions on untransformed (non-log) scale
    preds[,i] = exp(genSj(preds[,i]) + mu)
  }
  
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
         above7Prop=above7Prop, simCoords=simCoords))
}

#NOTE: Make sure to use joint model?
# mu97   sigmasq       phi     tausq 
# 1.5422466 0.1464766 0.1930667 0.0830306 
# > out$MLEs00$par
# mu00    sigmasq        phi      tausq 
# 0.72608418 0.15457486 0.27665466 0.01296299 
# > out$MLEsJoint$par
# mu97       mu00    sigmasq        phi      tausq 
# 1.55157286 0.72703505 0.13635714 0.30503468 0.05239249 
genPredictionsStandard = function(nsim=10000) {
  genPredictions(nsim=nsim, mu=1.55157286, sigma=sqrt(0.13635714), 
                 phi=0.30503468, tau=sqrt(0.05239249), nu=.5)
}

genPredictionsPreferential = function(nsim=10000) {
  genPredictions(nsim=nsim, mu=mu97MLE, sigma=sqrt(varianceMLE), 
                 phi=scaleMLE, tau=sqrt(nuggetMLE), nu=.5)
}

makePlots = function(predsPref, predsSD=NULL) {
  doSD = !is.null(predsSD)
  ##### first plot quantiles over space NOTE: domain over allCoords, or coords97?
  
  # get gridded coordinates that the quantiles are defined over (now unnecessary)
#   xrange = seq(win$xrange[1], win$xrange[2], l=predRes)
#   yrange = seq(win$yrange[1], win$yrange[2], l=predRes)
#   simCoords = make.surface.grid(list(x=xrange, y=yrange))
#   win.poly = matrix(unlist(win$bdry[[1]]), ncol=2)
#   simCoords = simCoords[in.poly(simCoords, win.poly),]
  
  # make quantile plots
  if(doSD) {
    quilt.plot(predsSD$simCoords, predsSD$quant5, zlim=c(0,11), xlab="", ylab="", main="")
    points(coords97)
    quilt.plot(predsSD$simCoords, predsSD$quant50, zlim=c(0,11), xlab="", ylab="", main="")
    points(coords97)
    quilt.plot(predsSD$simCoords, predsSD$quant95, zlim=c(0,31), xlab="", ylab="", main="")
    points(coords97)
  }
  
  quilt.plot(predsPref$simCoords, predsPref$quant5, zlim=c(0,11), xlab="", ylab="", main="")
  points(coords97)
  quilt.plot(predsPref$simCoords, predsPref$quant50, zlim=c(0,11), xlab="", ylab="", main="")
  points(coords97)
  quilt.plot(predsPref$simCoords, predsPref$quant95, zlim=c(0,31), xlab="", ylab="", main="")
  points(coords97)
  
  ##### make histograms of areal proportion above each level
  
  #level 3
  if(doSD) 
    sdAbove3 = hist(predsSD$above3Prop, plot=FALSE)
  prefAbove3 = hist(predsPref$above3Prop, plot=FALSE)
  plot(prefAbove3$mids, prefAbove3$density, type="l", xlab="Proportion", ylab="Density", 
       xlim=c(0,1), ylim=c(0, 6), 
       main=expression("Areal proportion over" ~ italic(3*mu*g) * "/(" * italic(g) ~ "dry weight)"))
  if(doSD) 
    lines(sdAbove3$mids, sdAbove3$density, lty=2)
  
  #level 5
  if(doSD) 
    sdAbove5 = hist(predsSD$above5Prop, plot=FALSE)
  prefAbove5 = hist(predsPref$above5Prop, plot=FALSE)
  plot(prefAbove5$mids, prefAbove5$density, type="l", xlab="Proportion", ylab="Density", 
       xlim=c(0,.9), ylim=c(0, 10), 
       main=expression("Areal proportion over" ~ italic(5*mu*g) * "/(" * italic(g) ~ "dry weight)"))
  if(doSD) 
    lines(sdAbove5$mids, sdAbove5$density, lty=2)
  
  #level 7
  if(doSD) 
    sdAbove7 = hist(predsSD$above7Prop, plot=FALSE)
  prefAbove7 = hist(predsPref$above7Prop, plot=FALSE)
  plot(prefAbove7$mids, prefAbove7$density, type="l", xlab="Proportion", ylab="Density", 
       xlim=c(0,.6), ylim=c(0, 30), 
       main=expression("Areal proportion over" ~ italic(7*mu*g) * "/(" * italic(g) ~ "dry weight)"))
  if(doSD) 
    lines(sdAbove7$mids, sdAbove7$density, lty=2)
#   hist(test$above3Prop, xlim=c(0, 1), freq = FALSE, xlab="Proportion", 
#   main=expression("Areal proportion over" ~ italic(3*mu*g) * "/(" * italic(g) ~ "dry weight)"))
#   hist(test$above5Prop, xlim=c(0, .7), freq = FALSE, xlab="Proportion", 
#   main=expression("Areal proportion over" ~ italic(5*mu*g) * "/(" * italic(g) ~ "dry weight)"))
#   hist(test$above7Prop, xlim=c(0, .3), freq = FALSE, xlab="Proportion", 
#   main=expression("Areal proportion over" ~ italic(7*mu*g) * "/(" * italic(g) ~ "dry weight)"))
}



# NOTE: for these, need histograms for pref and non-pref on same graphs
# hist(test$above3Prop, xlim=c(0, 1), freq = FALSE, xlab="Proportion", 
# main=expression("Areal proportion over" ~ italic(3*mu*g) * "/(" * italic(g) ~ "dry weight)"))
# hist(test$above5Prop, xlim=c(0, .7), freq = FALSE, xlab="Proportion", 
# main=expression("Areal proportion over" ~ italic(5*mu*g) * "/(" * italic(g) ~ "dry weight)"))
# hist(test$above7Prop, xlim=c(0, .3), freq = FALSE, xlab="Proportion", 
# main=expression("Areal proportion over" ~ italic(7*mu*g) * "/(" * italic(g) ~ "dry weight)"))


