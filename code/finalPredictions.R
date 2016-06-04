# Functions for generating the final predictions of the paper: Fig. 7-9
# NOTE: the figures are only for 1997
predRes <<- 64

genPredictions = function(nsim=10000, nx=predRes, ny=predRes, mu=mu97MLE, mu2=mu00MLE, sigma=sqrt(varianceMLE), 
                          phi=scaleMLE, tau=sqrt(nuggetMLE), nu=.5) {
  logAll = c(log97, log00)
  
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
  # predsAll = preds
  
  # n and N as defined in paper: n is number of observations, N is number of pts in grid
  n = length(log97)
  # nAll = n + length(log00)
  N = nrow(simCoords)
  # muAll = c(rep(mu, n), rep(mu2, nAll - n))
  
  # round true coordinates to the lattice coordinates
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords97[,1], xrange)
  roundY = roundToRange(coords97[,2], yrange)
  # roundXAll = roundToRange(allCoords[,1], xrange)
  # roundYAll = roundToRange(allCoords[,2], yrange)
  roundCoords97 = cbind(roundX, roundY)
  # roundCoordsAll =cbind(roundXAll, roundYAll)
  
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords97, simCoords)
  # indsAll = findIndex(roundCoordsAll, simCoords)
  
  #generate marginal covariance matrix of S, Sigma, for the given parameters
  distMat = rdist(simCoords)
  Sigma = Exp.cov(simCoords, theta=phi, distMat=distMat)*sigma^2
  
  #Now compute Sigma0, the covariance matrix of Y
  #Sigma0 = as.matrix(C %*% Sigma %*% Ct)
  Sigma0 = Sigma[inds, inds]
  Sigma10 = Sigma[,inds]
  # Sigma0All = Sigma[indsAll, indsAll]
  # Sigma10All = Sigma[,indsAll]
  
  # add nugget error to Sigma with fast internal function of fields in C
  invisible(.Call("addToDiagC", Sigma0, as.double(rep(tau^2, nrow(Sigma0))), nrow(Sigma0)))
  # invisible(.Call("addToDiagC", Sigma0All, as.double(rep(tau^2, nrow(Sigma0All))), nrow(Sigma0All)))
  
  #Cholesky decomposition of Sigma0 required for likelihood
  Sigma0U = chol(Sigma0)
  # Sigma0AllU = chol(Sigma0All)
  
  # compute muc, expectation of Sj based on this formula:
  # Sigma %*% Ct %*% Sigma0^-1 %*% (dat - mu)
  vec = log97 - mu
  z = backsolve(Sigma0U, vec, transpose=TRUE)
  x = cbind(backsolve(Sigma0U, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
  #muc = as.numeric(Sigma %*% Ct %*% x)
  muc = Sigma10 %*% x
  
#   vecAll = logAll - muAll
#   z = backsolve(Sigma0AllU, vecAll, transpose=TRUE)
#   x = cbind(backsolve(Sigma0AllU, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
#   #muc = as.numeric(Sigma %*% Ct %*% x)
#   mucAll = Sigma10All %*% x
  
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
  
  # same as genSj, but this also uses 2000 data and mu00
#   genSjJoint = function(S) {
#     Zs = rnorm(nAll)*tau
#     #vec = dat - mu97 + Zs - C %*% S
#     vec = logAll - muAll + Zs - S[indsAll]
#     z = backsolve(Sigma0AllU, vec, transpose=TRUE)
#     x = cbind(backsolve(Sigma0AllU, z)) #make sure x is a column vector TODO: MAKE SURE THIS WORKS *********
#     #Sj = as.numeric(S + Sigma %*% Ct %*% x)
#     Sj = S + Sigma10All %*% x
#     return(Sj)
#   }
  
  # for each simulation of S, generate the conditional simulation, Sj
  for(i in 1:nsim) {
    #get predictions on untransformed (non-log) scale
    preds[,i] = exp(genSj(preds[,i]) + mu)
    # predsAll[,i] = exp(genSjJoint(predsAll[,i]) + mu) #don't use muAll here because predictions are for 1997
  }
  
  # get 5%, 50%, and 95% quantiles
  quant5 = apply(preds, 1, quantile, probs=.05)
  quant50 = apply(preds, 1, quantile, probs=.5)
  quant95 = apply(preds, 1, quantile, probs=.95)
  
#   quant5All = apply(predsAll, 1, quantile, probs=.05)
#   quant50All = apply(predsAll, 1, quantile, probs=.5)
#   quant95All = apply(predsAll, 1, quantile, probs=.95)
  
  #get distributions for proportion larger than 3, 5, and 7 micrograms of lead
  propBigger = function(dat, thresh = 3) {
    mean(dat > thresh)
  }
  above3Prop = apply(preds, 2, propBigger, thresh=3)
  above5Prop = apply(preds, 2, propBigger, thresh=5)
  above7Prop = apply(preds, 2, propBigger, thresh=7)
  
#   above3PropAll = apply(predsAll, 2, propBigger, thresh=3)
#   above5PropAll = apply(predsAll, 2, propBigger, thresh=5)
#   above7PropAll = apply(predsAll, 2, propBigger, thresh=7)
  
#   return(list(simCoords=simCoords, 
#               quant5=quant5, quant50=quant50, quant95=quant95, above3Prop=above3Prop, above5Prop=above5Prop, 
#               above7Prop=above7Prop, 
#               quant5All=quant5All, quant50All=quant50All, quant95All=quant95All, above3PropAll=above3PropAll, above5PropAll=above5PropAll, 
#               above7PropAll=above7PropAll))
  return(list(simCoords=simCoords, 
              quant5=quant5, quant50=quant50, quant95=quant95, above3Prop=above3Prop, above5Prop=above5Prop, 
              above7Prop=above7Prop))
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

makePlotsJoint = function(predsPref, predsSD=NULL, zlim=c(0, 11)) {
  predsPrefAll = list(quant5=predsPref$quant5All, quant50=predsPref$quant50All, quant95=predsPref$quant95All, above3Prop=predsPref$above3PropAll, above5Prop=predsPref$above5PropAll, 
                      above7Prop=predsPref$above7PropAll, simCoords=predsPref$simCoords)
  predsPref = list(quant5=predsPref$quant5, quant50=predsPref$quant50, quant95=predsPref$quant95, above3Prop=predsPref$above3Prop, above5Prop=predsPref$above5Prop, 
                   above7Prop=predsPref$above7Prop, simCoords=predsPref$simCoords)
  if(!is.null(predsSD)) {
    predsSDAll = list(quant5=predsSD$quant5All, quant50=predsSD$quant50All, quant95=predsSD$quant95All, above3Prop=predsSD$above3PropAll, above5Prop=predsSD$above5PropAll, 
                      above7Prop=predsSD$above7PropAll, simCoords=predsSD$simCoords)
    predsSD = list(quant5=predsSD$quant5, quant50=predsSD$quant50, quant95=predsSD$quant95, above3Prop=predsSD$above3Prop, above5Prop=predsSD$above5Prop, 
                   above7Prop=predsSD$above7Prop, simCoords=predsSD$simCoords)
  }
  else{
    predsSD = NULL
    predsSDAll = NULL
  }
  
  makePlots(predsPref, predsSD)
  makePlots(predsPrefAll, predsSDAll, zlim=zlim)
}

makePlots = function(predsPref, predsSD=NULL, zlim=c(0, 31), addLegend=TRUE) {
  doSD = !is.null(predsSD)
  ##### first plot quantiles over space NOTE: domain over allCoords, or coords97?
  
  # get gridded coordinates that the quantiles are defined over (now unnecessary)
#   xrange = seq(win$xrange[1], win$xrange[2], l=predRes)
#   yrange = seq(win$yrange[1], win$yrange[2], l=predRes)
#   simCoords = make.surface.grid(list(x=xrange, y=yrange))
#   win.poly = matrix(unlist(win$bdry[[1]]), ncol=2)
#   simCoords = simCoords[in.poly(simCoords, win.poly),]
  
  # make quantile plots
  plotRes=62 # adjust plot resolution for best plotting
  quilt.plot(predsPref$simCoords, predsPref$quant5, main="Test")
  if(doSD) {
    quilt.plot(predsSD$simCoords, predsSD$quant5, zlim=zlim, xlab="", ylab="", main="", add.legend = addLegend, nx=plotRes, ny=plotRes)
    points(coords97, cex=.7)
    quilt.plot(predsSD$simCoords, predsSD$quant50, zlim=zlim, xlab="", ylab="", 
               main="Non-Preferential Predictions", add.legend = addLegend, nx=plotRes, ny=plotRes)
    points(coords97, cex=.7)
    quilt.plot(predsSD$simCoords, predsSD$quant95, zlim=zlim, xlab="", ylab="", main="", nx=plotRes, ny=plotRes)
    points(coords97, cex=.7)
  }
  
  quilt.plot(predsPref$simCoords, predsPref$quant5, zlim=zlim, xlab="", ylab="", main="", add.legend = addLegend, nx=plotRes, ny=plotRes)
  points(coords97, cex=.7)
  quilt.plot(predsPref$simCoords, predsPref$quant50, zlim=zlim, xlab="", ylab="", 
             main="Preferential Predictions", , add.legend = addLegend, nx=plotRes, ny=plotRes)
  points(coords97, cex=.7)
  quilt.plot(predsPref$simCoords, predsPref$quant95, zlim=zlim, xlab="", ylab="", main="", nx=plotRes, ny=plotRes)
  points(coords97, cex=.7)
  
  ##### make histograms of areal proportion above each level
  
  #level 3
  par(mar=c(5.5, 5, 3.1, 2.1))
  if(doSD) 
    sdAbove3 = hist(predsSD$above3Prop, plot=FALSE)
  prefAbove3 = hist(predsPref$above3Prop, plot=FALSE)
  plot(prefAbove3$mids, prefAbove3$density, type="l", xlab="Proportion", ylab="Density", 
       xlim=c(0,1), ylim=c(0, 12))#, 
       #main=expression("Areal proportion over" ~ italic(3*mu*g) * "/(" * italic(g) ~ "dry weight)"))
  if(doSD) 
    lines(sdAbove3$mids, sdAbove3$density, lty=2)
  
  #level 5
  if(doSD) 
    sdAbove5 = hist(predsSD$above5Prop, plot=FALSE)
  prefAbove5 = hist(predsPref$above5Prop, plot=FALSE)
  plot(prefAbove5$mids, prefAbove5$density, type="l", xlab="Proportion", ylab="Density", 
       xlim=c(0,.9), ylim=c(0, 10))#, 
       #main=expression("Areal proportion over" ~ italic(5*mu*g) * "/(" * italic(g) ~ "dry weight)"))
  if(doSD) 
    lines(sdAbove5$mids, sdAbove5$density, lty=2)
  
  #level 7
  if(doSD) 
    sdAbove7 = hist(predsSD$above7Prop, plot=FALSE)
  prefAbove7 = hist(predsPref$above7Prop, plot=FALSE)
  plot(prefAbove7$mids, prefAbove7$density, type="l", xlab="Proportion", ylab="Density", 
       xlim=c(0,.6), ylim=c(0, 30))#, 
       #main=expression("Areal proportion over" ~ italic(7*mu*g) * "/(" * italic(g) ~ "dry weight)"))
  if(doSD) 
    lines(sdAbove7$mids, sdAbove7$density, lty=2)
#   hist(test$above3Prop, xlim=c(0, 1), freq = FALSE, xlab="Proportion", 
#   main=expression("Areal proportion over" ~ italic(3*mu*g) * "/(" * italic(g) ~ "dry weight)"))
#   hist(test$above5Prop, xlim=c(0, .7), freq = FALSE, xlab="Proportion", 
#   main=expression("Areal proportion over" ~ italic(5*mu*g) * "/(" * italic(g) ~ "dry weight)"))
#   hist(test$above7Prop, xlim=c(0, .3), freq = FALSE, xlab="Proportion", 
#   main=expression("Areal proportion over" ~ italic(7*mu*g) * "/(" * italic(g) ~ "dry weight)"))
}

# predsSD = genPredictionsStandard()
# predsPref = genPredictionsPreferential()
# makePlots(predsPref, predsSD, zlim=c(0, 11.9), addLegend=FALSE)


# NOTE: for these, need histograms for pref and non-pref on same graphs
# hist(test$above3Prop, xlim=c(0, 1), freq = FALSE, xlab="Proportion", 
# main=expression("Areal proportion over" ~ italic(3*mu*g) * "/(" * italic(g) ~ "dry weight)"))
# hist(test$above5Prop, xlim=c(0, .7), freq = FALSE, xlab="Proportion", 
# main=expression("Areal proportion over" ~ italic(5*mu*g) * "/(" * italic(g) ~ "dry weight)"))
# hist(test$above7Prop, xlim=c(0, .3), freq = FALSE, xlab="Proportion", 
# main=expression("Areal proportion over" ~ italic(7*mu*g) * "/(" * italic(g) ~ "dry weight)"))


