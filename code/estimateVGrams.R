library(fields)
resGP <<- 64

# function to estimate empirical variogram of one PP
getVario = function(PP, GP, breaks=seq(0, 1, length=100), 
                    xrange=seq(0, 1, length=resGP), 
                    yrange=seq(0, 1, length=resGP)) {
  coords = cbind(PP$x, PP$y)
  
  #####get the associated marks for the PP at the GP
  
  # round PP coordinates in coords to be on GP coordinate grid
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords[,1], GP$xcol)
  roundY = roundToRange(coords[,2], GP$yrow)
  roundCoords = cbind(roundX, roundY)
  
  # find index of GP coords that are the same as the PP coords
  #GPCoords = attr(GP, "coords")
  GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords, GPCoords)
  
  #get PP marks from GP using the indices found above
  if(is.null(GP$noTranspose))
    PPmarks = t(matrix(GP$v, ncol=length(GP$xcol)))[inds]
  else
    PPmarks = matrix(GP$v, ncol=length(GP$xcol))[inds]
  
  # now compute the variogram (don't use rounded coordinates, though)
  return(vgram(roundCoords, PPmarks, breaks=breaks))
}

#function to estimate bias and standard deviations of variogram estimates 
#and make corresponding plots.  simulatePP.R must be sourced for this 
#function to work.  User must input observation locations (coords) and the 
#observation values themselves (ys) along with number of simulations (nsims)
varioAnalysis = function(nsims = 500) {
  mu=4
  
  # get bin breaks and centers for variograms
  breaks = seq(0, 1, length=40) #VG bin breaks
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  
  #generate GPs and PPs and calculate variograms for each PP/GP pair (also check number of observations in each bin)
  print("Calculating variograms")
  unifVGs = matrix(nrow=nsims, ncol=length(centers))
  clustVGs = matrix(nrow=nsims, ncol=length(centers))
  prefVGs = matrix(nrow=nsims, ncol=length(centers))
  unifDs = matrix(nrow=nsims, ncol=length(centers))
  clustDs = matrix(nrow=nsims, ncol=length(centers))
  prefDs = matrix(nrow=nsims, ncol=length(centers))
  unifNs = matrix(nrow=nsims, ncol=length(centers))
  clustNs = matrix(nrow=nsims, ncol=length(centers))
  prefNs = matrix(nrow=nsims, ncol=length(centers))
  
  # generate PPs and GPs
  unifPPs = genUnifPP(numSamples = nsims)
  clustPPs = genClusterPP2(numSamples = nsims)
  prefPPs = genPreferentialPP2(numSamples = nsims)
  GPs = prefPPs$GPs
  # GPsTest = genGPTest2(nsim=nsims)
  unifVGrams = list()
  
  for(i in 1:nsims) {
    print(paste0("generating ", i, "th set of simulations"))
    
    #get PP i of each type
    unifPP = unifPPs[[i]]
    clustPP = clustPPs$sims[[i]]
    prefPP = prefPPs$sims[[i]]
    
    #get corresponding GP
    GP = GPs[[i]]
    
    #compute empirical variograms
    unifVGrams[[i]] = getVario(unifPP, GP, breaks=breaks)
    unifVG = unifVGrams[[i]]
    clustVG = getVario(clustPP, GP, breaks=breaks)
    prefVG = getVario(prefPP, GP, breaks=breaks)
    unifVGStats = getVGMean(unifVG, breaks=breaks)
    clustVGStats = getVGMean(clustVG, breaks=breaks)
    prefVGStats = getVGMean(prefVG, breaks=breaks)
    unifVGs[i,] = unifVGStats$ys
    clustVGs[i,] = clustVGStats$ys
    prefVGs[i,] = prefVGStats$ys
    unifDs[i,] = unifVGStats$meanCenters
    clustDs[i,] = clustVGStats$meanCenters
    prefDs[i,] = prefVGStats$meanCenters
    unifNs[i,] = getVGMean(unifVG, breaks=breaks, statFun="length", statArgs=NULL)$ys
    clustNs[i,] = getVGMean(clustVG, breaks=breaks, statFun="length", statArgs=NULL)$ys
    prefNs[i,] = getVGMean(prefVG, breaks=breaks, statFun="length", statArgs=NULL)$ys
  }
  
  # get weights for averages for any given bin
  unifTotNs = colSums(unifNs)
  clustTotNs = colSums(clustNs)
  prefTotNs = colSums(prefNs)
  unifWeights = sweep(unifNs, 2, unifTotNs, "/")
  clustWeights = sweep(clustNs, 2, clustTotNs, "/")
  prefWeights = sweep(prefNs, 2, prefTotNs, "/")
  
  #compute average centers
  unifMeanCenters = colSums(unifDs * unifWeights, na.rm=TRUE)
  clustMeanCenters = colSums(clustDs * clustWeights, na.rm=TRUE)
  prefMeanCenters = colSums(prefDs * prefWeights, na.rm=TRUE)
  
  #get theoretical variogram function based on true Matern parameters
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  # exp(mu + sigma^2/2)
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 2
  # varEst = beta^2*sigmasq
  theorVG = function(x) {
    sigmasq*(1 - Matern(x, range=phi, nu=kappa)) #since tausq = 0
  }
  trueYs = theorVG(centers)
  
  # calculate bias of empirical VGs
  print("Computing variogram bias and standard deviations")
  unifBias = apply(sweep(unifVGs, 2, theorVG(unifMeanCenters)), 2, mean, na.rm=TRUE)
  clustBias = apply(sweep(clustVGs, 2, theorVG(clustMeanCenters)), 2, mean, na.rm=TRUE)
  prefBias = apply(sweep(prefVGs, 2, theorVG(prefMeanCenters)), 2, mean, na.rm=TRUE)
  
  # also calculate bias using weighted average
#   unifWeightedBias = colSums(sweep(unifVGs, 2, trueYs)*unifWeights)
#   clustWeightedBias = colSums(sweep(clustVGs, 2, trueYs)*clustWeights)
#   prefWeightedBias = colSums(sweep(prefVGs, 2, trueYs)*prefWeights)
  unifWeightedBias = colSums(sweep(unifVGs, 2, theorVG(unifMeanCenters))*unifWeights, na.rm=TRUE)
  clustWeightedBias = colSums(sweep(clustVGs, 2, theorVG(clustMeanCenters))*clustWeights, na.rm=TRUE)
  prefWeightedBias = colSums(sweep(prefVGs, 2, theorVG(prefMeanCenters))*prefWeights, na.rm=TRUE)
  
  # calculate standard deviation of empirical VGs
  unifSDs = apply(unifVGs, 2, sd, na.rm=TRUE)
  clustSDs = apply(clustVGs, 2, sd, na.rm=TRUE)
  prefSDs = apply(prefVGs, 2, sd, na.rm=TRUE)
  
  ##### make plot emulating Fig. 2 from paper
  #plot of bias +- 1.96 SEs
  print("Generating plots")
  par(mfrow=c(1,2))
  plot(unifMeanCenters, unifBias + 1.96*unifSDs/sqrt(nsims), main="", ylab="Bias", type="l", xlab="u", ylim=c(-1, .5))
  lines(unifMeanCenters, unifBias - 1.96*unifSDs/sqrt(nsims))
  lines(clustMeanCenters, clustBias + 1.96*clustSDs/sqrt(nsims), lty=2)
  lines(clustMeanCenters, clustBias - 1.96*clustSDs/sqrt(nsims), lty=2)
  lines(prefMeanCenters, prefBias + 1.96*prefSDs/sqrt(nsims), lty=3)
  lines(prefMeanCenters, prefBias - 1.96*prefSDs/sqrt(nsims), lty=3)
  abline(h=0)
  
  #plot of SDs
  plot(centers, unifSDs, main="", ylab="SD", xlab="u", ylim=c(0, 2), type="l")
  lines(centers, clustSDs, lty=2)
  lines(centers, prefSDs, lty=3)
  
  # do same plots but with weighted bias
#   plot(centers, unifWeightedBias + 1.96*unifSDs/sqrt(unifTotNs), main="", ylab="Bias", type="l", xlab="u", ylim=c(-1, .5))
#   lines(centers, unifWeightedBias - 1.96*unifSDs/sqrt(unifTotNs))
#   lines(centers, clustWeightedBias + 1.96*clustSDs/sqrt(clustTotNs), lty=2)
#   lines(centers, clustWeightedBias - 1.96*clustSDs/sqrt(clustTotNs), lty=2)
#   lines(centers, prefWeightedBias + 1.96*prefSDs/sqrt(prefTotNs), lty=3)
#   lines(centers, prefWeightedBias - 1.96*prefSDs/sqrt(prefTotNs), lty=3)
  
  #plot of SDs
#   plot(centers, unifSDs, main="", ylab="SD", xlab="u", ylim=c(0, 2), type="l")
#   lines(centers, clustSDs, lty=2)
#   lines(centers, prefSDs, lty=3)
#   par(mfrow=c(1,1))
  
  ##### now recreate Fig. 1 from paper
  GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
  par(mfrow=c(1,1), family="serif")
  png("unifScheme.png", width=500, height=500)
  quilt.plot(GPCoords, GP$v, main="Uniform Sample", nx=resGP, ny=resGP)
  axis(1, at=seq(0, 1, l=3))
  axis(2, at=seq(0, 1, l=3))
  plot(unifPP, add=TRUE, pch=19, cex=.5)
  dev.off()
  
#   quilt.plot(GPCoords, GP$v, main="Uniform Sample", nx=resGP, ny=resGP, cex.main=3, cex.axis=2, legend.cex=2)
#   axis(1, at=seq(0, 1, l=3))
#   axis(2, at=seq(0, 1, l=3))
#   plot(unifPP, add=TRUE, pch=19, cex=1)
  
  png("clustScheme.png", width=500, height=500)
  quilt.plot(GPCoords, GP$v, main="Clustered Sample", nx=resGP, ny=resGP)
  points(clustPP, pch=19, cex=.5)
  dev.off()
  
  png("prefScheme.png", width=500, height=500)
  quilt.plot(GPCoords, GP$v, main="Preferential Sample", nx=resGP, ny=resGP)
  # plot(prefPP, pch=19, cex=.5)
  points(prefPP$x, prefPP$y, pch=19, cex=.5)
  dev.off()
  
  par(mfrow=c(1,1))
  
  return(list(centers=centers, breaks=breaks, unifNs=unifNs, clustNs=clustNs, prefNs=prefNs,
              unifSDs=unifSDs, clustSDs = clustSDs, prefSDs=prefSDs, unifBias=unifBias, 
              clustBias=clustBias, prefBias=prefBias, unifVGs=unifVGs, clustVGs=clustVGs, 
              prefVGs=prefVGs, trueYs=trueYs, nsims=nsims))
}

#recreate figure 5
makeVGPlot = function(coords, logDat, tau=sqrt(nuggetMLE), sigma=sqrt(sillMLE), kappa=0.5, phi=scaleMLE, 
                      main="Variogram", ylim=c(0, .3), breaks=seq(0, 1.2, l=25)) {
  
  #theoretical variogram given x (distance)
  trueVG = function(x) {
    sigma^2*(1 - Matern(x, nu=kappa, range=phi)) + tau^2
  }
  
  #get empirical variogram
  empVG = vgram(coords, logDat)
  
  #make plots
  VGMean = getVGMean(empVG, breaks=breaks)
  centers = VGMean$centers
  ys = VGMean$ys
  plot(centers, ys, main=main, ylim=ylim, pch=19, cex=.75, ylab="Variance", xlab="Distance (100's of km)")
  xs = seq(min(breaks), max(breaks), l=100)
  lines(xs, trueVG(xs))
}
# makeVGPlot(coords97, log97, tau=sqrt(out$MLEs97$par[4]), sigma=sqrt(out$MLEs97$par[2]), phi=out$MLEs97$par[3], main="1997 Variogram")
# makeVGPlot(coords00, log00, tau=sqrt(out$MLEs00$par[4]), sigma=sqrt(out$MLEs00$par[2]), phi=out$MLEs00$par[3], main="2000 Variogram")

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
#extensions to fields package making sure variogram empirical means are 
#plotted rather than medians

meanVG = function(VG, minD=-Inf, maxD=Inf, statFun="mean", ...) {
  ind = (VG$d > minD) & (VG$d < maxD)
  do.call(statFun, c(list(VG$vgram[ind]), list(...)))
}

getVGMean = function(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1), 
                     statFun="mean", statArgs=list(na.rm=TRUE)) 
{
  meansFromBreak = function(breakBounds = c(-Inf, Inf)) {
    do.call("meanVG", c(list(x, minD=breakBounds[1], maxD=breakBounds[2], statFun = statFun), statArgs))
  }
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  ys = apply(breakBounds, 1, meansFromBreak)
  x$vgram = x$d
  meanCenters = apply(breakBounds, 1, meansFromBreak)
  
  return(list(centers=centers, meanCenters=meanCenters, ys=ys, type=x$type))
}

plotVGStat = function(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1), 
                      add = FALSE, statFun="mean", statArgs=list(na.rm=TRUE), ...) 
{
  otherArgs = list(...)
  type = x$type
  if (is.null(otherArgs$ylab)) {
    if (type == "variogram")
      ylab = "sqrt(Variance)"
    else if (type == "covariogram" || type == "cross-covariogram") 
      ylab = "Covariance"
    else if (type == "correlogram" || type == "cross-correlogram") 
      ylab = "Correlation"
    else stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  else {
    ylab = otherArgs$ylab
    otherArgs$ylab = NULL
  }
  if (is.null(otherArgs$xlab)) 
    xlab = "Distance (Miles)"
  else {
    xlab = otherArgs$xlab
    otherArgs$xlab = NULL
  }
  if (is.null(otherArgs$main)) {
    if (type == "variogram") 
      main = "Empirical Variogram"
    else if (type == "covariogram") 
      main = "Empirical Covariogram"
    else if (type == "correlogram") 
      main = "Empirical Correlogram"
    else if (type == "cross-covariogram") 
      main = "Empirical Cross-Covariogram"
    else if (type == "cross-correlogram") 
      main = "Empirical Cross-Correlogram"
    else stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  else {
    main = otherArgs$main
    otherArgs$main = NULL
  }
  if (is.null(otherArgs$ylim)) {
    if (type == "correlogram" || type == "cross-correlogram") 
      ylim = c(-1, 1)
    else ylim = NULL
  }
  else {
    ylim = otherArgs$ylim
    otherArgs$ylim = NULL
  }
  if (is.null(otherArgs$type)) 
    type = "o"
  else {
    type = otherArgs$type
    otherArgs$type = NULL
  }
  
  # get variogram statistic we're interested in (usually mean of point cloud in bin)
  out = getVGMean(x, breaks=breaks, statFun=statFun, statArgs=statArgs)
  ys = out$ys
  centers = out$centers
#   if(x$type == "variogram")
#     ys=sqrt(ys)
  notNas = !is.na(ys)
  centers = centers[notNas]
  ys = ys[notNas]
  if (!add) 
    do.call(plot, c(list(centers, ys, main = main, xlab = xlab, 
                         ylab = ylab, type = type, ylim = ylim), otherArgs))
  else do.call(lines, c(list(centers, ys, main = main, xlab = xlab, 
                             ylab = ylab, type = type, ylim = ylim), otherArgs))
  
  return(invisible(list(centers=centers, ys=ys, type=type)))
}

genGPPPTest = function(mu=4, sigma=sqrt(1.5), phi=.15, kappa=1, nsim=1, npts=100) {
  xrange=seq(0, 1, length=resGP)
  yrange=seq(0, 1, length=resGP)
  
  # get bin breaks and centers for variograms
  breaks = seq(0, 1, length=40) #VG bin breaks
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  
  unifPP = genPPTest(npts=npts)
  GPs = genGPTest(unifPP$coords, mu=mu, sigma=sigma, phi=phi, kappa=kappa, nsim=nsim)
  
  unifVGs = matrix(nrow=nsim, ncol=length(centers))
  unifDs = matrix(nrow=nsim, ncol=length(centers))
  unifNs = matrix(nrow=nsim, ncol=length(centers))
  for(i in 1:nsim) {
    print(paste0("generating ", i, "th set of simulations"))
    GP = list(xcol=xrange, yrow=yrange, v=GPs[,i])
    PP = unifPP
    
    #compute empirical variograms
    unifVG = getVarioTest(unifPP, GP, breaks=breaks)
    unifVGStats = getVGMean(unifVG, breaks=breaks)
    unifVGs[i,] = unifVGStats$ys
    unifDs[i,] = unifVGStats$meanCenters
    unifNs[i,] = getVGMean(unifVG, breaks=breaks, statFun="length", statArgs=NULL)$ys
  }
  
  # get weights for averages for any given bin
  unifTotNs = colSums(unifNs)
  unifWeights = sweep(unifNs, 2, unifTotNs, "/")
  
  #compute average centers
  unifMeanCenters = colSums(unifDs * unifWeights, na.rm=TRUE)
  
  #get theoretical variogram function based on true Matern parameters
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  # exp(mu + sigma^2/2)
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 2
  # varEst = beta^2*sigmasq
  theorVG = function(x) {
    sigmasq*(1 - Matern(x, range=phi, nu=kappa)) #since tausq = 0
  }
  trueYs = theorVG(centers)
  
  # calculate bias of empirical VGs
  print("Computing variogram bias and standard deviations")
  unifBias = apply(sweep(unifVGs, 2, theorVG(unifMeanCenters)), 2, mean, na.rm=TRUE)
  
  # calculate standard deviation of empirical VGs
  unifSDs = apply(unifVGs, 2, sd, na.rm=TRUE)
  
  ##### make plot emulating Fig. 2 from paper
  #plot of bias +- 1.96 SEs
  print("Generating plots")
  par(mfrow=c(1,2))
  plot(centers, unifBias + 1.96*unifSDs/sqrt(nsim), main="", ylab="Bias", type="l", xlab="u", ylim=c(-1, .5))
  lines(centers, unifBias - 1.96*unifSDs/sqrt(nsim))
  
  #plot of SDs
  plot(centers, unifSDs, main="", ylab="SD", xlab="u", ylim=c(0, 2), type="l")
  par(mfrow=c(1,1))
}

genGPTest2 = function(nsim=1, nx=resGP, ny=resGP, mu=4, sigma=sqrt(1.5), phi=.15, kappa=1, beta=2, tausq=0) {
  xs = seq(0, 1, length=nx)
  ys = seq(0, 1, length=ny)
  coords = make.surface.grid(list(x=xs, y=ys))
  
  Sigma = stationary.cov(coords, Covariance="Matern", theta=phi, nu=kappa)*sigma^2
  L = t(chol(Sigma))
  Z = matrix(rnorm(nrow(L)*nsim), nrow=nrow(L), ncol=nsim)
  simMat = L %*% Z + mu
  
  sims = list()
  for(i in 1:nsim) {
    sim = simMat[,i]
    sims[[i]] = list(xcol=xs, yrow=ys, v=sim)
  }
  
  return(sims)
}

genGPTest = function(coords, mu=4, sigma=sqrt(1.5), phi=.15, kappa=1, nsim=1) {
  Sigma = stationary.cov(coords, Covariance="Matern", theta=phi, nu=kappa)*sigma^2
  L = t(chol(Sigma))
  Z = matrix(rnorm(nrow(L)*nsim), nrow=nrow(L), ncol=nsim)
  return(L %*% Z + mu)
}

genPPTest = function(npts=100) {
  coords = matrix(runif(2*npts), ncol=2)
  list(x=coords[,1], y=coords[,2], coords=coords)
}

# function to estimate empirical variogram of one PP
getVarioTest = function(PP, GP, breaks=seq(0, 1, length=100), 
                    xrange=seq(0, 1, length=resGP), 
                    yrange=seq(0, 1, length=resGP)) {
  coords = cbind(PP$x, PP$y)
  
  #####get the associated marks for the PP at the GP
  
  # round PP coordinates in coords to be on GP coordinate grid
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords[,1], GP$xcol)
  roundY = roundToRange(coords[,2], GP$yrow)
  roundCoords = cbind(roundX, roundY)
  
  # find index of GP coords that are the same as the PP coords
  #GPCoords = attr(GP, "coords")
  GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords, GPCoords)
  
  # now compute the variogram (don't use rounded coordinates, though)
  return(vgram(roundCoords, GP$v, breaks=breaks))
}

#function for testing GP variance
getVar = function(GP) {
var(c(GP$v)) }
getTotVar = function(GPs) {
  vals = c()
  for(i in 1:length(GPs)) {
    vals = c(vals, GPs[[i]]$v)
  }
  var(vals)
}

plotSquareVG = function(VG, ...) {
  VG$vgram=VG$vgram^2
  do.call("plot.vgram", c(list(VG), list(...)))
}

concatVGs = function(VGs) {
  VG = VGs[[1]]
  for(i in 2:length(VGs)) {
    tmp = VGs[[i]]
    VG$d = c(VG$d, tmp$d)
    VG$vgram = c(VG$vgram, tmp$vgram)
  }
  return(VG)
}

fullGeoRTest = function(PPs, GPs) {
  n = length(PPs)
  nbins=39
  
  ##### concatenate VGrams
  VG = geoRVGramTest(1, PPs, GPs)
  for(i in 2:n) {
    thisVG = geoRVGramTest(i, PPs, GPs)
    for(b in 1:nbins) {
      VG$bin.cloud[[b]] = c(VG$bin.cloud[[b]], thisVG$bin.cloud[[b]])
    }
  }
  VG$v = sapply(VG$bin.cloud, mean, na.rm=TRUE)
  
  ##### plot results and return
  plot(VG)
  
  return(VG)
}

geoRVGramTest = function(i, PPs, GPs) {
  PP = PPs[[i]]
  GP = GPs[[i]]
  coords = cbind(PP$x, PP$y)
  
  #####get the associated marks for the PP at the GP
  
  # round PP coordinates in coords to be on GP coordinate grid
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords[,1], GP$xcol)
  roundY = roundToRange(coords[,2], GP$yrow)
  roundCoords = cbind(roundX, roundY)
  
  # find index of GP coords that are the same as the PP coords
  #GPCoords = attr(GP, "coords")
  GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords, GPCoords)
  
  #get PP marks from GP using the indices found above
  PPmarks = t(matrix(GP$v, ncol=length(GP$xcol)))[inds]
  
  ##### get and plot variogram
  breaks = seq(0, 1, length=40) #VG bin breaks
  gDat = as.geodata(cbind(coords, PPmarks))
  
  VG = variog(gDat, bin.cloud=TRUE, breaks=breaks)
#   plot(VG, ylim=c(0, 2))
#   abline(h=sqrt(1.5))
  return(VG)
}

getPPMarks = function(PP, GP) {
  coords = cbind(PP$x, PP$y)
  
  #####get the associated marks for the PP at the GP
  
  # round PP coordinates in coords to be on GP coordinate grid
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords[,1], GP$xcol)
  roundY = roundToRange(coords[,2], GP$yrow)
  roundCoords = cbind(roundX, roundY)
  
  # find index of GP coords that are the same as the PP coords
  #GPCoords = attr(GP, "coords")
  GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  inds = findIndex(roundCoords, GPCoords)
  
  #get PP marks from GP using the indices found above
  PPmarks = t(matrix(GP$v, ncol=length(GP$xcol)))[inds]
  return(PPmarks)
}

getAllGPVals = function(GPs) {
  marks = c()
  for(i in 1:length(GPs)) {
    marks = c(marks, GPs[[i]]$v)
  }
  return(marks)
}

getGPMeans = function(GPs) {
  marks = c()
  for(i in 1:length(GPs)) {
    marks = c(marks, mean(GPs[[i]]$v))
  }
  return(marks)
}

getAllPPMarks = function(PPs, GPs) {
  marks = c()
  for(i in 1:length(PPs)) {
    marks = c(marks, getPPMarks(PPs[[i]], GPs[[i]]))
  }
  return(marks)
}
