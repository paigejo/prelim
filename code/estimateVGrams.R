library(fields)

# function to estimate empirical variogram of one PP
getVario = function(PP, GP, breaks=seq(0, 1, length=100), 
                    xrange=seq(0, 1, length=128), 
                    yrange=seq(0, 1, length=128)) {
  coords = cbind(PP$x, PP$y)
  
  #####get the associated marks for the PP at the GP
  
  # round PP coordinates in coords to be on GP coordinate grid
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundX = roundToRange(coords[,1], xrange)
  roundY = roundToRange(coords[,2], yrange)
  roundCoords = cbind(roundX, roundY)
  
  # find index of GP coords that are the same as the PP coords
  #GPCoords = attr(GP, "coords")
  GPCoords = make.surface.grid(list(x=GP$xcol, y=GP$yrow))
  findIndex = function(coords, len=128) {
    #convert coords to x and y grid indices
    coords = coords*(len-1) + 1
    xInd = coords[1]
    yInd = coords[2]
    return((yInd-1)*len + xInd)
  }
  inds = apply(roundCoords, 1, findIndex, len=length(xrange))
  
  #get PP marks from GP using the indices found above
  PPmarks = t(matrix(GP$v, ncol=length(GP$xcol)))[inds]
  
  # now compute the variogram (don't use rounded coordinates, though)
  return(vgram(coords, PPmarks, breaks=breaks))
}

#function to estimate bias and standard deviations of variogram estimates 
#and make corresponding plots.  simulatePP.R must be sourced for this 
#function to work.  User must input observation locations (coords) and the 
#observation values themselves (ys) along with number of simulations (nsims)
varioAnalysis = function(nsims = 500) {
  mu=4
  
  ##### Generate nsims (500 in paper) PPs of each type with corresponding GPs
  print("Generating GPs and PPs")
  unifPPs = genUnifPP(nsims)
  clustPPs = genClusterPP(nsims)
  prefPPs = genPreferentialPP(nsims)
  GPs = prefPPs$GPs
  
  # get bin breaks and centers for variograms
  breaks = seq(0, 1, length=20) #VG bin breaks
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  
  #calculate variograms for each PP
  print("Calculating variograms")
  unifVGs = matrix(nrow=nsims, ncol=length(centers))
  clustVGs = matrix(nrow=nsims, ncol=length(centers))
  prefVGs = matrix(nrow=nsims, ncol=length(centers))
  for(i in 1:nsims) {
    #get PP i of each type
    unifPP = unifPPs[[i]]
    clustPP = clustPPs$sims[[i]]
    prefPP = prefPPs$sims[[i]]
    
    #get corresponding GP
    GP = GPs[[i]]
    
    #compute empirical variograms
    unifVGs[i,] = getVGMean(getVario(unifPP, GP, breaks=breaks), breaks=breaks)$ys
    clustVGs[i,] = getVGMean(getVario(clustPP, GP, breaks=breaks), breaks=breaks)$ys
    prefVGs[i,] = getVGMean(getVario(prefPP, GP, breaks=breaks), breaks=breaks)$ys
  }
  
  #get theoretical variogram function based on true Matern parameters
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  # exp(mu + sigma^2/2)
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 2
  varEst = beta^2*sigmasq
  theorVG = function(x) {
    sigmasq*(1 - Matern(x, range=phi, nu=kappa)) #since tausq = 0
  }
  trueYs = theorVG(centers)
  
  # calculate bias of empirical VGs
  print("Computing variogram bias and standard deviations")
  unifBias = apply(sweep(unifVGs, 2, trueYs), 2, mean, na.rm=TRUE)
  clustBias = apply(sweep(clustVGs, 2, trueYs), 2, mean, na.rm=TRUE)
  prefBias = apply(sweep(prefVGs, 2, trueYs), 2, mean, na.rm=TRUE)
  
  # calculate standard deviation of empirical VGs
  unifSDs = apply(unifVGs, 2, sd)
  clustSDs = apply(clustVGs, 2, sd)
  prefSDs = apply(prefVGs, 2, sd)
  
  ##### make plot emulating Fig. 2 from paper
  #plot of bias +- 1.96 SEs
  print("Generating plots")
  par(mfrow=c(1,2))
  plot(centers, unifBias + 1.96*unifSDs/sqrt(nsims), main="", ylab="Bias", type="l", xlab="u", ylim=c(-1, .5))
  lines(centers, unifBias - 1.96*unifSDs/sqrt(nsims))
  lines(centers, clustBias + 1.96*clustSDs/sqrt(nsims), lty=2)
  lines(centers, clustBias - 1.96*clustSDs/sqrt(nsims), lty=2)
  lines(centers, prefBias + 1.96*prefSDs/sqrt(nsims), lty=3)
  lines(centers, prefBias - 1.96*prefSDs/sqrt(nsims), lty=3)
  
  #plot of SDs
  plot(centers, unifSDs, main="", ylab="SD", xlab="u", ylim=c(0, 2), type="l")
  lines(centers, clustSDs, lty=2)
  lines(centers, prefSDs, lty=3)
  
  invisible(NULL)
}

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

getVGMean = function(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1)) 
{
  meansFromBreak = function(breakBounds = c(-Inf, Inf)) {
    meanVG(x, minD=breakBounds[1], maxD=breakBounds[2], na.rm=TRUE)
  }
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  ys = apply(breakBounds, 1, meansFromBreak)
  if(x$type == "variogram")
    ys=sqrt(ys)
  
  return(list(centers=centers, ys=ys, type=x$type))
}

plotVGMean = function(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1), 
                      add = FALSE, ...) 
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
  meansFromBreak = function(breakBounds = c(-Inf, Inf)) {
    meanVG(x, minD=breakBounds[1], maxD=breakBounds[2], na.rm=TRUE)
  }
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  ys = apply(breakBounds, 1, meansFromBreak)
  if(x$type == "variogram")
    ys=sqrt(ys)
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