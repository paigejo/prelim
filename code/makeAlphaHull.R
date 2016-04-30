library(alphahull)

#make concave hull around coordinates
makePredHull = function(coords, alpha=1) {
  ashape = ahull(coords, alpha=alpha)
  
  indx=ashape$arcs[,"end1"]  
  hullPts <- coords[indx,]
  
  hullPoly = rbind(hullPts, hullPts[1,])
  
  return(list(ashape=ashape, hullPts=hullPts, hullPoly=hullPoly))
}

#plot the concave hull and the data coordinates
plotAlphaHull = function(hullPts, coords=NULL) {
  
  #plot prediction points
  plotSubPoly(rbind(hullPts, hullPts[1,]), coords)
  
  invisible(NULL)
}

plotSubPoly = function(bds, coords, range=NULL) {
  if(! is.null(range))
    bds = rbind(bds[range,], bds[range[1],])
  quilt.plot(coords, rep(1, nrow(coords)))
  polygon(x=bds[,1], y=bds[,2], col=rgb(.5, .1, .1, .5))
}