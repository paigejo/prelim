# function to estimate empirical variogram of one PP
getVario = function(PP, GP, 
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
  GPCoords = attr(GP, "coords")
  findIndex = function(coords, len=128) {
    #convert coords to x and y grid indices
    coords = coords*(len-1) + 1
    xInd = coords[1]
    yInd = coords[2]
    return((yInd-1)*len + xInd)
  }
  inds = apply(roundCoords, 2, findIndex, len=length(xrange))
  
  #get PP marks from GP using the indices found above
  PPmarks = GP$variable1[inds]
  
}