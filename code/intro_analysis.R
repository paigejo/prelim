#set working directory
#setwd("~/Google Drive/UW/coursework/2016_SPRING/prelim/code/")

# read in data
setwd("~/git/prelim/code")
source("loadData.R")
source("support_funs.R")

#plot data
library(fields)
library(latex2exp)
pdf("data_1997.pdf", width=5, height=6)
quilt.plot(coords97, data97, main=TeX("Lead Concentration ($\\mu g/g$), 1997"), 
           xlab="Easting (100 km)", ylab="Northing (100 km)", asp=1)
dev.off()

pdf("data_2000.pdf", width=5, height=6)
quilt.plot(coords00, data00, main=TeX("Lead Concentration ($\\mu g/g$), 2000"), 
           xlab="Easting (100 km)", ylab="Northing (100 km)", asp=1)
dev.off()

##### find outliers in data: points 32, 73.  Replace with mean of rest of data
text(coords00[,1], coords00[,2], 1:nrow(coords00))
outs = c(32, 73)
data00[outs] = mean(data00[-outs])

#plot resulting field
quilt.plot(coords00, data00, main=TeX("Lead Concentration ($\\mu g/g$), 2000"))

#plot both data sets simultaneously
pdf("sampling_locations.pdf", width=5, height=6)
plot(coords00[,1], coords00[,2], cex=.5, xlab="Easting (100 km)", asp=1, 
     ylab="Northing (100 km)", main="Lead Concentration Sampling locations")
points(coords97[,1], coords97[,2], pch=19, cex=.5)
legend("topleft", c("1997", "2000"), pch=c(19, 1))
dev.off()

##### make concave hull around data from prediction
library(alphahull)
source("~/git/prelim/code/makeAlphaHull.R")

# make the hull
allCoords = rbind(coords97, coords00)
hull = makePredHull(allCoords, alpha=5) #whoops, I guess I may has well have used a convex hull here...

#reverse order of hull polygon so it's compatible with owin() function
hull$hullPoly = hull$hullPoly[nrow(hull$hullPoly):1,]

# plot the hull and prediction points
plotAlphaHull(hull$hullPts, allCoords)

##### code for simulating log-Gaussian Doubly Stochastic Poisson Process
library(spatstat)
library(RandomFields)

#set parameters for point process (as in equation 5 from paper)
alpha = 0
beta = 2

#matern parameters:
#nu: smoothness
#var (or variance?): (marginal?) variance (I guess the sill)
#scale
#nugget?

#make window object by calculating range of data
range97 = apply(coords97, 2, range)
range00 = apply(coords00, 2, range)
ranges = apply(rbind(range97, range00), 2, range)
rangeX = ranges[,1]
rangeY = ranges[,2]
win = owin(poly=hull$hullPoly)

#calculate area of domain in (100km)^2 = 10,000 km^2
library(pracma)
area = polyarea(hull$hullPoly[,1], hull$hullPoly[,2])
lambda97 = nrow(coords97)/area

#####test simulation with MLEs from paper (p. 198)
#matern MLEs
nuMLE=.5 #implies exponential covariance
varianceMLE=.138
sillMLE=varianceMLE
scaleMLE=.313 #(in 100's of km)
nuggetMLE=.059

#other MLEs
muMLE=1.515 #overall mean (I think this is for the actual lead data rather than 
#the LGCP parameters?)
betaMLE=-2.198 #preferentiality parameter

#run simulation (input Gaussian Process has additional variance due to 
#being multiplied by beta in eq. 5 from paper)
#NOTE: I'm assuming I don't need a nugget to simulate a true random field, here
nsims=4
out = rLGCP("exp", muMLE, var=betaMLE^2*varianceMLE, scale=scaleMLE, 
            nsim=nsims, saveLambda=TRUE, win=win)
#exp(mu + sill/2) is new expectation
plot(out)
#plot(attr(out, "Lambda"))
#points(out)