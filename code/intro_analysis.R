#plot data
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

##### Make table 2 from paper:
loadRawData()
library(xtable)
tab = matrix(nrow=8, ncol=5)
tab[1,2] = "Levels ($mu$($g$ dray weight$)^{-1}$) for the following scales and years:)"
tab[2,2] = "Untransformed"
tab[2,4] = "Log-transformed"
tab[3,2] = 1997
tab[3,3] = 2000
tab[3,4] = 1997
tab[3,5] = 2000
tab[4,1] = "Number of locations"
tab[4,2] = length(log97)
tab[4,3] = length(log00)
tab[4,4] = length(log97)
tab[4,5] = length(log00)
tab[5,1] = "Mean"
tab[5,2] = mean(data97)
tab[5,3] = mean(data00)
tab[5,4] = mean(log97)
tab[5,5] = mean(log00)
tab[6,1] = "Standard deviation"
tab[6,2] = sd(data97)
tab[6,3] = sd(data00)
tab[6,4] = sd(log97) 
tab[6,5] = sd(log00)
tab[7,1] = "Minimum"
tab[7,2] = min(data97)
tab[7,3] = min(data00)
tab[7,4] = min(log97)
tab[7,5] = min(log00)
tab[8,1] = "Maximum"
tab[8,2] = max(data97)
tab[8,3] = max(data00)
tab[8,4] = max(log97)
tab[8,5] = max(log00)
xtable(tab, digits=2)
loadCorrectedData()

##### recreate Figure 4 from paper (ecdf of data)
pdf("ecdf.pdf", width=7, height=5)
parSet = par()
par(mar=c(5.1, 5, 4.1, 2.1))
plot(ecdf(log97), lty=1, pch=".", xlim=c(-.5, 2.5), main="Empirical CDF of data", 
     xlab=expression(italic(x)), ylab=expression(italic(hat(F)(x))))
plot(ecdf(log00), lty=1, pch=".", col="blue", add=TRUE)
legend("left", legend=c(1997, 2000), col=c("black", "blue"), lty=1)
suppressWarnings(par(parSet))
dev.off()
