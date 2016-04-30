#set working directory
#setwd("~/Google Drive/UW/coursework/2016_SPRING/prelim/code/")
setwd("~/git/prelim/code")

# read in data
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

# find outliers in data: points 32, 73
text(coords00[,1], coords00[,2], 1:nrow(coords00))
outs = c(32, 73)

#replace outlier values with mean of the rest of the data
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


##### K function estimation
# lgcp.estK, use Kest for estimating Ripley's reduced second moment K function 
#from given point pattern

#get empirical K function from 1997 lead data
#NOTE: Ripley isotropic correction is used, since it is claimed to be best for 
#polygonal boundary windows in Kest documentation
xs = coords97[,1]
ys = coords97[,2]
zs = data97
ppp97 = ppp(xs, ys, marks=zs, window=win)
KEmp97Test = Kest(ppp97, correction="isotropic") #as in the paper, only go out 25km
plot(KEmp97Test)
env97Test = envelope(ppp97, correction="isotropic", fix.n=TRUE)
plot(env97Test)

KEmp97 = Kest(ppp97, rmax=.25) #as in the paper, only go out 25km
plot(KEmp97)
env97 = envelope(ppp97, rmax=.25, fix.n=TRUE)
plot(env97)

rs = seq(.001, .25, length=100)
Krs = sapply(rs, K)
lines(rs, K(rs), col="green")

# since 2000 data is in lattice, K function looks terrible
# xs = coords00[,1]
# ys = coords00[,2]
# zs = data00
# ppp00 = ppp(xs, ys, marks=zs, window=win)
# KEmp00 = Kest(ppp00, rmax=.25) #as in the paper, only go out 25km
# plot(KEmp00)
# env00 = envelope(ppp00, rmax=.25)
# plot(env00)


# NOTE: Kinhom is like Kest but for anisotropic data.  Why wasn't this used in paper?
# maybe the differences between estimate and predicted K func in paper indicated 
# clustering

#in order to simulate log gaussian cox process, we increase alpha parameter from
#paper, which is equivalent to increasing mu in rLGCP function, so as to generate 
#a large number of points.  We then randomly thin until we reach correct number of 
#points.  A paper (http://people.math.aau.dk/~jm/teheran.pdf) notes on page 3 
#that this process results in a LGCP with same beta, just different alpha.

simExp = expression(simFun97())
testEnv97 = envelope(ppp97, fun=Kinhom, rmax=.25, simulate=simExp, clamp=FALSE, 
                     savefuns=TRUE, savepatterns = TRUE)
legendargs = list(x="topleft", 
                  legend=expression(italic(hat(K)[obs](r)), italic(bar(K)[sim](r)), 
                                    italic(K[theo](r))), 
                  lty=c(1, 2, 1), col=c("black", "red", "blue"), y.intersp=1.25)
# legendargs = list(y.intersp=.3)
pdf("envelopeInhom.pdf", width=7, height=5)
plot(testEnv97, legendargs=legendargs, xlim=c(0, .25), ylim=c(0, .3), 
     main=expression("Estimated, simulated, and theoretical" ~ italic(K) ~ 
                       "functions (not assuming homogeneity)"))
rs = seq(.001, .25, length=100)
Krs = sapply(rs, K)
lines(rs, Krs, col="blue")
dev.off()

simKs97 = attr(testEnv97, "simfuns")
class(simKs97) = "data.frame"
r = simKs97[,1]
simKs97 = simKs97[,-1]
simKSD97 = apply(test, 1, sd)
Krs = sapply(r, K)
numSDs = qnorm(.975) #for 95% confidence interval
uncertainties = simKSD97*numSDs
minRange = min(Krs - uncertainties)
minRange = min(c(minRange, testEnv97$obs))
maxRange = max(Krs + uncertainties)
maxRange = max(c(maxRange, testEnv97$obs))
polyX = c(r, rev(r), r[1])
polyY = c(Krs + uncertainties, rev(Krs - uncertainties), Krs[1] + uncertainties[1])
plot(r, Krs, type="n", lwd=2, col="blue", 
     main=expression("Observed and theoretical" ~ italic(K) ~ 
                       "functions"), 
     xlim=c(0, .25), ylim=c(0, .4), xlab=expression(italic(r)), 
     ylab=expression(italic(K(r))))
polygon(polyX, polyY, col=rgb(.4,.4,.4,.5), border=rgb(.5,.5,.5,.5))
lines(r, Krs, lwd=2, col=rgb(.4,.4,.4), lty=2)
# plot(r, Krs, type="l", ylim=c(minRange, maxRange), lwd=2, col="blue", 
#      main=expression("Observed and theoretical" ~ italic(K) ~ 
#                      "functions"))
# lines(r, Krs - uncertainties, lty=2, col=rgb(.4,.4,.4))
# lines(r, Krs + uncertainties, lty=2, col=rgb(.4,.4,.4))
lines(r, testEnv97$obs, lwd=2)
legend(x="topleft", lty=c(1, 2), col=c("black", rgb(.4,.4,.4)), y.intersp=1, 
       legend=expression(italic(hat(K)[obs](r)), italic(K[theo](r))), 
       lwd=2)

#plot of just theoretical K function:
r = seq(0, .25, length=100)
Krs = K(r)

pdf("Ktheo.pdf", width=7, height=5)
plot(r, Krs, main="Theoretical K function using MLEs", xlab=expression(italic(r)), 
     ylab=expression(italic(K(r))), type="l", lwd=2, col="blue")
dev.off()

#questions:
#1) K function: what is it?
#2) Why is my K function not the same as the one in the paper?
#3) How is the test statistic given on page 205 turned into a p-value?
#4) How does one typically fit a quadratic to a log-likelihood surface?
#5) Can I just take an independent thinning of a simulated LGCP?
