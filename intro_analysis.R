#set working directory
setwd("~/Google Drive/UW/coursework/2016_SPRING/prelim/code/")

# read in data
source("loadData.R")

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
