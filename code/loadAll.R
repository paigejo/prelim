
# load in all my functions, used libraries, and the corrected data
setwd("~/git/prelim/code")
library(Matrix)
library(fields)
library(spatstat)
library(RandomFields)
library(snow)
#library(alphahull) #not really necessary, but for creating spatial domain polygon
library(pracma) #for calculating domain area
source("loadData.R")
source("support_funs.R")
#source("makeAlphaHull.R")
source("simulatePP.R")
source("simKriging.R")
source("estimateVGrams.R")
source("ML.R")
#source("estimateKfun.R") this is a script, not a set of functions
loadCorrectedData()