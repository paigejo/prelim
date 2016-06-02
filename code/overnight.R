# overnight code run script:
source("estimateVGrams.R")
source("simKriging.R")
varioResult = varioAnalysis(500)
krigingResult = unitSquareKriging(500)

source("estimateKfun.R")