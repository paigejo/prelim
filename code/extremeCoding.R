# coding tasks:
#
# 1: Make functions to generate 3 types of random point processes on unit square:
#      1) Completely random (uniform)
#      2) Preferential sample
#      3) Clustered sample
#    Note: mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
# 
# 2: Make function to estimate variograms from d=0 to 1
# 
# 3: Plot bias and standard deviation as in Fig. 2:
#      1) For bias, plot estimated bias +- 2 sd's vs dist for each sampling method
#      2) For SD, plot SD vs dist for each sampling method
#      3) Make sure to save the plots
#
# If there's more time, do the following tasks as well:
# 
# 4: Remember to push to git after each major step!
#
# 5: Miscellaneous:
#      1) Remove the 2 outliers in the 2000 data, replace with average of other points
#      2) Create log-transformed data
#
# 6: Fit variograms to each year's (log-transformed) data using maximum likelihood grid search
#    with traditional assumptions (assuming no preferential sampling), as in Fig. 5:
#      1) Make function for maximum likelihood or use function from research
#      2) Assume kappa=0.5, make grid on log scale centered on final values:
#           - mu97 = 1.44
#           - mu00 = 0.66
#           - sigma^2 = .45
#           - phi = exp(-1.163) = 0.3125471
#           - tau^2 = exp(-1.419)^2 = 0.05854263
#
# 7: Make a function that, for a given set of parameters, simulates S conditional on Y
#    as described on page 199.
#      - First draw S ~ MVN(0, Sigma(theta))
#      - Then draw iid Zi ~ N(0, tausq)
#      - Round n data locations to a big grid with N locations in the domain
#      - Define C (n x N) matrix given location of n obsevation in grid
#      - Define Sigma0 = C Sigma C' + tau^2 I
#    Then 
#      S + Sigma C' Sigma0'(y - mu + Z - C S)
#    is a draw from S|Y
#    NOTE: don't do any of this, try to use RPcirculant with RMmatern








