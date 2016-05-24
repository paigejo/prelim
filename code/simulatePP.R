library(RandomFields)
library(spatstat)

# genMaternGP generates a Matern covariance GP on the unit square
# with the parameters given the in text
genMaternGP = function(nsim=1, nx=resGP, ny=resGP, mu=4, sigmasq=1.5, phi=.15, kappa=1, beta=2, tausq=0) {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  
  # use RFmatern and RFsimulate
  obj = RMmatern(nu=kappa, var=sigmasq, scale=phi)
  
  xs = seq(0, 1, length=nx)
  ys = seq(0, 1, length=ny)
  coords = make.surface.grid(list(x=xs, y=ys))
  
  sims = list()
  for(i in 1:nsim) {
    sims = c(sims, RFsimulate(obj, x=coords[,1], y=coords[,2]))
  }
  
  return(sims)
}


# get a uniformly sampled point process with expected number of points
# equal to nPts
genUnifPP = function(numSamples=1, npts=100) {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  
  mu = 4
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 0
  tausq = 0
  
  #overallMu = exp(mu + sigmasq/2)
  overallMu = 1000 #overgenerate points
  
  sims = list()
  i = 1
  while(length(sims) < numSamples) {
    # generate PP
    sim = rpoispp(lambda=overallMu)
    n = sim$n
    if(n < npts)
      next
    
    #perform iid thinning to get 100 points
    samp = sample(1:n, npts)
    sim$x = sim$x[samp]
    sim$y = sim$y[samp]
    sim$n = npts
    
    #save results
    sims[[i]] = sim
    i = i+1
  }
  
  return(sims)
}

# get a a clustered point process using log-Gaussian Cox Process
# model
genPreferentialPP = function(numSamples=1, xyRes = resGP, 
                             mu=4, sigmasq=1.5, phi=.15, kappa=1, beta=2, tausq=0, 
                             npts=100) {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  # exp(mu + sigma^2/2)
  varEst = beta^2*sigmasq
#   overallMuEst = log(npts) - varEst/2
#   overallMuEst = overallMuEst + log(3) #overgen points
  
  #now generate point processes clustered based on these GPs
  #make sure to save random rate Lambda and the GPs themselves
  sims = list()
  Lambdas = list()
  GPs = list()
  i = 1
  while(length(sims) < numSamples) {
    # generate PP simulation
    sim = rLGCP("matern", mu, var=varEst, scale=phi, dimyx=xyRes, 
                      nu=kappa, nsim=1, drop=TRUE, saveLambda=TRUE)
    n = sim$n
    if(n < npts)
      next
    
    #perform independent thinning
    samp = sample(1:sim$n, npts)
    sim$n=npts
    sim$x = sim$x[samp]
    sim$y = sim$y[samp]
    attr(sim, "Lambda") = attr(sim, "Lambda")*(npts/n)
    
    #save results
    sims[[i]] = sim
    Lambdas[[i]] = attr(sim, "Lambda")
    GPs[[i]] = (log(Lambdas[[i]]) - mu)*(1/beta) + mu #transform Lambda to GP with mean mu, sd sigma
    i = i + 1
  }
  
  return(list(sims=sims, Lambdas=Lambdas, GPs=GPs))
}

#generate clustered point processes (preferentially sampled from 
# a different GP)
genClusterPP = function(numSamples=1, nx=resGP, ny=resGP, 
                        mu=4, sigmasq=1.5, phi=.15, kappa=1, beta=2, tausq=0, 
                        npts=100) {
  # generate pereferential PPs, but associate with different
  # GPs to make clustered, non-preferential sampling
  
  preferentialPPs = genPreferentialPP(numSamples, mu=mu, sigmasq=sigmasq, phi=phi, 
                                      kappa=kappa, beta=beta, tausq=tausq, npts=npts)
  newGPs = genMaternGP(nsim=numSamples, nx=nx, ny=nx, mu=mu, sigmasq=sigmasq, phi=phi, 
                       kappa=kappa, beta=beta, tausq=tausq)
  preferentialPPs$GPs = newGPs
  preferentialPPs$Lambda = NULL
  
  return(preferentialPPs)
}