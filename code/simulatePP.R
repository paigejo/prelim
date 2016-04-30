# genMaternGP generates a Matern covariance GP on the unit square
# with the parameters given the in text
genMaternGP = function(nsim=1, nx=128, ny=128) {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  
  mu = 4
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 0
  tausq = 0
  
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
genUnifPP = function(numSamples=1, nPts = 100) {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  
  mu = 4
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 0
  tausq = 0
  
  sims = list()
  for(i in 1:numSamples) {
    sims[[i]] = rpoispp(lambda=100)
  }
  
  return(sims)
}

# get a a clustered point process using log-Gaussian Cox Process
# model
genPreferentialPP = function(numSamples=1, nPts = 100) {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  # exp(mu + sigma^2/2)
  mu = 4
  sigmasq = 1.5
  phi= 0.15
  kappa=1
  beta = 2
  tausq = 0
  varEst = beta^2*sigmasq
  overallMuEst = log(nPts) - varEst/2
  
  #now generate point processes clustered based on these GPs
  #make sure to save random rate Lambda and the GPs themselves
  sims = list()
  Lambdas = list()
  GPs = list()
  for(i in 1:numSamples) {
    sim = rLGCP("matern", overallMuEst, var=varEst, scale=phi, 
                      nu=kappa, nsim=1, drop=TRUE, saveLambda=TRUE)
    sims[[i]] = sim
    Lambdas[[i]] = attr(sim, "Lambda")
    GPs[[i]] = log(Lambdas[[i]])
  }
  
  return(list(sims=sims, Lambdas=Lambdas, GPs=GPs))
}

#generate clustered point processes (preferentially sampled from 
# a different GP)
genClusterPP = function(numSamples=1, nPts = 100) {
  # generate pereferential PPs, but associate with different
  # GPs to make clustered, non-preferential sampling
  
  preferentialPPs = genPreferentialPP(numSamples, nPts)
  newGPs = genMaternGP(nsim=numSamples)
  preferentialPPs$GPs = newGPs
  preferentialPPs$Lambda = NULL
  
  return(preferentialPPs)
}