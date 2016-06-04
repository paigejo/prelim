library(RandomFields)
library(spatstat)

# genMaternGP generates a Matern covariance GP on the unit square
# with the parameters given the in text
genMaternGP = function(nsim=1, nx=resGP, ny=resGP, mu=4, sigmasq=1.5, phi=.15, kappa=1, asList=TRUE, coords=NULL, method="instrinsic") {
  #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
  
  # use RFmatern and RFsimulate
  # obj = RMmatern(nu=kappa, var=sigmasq, scale=phi)
  obj = RMwhittle(nu=kappa, var=sigmasq, scale=phi)
  
  if(is.null(coords)) {
    coordsSet=TRUE
    xs = seq(0, 1, length=nx)
    ys = seq(0, 1, length=ny)
    coords = make.surface.grid(list(x=xs, y=ys))
  }
  
  if(method == "instrinsic")
    sims = as.matrix(RFsimulate(RPintrinsic(obj), x=coords[,1], y=coords[,2], n=nsim)) + mu
  else if(method == "circulant")
    sims = as.matrix(RFsimulate(RPcirculant(obj), x=coords[,1], y=coords[,2], n=nsim)) + mu
  else if(method == "cutoff")
    sims = as.matrix(RFsimulate(RPcutoff(obj), x=coords[,1], y=coords[,2], n=nsim)) + mu
  
  if(asList) {
    listSims = list()
    for(i in 1:nsim) {
      # GPs computed from rLGCP must be transposed in VG computation, but these don't
      if(coordsSet)
        listSims[[i]] = list(v=cbind(sims)[,i], xcol=xs, yrow=ys, noTranspose=TRUE)
      else
        listSims[[i]] = list(v=cbind(sims)[,i], noTranspose=TRUE)
    }
    return(listSims)
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
  tmpMu = log(npts) - varEst/2
  tmpMu = tmpMu + log(3) #overgen points
  
  #now generate point processes clustered based on these GPs
  #make sure to save random rate Lambda and the GPs themselves
  sims = list()
  Lambdas = list()
  GPs = list()
  test = c()
  testLambdas = list()
  i = 1
  while(length(sims) < numSamples) {
    # generate PP simulation
#     sim = rLGCP("matern", mu, var=varEst, scale=phi, dimyx=xyRes, 
#                       nu=kappa, nsim=1, drop=TRUE, saveLambda=TRUE)
    sim = rLGCP("whittle", tmpMu, var=varEst, scale=phi, dimyx=xyRes, 
                nu=kappa, nsim=1, drop=TRUE, saveLambda=TRUE)
    n = sim$n
    if(n < npts)
      next
    
    #perform independent thinning
    samp = sample(1:sim$n, npts)
    sim$n=npts
    sim$x = sim$x[samp]
    sim$y = sim$y[samp]
    
    #save results
    testLambdas[[i]] = attr(sim, "Lambda")
    GPs[[i]] = (log(attr(sim, "Lambda")) + beta*mu - tmpMu)/beta #transform Lambda to S
    test[i] = var(c(log(attr(sim, "Lambda")$v)))
    attr(sim, "Lambda") = attr(sim, "Lambda")*(npts/n) #modify Lambda to account for random thinning
    Lambdas[[i]] = attr(sim, "Lambda")
    sims[[i]] = sim
    i = i + 1
  }
  
  return(list(sims=sims, Lambdas=Lambdas, GPs=GPs, test=test, testLambdas=testLambdas, tmpMu=tmpMu))
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
  preferentialPPs$GPs = NULL
  preferentialPPs$Lambda = NULL
  
  return(preferentialPPs)
}



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
##### My versions of PP generating functions
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################




# get a a clustered point process using log-Gaussian Cox Process
# model
genPreferentialPP2 = function(numSamples=1, xyRes=resGP, 
                             mu=4, sigmasq=1.5, phi=.15, kappa=1, beta=2, tausq=0, 
                             npts=100, GPcoords=NULL, method="instrinsic", replace=TRUE) {
  nx=xyRes
  ny=xyRes
  
  ##### generate GPs
  GPs = genMaternGP(nsim=numSamples, nx=resGP, ny=resGP, mu=mu, sigmasq=sigmasq, phi=phi, kappa=kappa, coords=GPcoords, method=method)
  
  #add nugget noise to GPs
  addNugget = function(GP) {
    GP$v = GP$v + rnorm(length(GP$v), sd=sqrt(tausq))
    return(GP)
  }
  GPs = lapply(GPs, addNugget)
  
  #compute GP coordinates
  if(is.null(GPcoords)) {
    xs = seq(0, 1, length=nx)
    ys = seq(0, 1, length=ny)
    GPcoords = make.surface.grid(list(x=xs, y=ys))
  }
  
  ##### convert to Lambdas
  GP2Lambda = function(GP) {
    return(exp(beta*GP$v))
  }
  Lambdas = lapply(GPs, GP2Lambda)
  
  ##### generate PPs using Lambdas
  Lambda2PP = function(i) {
    Lambda = Lambdas[[i]]
    Lamda = Lambda/sum(Lambda)
    inds = sample(1:length(Lambda), size=npts, prob=Lambda, replace=replace)
    PPcoords = GPcoords[inds,]
    return(list(x=PPcoords[,1], y=PPcoords[,2], n=nrow(PPcoords), marks=GPs[[i]][inds]))
  }
  PPs = lapply(1:numSamples, Lambda2PP)
  
  return(list(sims=PPs, Lambdas=Lambdas, GPs=GPs))
}

# get a a clustered point process using log-Gaussian Cox Process
# model
genClusterPP2 = function(numSamples=1, xyRes=resGP, 
                              mu=4, sigmasq=1.5, phi=.15, kappa=1, beta=2, tausq=0, 
                              npts=100, method="instrinsic") {
  
  # first generate preferentially
  prefPPs = genPreferentialPP2(numSamples, xyRes, mu, sigmasq, phi, kappa, beta, tausq, npts, method=method)
  
  prefPPs$Lambdas=NULL
  prefPPs$GPs=NULL
  prefPPs$noTranspose = TRUE # GPs computed from rLGCP must be transposed in VG computation, these don't
  
  return(prefPPs)
}