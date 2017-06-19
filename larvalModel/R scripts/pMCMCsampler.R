require(boot)
require(deSolve)
require(ellipse)
require(coda)
require(mnormt)
require(emdbook)
library(odin)
library(stringr)
library(ggplot2)
library(lattice)
library(growthrates)
library(FME)
library(lubridate)
library(VGAM)
library(parallel)



##model step function - runs model in steps taking the liklihood of the subsequent observed data point
modStep3<-function(weightInput){
  
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  initialState<-as.numeric(dataInput[,c(1:4)])
  currentTime<-as.numeric(dataInput[,c(5)]*10)
  nextTime<-as.numeric(dataInput[,c(6)]*10)
  pr<-as.numeric(dataInput[,c(7:12)])
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],sf=pr[6])
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(currentTime, nextTime, length.out=nextTime-currentTime)))
  res2<-(simDat[simDat$step == nextTime,])
  return(as.data.frame(rbind(res2[,c(8:11)])))
}


## Log-Prior 
lprior <- function(parms) {
  uoEprior<-dnorm(parms[1],mean=0.035,sd=0.00485,log=T)
  uoLprior<-dnorm(parms[2],mean=0.035,sd=0.00485,log=T)
  uPprior<-dnorm(parms[3],mean=0.25,sd=0.0357,log=T)
  Yprior<-dnorm(parms[4],mean=13.06,sd=3.53,log=T)
  sfprior<-dunif(parms[6],min=0,max=100,log=T)
  priorSum<-as.vector(uoEprior+uoLprior+uPprior+Yprior+sfprior)
  return(priorSum)
}

##function that sums log-likelihood & log-prior inside MCMC sampler - use for when submitting to cluster at each p filter step
llikePriorLocal <- function(fit.params=NULL, ## parameters to fit
                            ref.params = mosParams(), ## reference parameters
                            obsDat=myDat) { ## observed data
  parms <- within(ref.params, { ## switch out old parameters in mos_params for new ones, keeping fixed parameters in place
    for(nm in names(fit.params)) assign(nm, as.numeric(fit.params[nm]))
    rm(nm)
  })
  print(fit.params)
pFilt(96,simx0,0,modStep3,dataLik2,garkiObs,pr=fit.params)+ lprior(fit.params)

}



## functions for logging and unlogging
logParms <- function(fit.params) {
  fit.params <- log(fit.params)
  names(fit.params) <- paste0('log',names(fit.params))
  return(fit.params)
}

unlogParms <- function(fit.params) {
  fit.params <- exp(fit.params)
  names(fit.params) <- sub('log','', names(fit.params))
  return(round(fit.params,8))
}



# set bounds on initial parameter guesses
initBounds <- data.frame(rbind( ## for initial conditions
  c(log(0.01),log(0.03)), ## uoE
  c(log(0.01),log(0.03)), ## uoL
  c(log(0.2),log(0.3)), ## uP
  c(log(10),log(15)), ## Y
  c(log(10),log(10)),##n
  c(log(1),log(10)),##sf
  c(log(0.1),log(1)))) ## p0

colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('loguoE','loguoL','loguP','logY','logn','logsf','logp0')
class(initBounds[,2]) <- class(initBounds[,1]) <- 'numeric'
initBounds


##  randomly select a value that is uniformly distributed between these bounds
initRand <- function(fit.params) {
  fit.params <- logParms(fit.params)
  tempnm <- names(fit.params)
  for(nm in tempnm) fit.params[nm] <- runif(1, min = initBounds[rownames(initBounds)==nm, 'lower'], 
                                            max =  initBounds[row.names(initBounds)==nm, 'upper'])
  return(unlogParms(fit.params))
}


## Sequential proposal function: Propose one parameter at a time
sequential.proposer <- function(sdProps) {
  nfitted <- length(sdProps)
  on <- 0
  return(list(sdProps = sdProps, type = 'sequential',
              fxn = function(current) {
                proposal <- current
                proposal[on + 1] <- proposal[on + 1] + rnorm(1, mean = 0, sd = sdProps[on + 1])
                on <<- (on+1) %% nfitted
                proposal
              }))
}





## Metropolis-Hastings Sampler
mcmcSampler <- function(init.params, ## initial parameter guess
                        randInit = T, ## if T then randomly sample initial parameters instead of above value
                        seed = 1, ## RNG seed
                        ref.params=mosParamsP(), ## fixed parameters
                        obsDat = myDat, ## data
                        proposer = sequential.proposer(sdProps=sdProps), ## proposal distribution
                        niter = 100, ## MCMC iterations
                        nburn = 0, ## iterations to automatically burn
                        adaptiveMCMC = F, ## adapt proposal distribution?
                        startAdapt = 150, ## start adapting at what iteration?
                        adptBurn = 200, ## ignore first so many iterations for adapting posterior
                        verbose=0, ## if >2 browses, if >1 prints progress
                        tell = 100) { ## how often to print progress
  if(verbose>2) browser()
  if(randInit) init.params <- initRand(init.params)
  current.params <- init.params
  nfitted <- length(current.params) ## number fitted parameters
  vv <- 2 ## mcmc iteration (started at 1 so we're already on 2
  accept <- 0 ## initialize proportion of iterations accepted
  ## Calculate log(likelihood X prior) for first value
  curVal <- llikePriorLocal(current.params, ref.params = ref.params, obsDat=obsDat)
  print(curVal)
  ## Initialize matrix to store MCMC chain
  out <- matrix(NA, nr = niter, nc=length(current.params)+1)
  out[1,] <- c(current.params, ll = -curVal) ## add first value
  colnames(out) <- c(names(current.params), 'll') ## name columns
  ## Store original covariance matrix
  if(proposer$type=='block') originalCovar <- get('covar', envir = environment(proposer$fxn)) 
  while(vv <= niter) {
    if ((verbose > 1) || (verbose && (vv%%tell == 0))) print(paste("on iteration",vv,"of", niter + 1))
    ## Adaptive MCMC: adapt covariance every 50 iterations (don't
    ## do it more often because it adds to coputational burden.
    if(adaptiveMCMC & proposer$type=='block' & vv > startAdapt & vv %% 50 == 0) {
      adptBurn <- min((startAdapt-50), adptBurn)
      ## Below equation gives ideal covariance-variance matrix based on posterior
      adaptedCovar <- 2.38^2 / nfitted * cov.wt(log(out[adptBurn:(vv-1),1:nfitted]))$cov
      adaptedCovar <- adaptedCovar*.95 + originalCovar*.05 ## 95% adapted & 5% original
      rownames(adaptedCovar) <- colnames(adaptedCovar) <- names(current.params)
      assign('covar', adaptedCovar, envir = environment(proposer$fxn))
    }
    proposal <- proposer$fxn(logParms(current.params))
    proposal <- unlogParms(proposal)
    propVal <- llikePriorLocal(proposal, ref.params = ref.params, obsDat=obsDat)
    lmh <- propVal - curVal ## likelihood ratio = log likelihood difference
    if (is.na(lmh)) { ## if NA, print informative info but don't accept it
     # print(list(lmh=lmh, proposal=exp(proposal), vv=vv, seed=seed))
    } else { ## if it's not NA then do acception/rejection algorithm
      if (verbose > 1) print( c(lmh=lmh, propVal=propVal, curVal=curVal) )
      ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
      if ( (lmh >= 0) | (runif(1,0,1) <= exp(lmh)) ) {
        current.params <- proposal
        if (vv>nburn) accept <- accept + 1 ## only track acceptance after burn-in
        curVal <- propVal
      }
    }
    out[vv, ] <- c(current.params, ll=curVal)
    vv <- vv+1
    aratio <- accept/((vv-nburn))
  }
  colnames(out) <- c(names(current.params), 'll')
  samp <- as.mcmc(out[1:nrow(out)>(nburn+1),])
  return(list(ref.params=ref.params
              , seed = seed
              , init.params = init.params
              , aratio = aratio
              , samp = samp
  ))
}



set.seed(44)
runX200 <- mcmcSampler(init.params = c(uoE=0.03,uoL=0.03,uP=0.25,Y=13,n=15,sf=19,p0=0.1)
                    , seed = 1
                    ,nburn=2000
                    , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.0,0.1,0.01))
                    , randInit = T
                    ,verbose=2
                    , niter = 50000)

write.table(runX200,"C:\\res.csv")


test.params = c(uoE=0.001729949,uoL=0.001400767,uP=0.23779609,Y=13.24491630,n=1,sf=1,p0=0.1)
set.seed(10)
res3<-NULL
for (i in 1:10){
  ss<-pFilt(50,simx0,0,modStep3,dataLik2,garkiObs,pr=test.params)#+ lprior(test.params)
  res3<-rbind(ss,res3)
  print(ss)

}


res4<-NULL
for (i in 1:100){
  
  ss<-pFilt(5,simx0,0,modStep2,dataLik2,garkiObs,pr=test.params)+ lprior(test.params)
  res4<-rbind(ss,res4)
  print(ss)
  
}