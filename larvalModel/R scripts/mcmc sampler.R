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
#setwd("C:\\Users\\Aaron\\Documents\\Imperial\\larval model\\Data")
#setwd("C:\\ImperialMalaria\\larval model\\Data")




## Log-Prior 
lprior <- function(parms) with(parms, {
  uoEprior<-dnorm(parms$uoE,mean=0.035,sd=0.00485,log=T)
  uoLprior<-dnorm(parms$uoL,mean=0.035,sd=0.00485,log=T)
  uPprior<-dnorm(parms$uP,mean=0.25,sd=0.0357,log=T)
  Yprior<-dnorm(parms$Y,mean=13.06,sd=3.53,log=T)
  sfprior<-dunif(parms$sf,min=0,max=100,log=T)
  priorSum<-(uoEprior+uoLprior+uPprior+Yprior+sfprior)
  return(priorSum)
})


##function that sums log-likelihood & log-prior inside MCMC sampler
llikePrior <- function(fit.params=NULL, ## parameters to fit
                       ref.params = mosParamsP(), ## reference parameters
                       obsDat=myDat) { ## observed data
    parms <- within(ref.params, { ## switch out old parameters in mos_params for new ones, keeping fixed parameters in place
    for(nm in names(fit.params)) assign(nm, as.numeric(fit.params[nm]))
    rm(nm)
  })
    
    print(fit.params)
##apply particle filter function on cluster
   particleTemp1<- obj$enqueue(parRunSAll(runs=1,particles=400,theta=fit.params),name="particle MCMC")

   
   pt1<-particleTemp1$wait(Inf)


   px<-c(pt1)
   
   mean(as.numeric(unlist(px))) + lprior(parms)
}


##function that sums log-likelihood & log-prior inside MCMC sampler - use for when submitting to cluster at each p filter step
llikePriorBigClust <- function(fit.params=NULL, ## parameters to fit
                       ref.params = mosParams(), ## reference parameters
                       obsDat=myDat) { ## observed data
  parms <- within(ref.params, { ## switch out old parameters in mos_params for new ones, keeping fixed parameters in place
    for(nm in names(fit.params)) assign(nm, as.numeric(fit.params[nm]))
    rm(nm)
  })
  print(fit.params)
  
  particleFilterMCMC(larvalModP, ###particle filter
                     theta=fit.params,
                     init.state = init.state,
                     data = garkiObs,
                     nParticles = 20) + lprior(parms)
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
  
  pfMLLik(50,simx0,0,modStep2,dataLik2,garkiObs,pr=fit.params)+ lprior(parms)
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
  c(log(0.01),log(0.06)), ## uoE
c(log(0.01),log(0.06)), ## uoL
c(log(0.1),log(0.99)), ## uP
c(log(1),log(20)), ## Y
c(log(7),log(7)),##n
c(log(1),log(50)))) ## sf

colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('loguoE','loguoL','loguP','logY','logn','logsf')
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
                        ref.params=mosParams(), ## fixed parameters
                        obsDat = myDat, ## data
                        proposer = sequential.proposer(sdProps=sdProps), ## proposal distribution
                        niter = 100, ## MCMC iterations
                        nburn = 1, ## iterations to automatically burn
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
   # obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",cores = 16,home="//fi--san02/homes/alm210",credentials = creds))
    
    if ((verbose > 1) || (verbose && (vv%%tell == 0))) print(paste("on iteration",vv,"of", niter + 1))
    ## Adaptive MCMC: adapt covariance every 50 iterations (don't
    ## do it more often because it adds to coputational burden.
    if(adaptiveMCMC & proposer$type=='block' & vv > startAdapt & vv %% 50 == 0) {
      adptBurn <- min((startAdapt-50), adptBurn)
      ## Below equation gives ideal covariance-variance matrix based on posterior
      adaptedCovar <- 2.38^2 / nfitted * cov.wt(log(out[adptBurn:(vv-1),1:nfitted]))$cov
      ## Take a weighted average of the original & the empirical cov-var matrices to ensure
      ## that we never let the matrix collapse to zero (ie if the empirical one is zero
      ## because we haven't accepted anything yet)
      adaptedCovar <- adaptedCovar*.95 + originalCovar*.05 ## 95% adapted & 5% original
      rownames(adaptedCovar) <- colnames(adaptedCovar) <- names(current.params)
      assign('covar', adaptedCovar, envir = environment(proposer$fxn))
    }
    proposal <- proposer$fxn(logParms(current.params))
    proposal <- unlogParms(proposal)
    propVal <- llikePriorLocal(proposal, ref.params = ref.params, obsDat=obsDat)
    lmh <- propVal - curVal ## likelihood ratio = log likelihood difference
    print(propVal)
    print("curval")
    print(curVal)
    if (is.na(lmh)) { ## if NA, print informative info but don't accept it
     # print(list(lmh=lmh, proposal=exp(proposal), vv=vv, seed=seed))
    } else { ## if it's not NA then do acception/rejection algorithm
      if (verbose > 1) print( c(lmh=lmh, propVal=propVal) )
      ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
      if ( (lmh >= 0) | (runif(1,0,1) <= exp(lmh)) ) {
        current.params <- proposal
        if (vv>nburn) accept <- accept + 1 ## only track acceptance after burn-in
        curVal <- propVal
        
        out[vv, ] <- c(current.params, ll=curVal)
        vv <- vv+1
        aratio <- accept/((vv-nburn))
        
        colnames(out) <- c(names(current.params), 'll')
        outT<<-out
        
        samp <- as.mcmc(out[1:nrow(out)>(nburn+1),])

        resTemp<<-(list(ref.params=ref.params
                       , seed = seed
                       , init.params = init.params
                       , aratio = aratio
                       , samp = samp
        ))
        
      }
    }

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



set.seed(10)
run <- mcmcSampler(init.params = c(uoE=0.02507989,uoL=0.01657317,uP=0.28016918,Y=7.37661308,n=2.367782,sf=2.28966396)
                        , seed = 1
                        ,nburn=1
                        , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.1))
                        , randInit = T
                        , niter = 100000)






par(oma = c(0,0,2,0), bty='n', 'ps' = 18)
plot(run$samp)
mtext('mcmc', side = 3, outer = T, line = 0)



class(run$samp)
median(as.vector(run$samp[,c('uoE')][[1]]))
median(as.vector(run$samp[,c('uoL')][[1]]))
median(as.vector(run$samp[,c('uP')][[1]]))
median(as.vector(run$samp[,c('Y')][[1]]))
median(as.vector(run$samp[,c('sf')][[1]]))
median(as.vector(run$samp[,c('n')][[1]]))
run$aratio
