
##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    pMCMC functions                                                                             #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

##function that sums log-likelihood & log-prior inside MCMC sampler
llfnc <- function(fitParams=NULL, ## parameters to fit
                            particles) { ## observed data
  #run particle filter and return ll for global and local parameters
  ##put each village data into a column, then loop through each column doing the below
  pX<-NULL
  for (i in 1:4){
  globalParms<-fitParams[c(1:8)]
  globalParms[8]<-fitParams[c(i+7)]#i+6 to fit the scaling factor specific for each village
  p1Parms<-globalParms
garkDat<-garkiObsX[,c(1,i+1)]#i+1 as first column is "time"

p1<-pFilt(particles,simx0,modStep3,dataLikFunc,garkDat,pr=p1Parms)+lprior(p1Parms)

pX<-rbind(pX,p1)
}
return(mean(pX))
}


## Log-Prior function
lprior <- function(parms) {
  uoEprior<-dnorm(parms[1],mean=0.035,sd=0.00485,log=T)
  uoLprior<-dnorm(parms[2],mean=0.035,sd=0.00485,log=T)
  uPprior<-dnorm(parms[3],mean=0.25,sd=0.0357,log=T)
  Yprior<-dnorm(parms[4],mean=13.06,sd=3.53,log=T)
  #sfprior<-dunif(parms[6],min=0,max=100,log=T)
  p0prior<-dnorm(parms[6],mean=0.5,sd=0.015,log=T)
  priorSum<-as.vector(uoEprior+uoLprior+uPprior+Yprior+p0prior)
  return(priorSum)
}


# set bounds on initial parameter guesses
initBounds <- data.frame(rbind( ## for initial conditions
  c(0.03,0.04), ## uoE
  c(0.03,0.04), ## uoL
  c(0.2,0.3), ## uP
  c(8,15), ## Y
  c(25,25),##n
  c(0.4,0.6),## p0
  c(0.01,0.9),##o
  c(2,6),
  c(2,6),
  c(2,6),
  c(2,6)))##sf


colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('uoE','uoL','uP','Y','n','p0','logo','sf1','sf2','sf3','sf4')
class(initBounds[,2]) <- class(initBounds[,1]) <- 'numeric'
initBounds


##  randomly select a value that is uniformly distributed between these bounds
initRand <- function(fitParams) {
  fitParams <- logParms(fitParams)
  tempnm <- names(fitParams)
  for(nm in tempnm) fitParams[nm] <- runif(1, min = initBounds[rownames(initBounds)==nm, 'lower'], 
                                            max =  initBounds[row.names(initBounds)==nm, 'upper'])
  return(unlogParms(fitParams))
}



## Sequential proposal function: Propose one parameter at a time
sequential.proposer <- function(sdProps) {
  nfitted <- length(sdProps)
  on <- 0
  return(list(sdProps = sdProps, type = 'sequential',
              fxn = function(current) {
                proposal <- current
                if(sdProps[on+1]!=0){##adds clause that if sd is 0 then skip the paramter
                  propVal<-proposal[on + 1] +rnorm(1, mean = 0, sd = sdProps[on + 1])
                proposal[on + 1] <- if(propVal<0) 0 else propVal
                on <<- (on+1) %% nfitted
                return(proposal)}
                else  proposal[on + 2] <- proposal[on + 2] + rnorm(1, mean = 0, sd = sdProps[on + 2])
                on <<- (on+2) %% nfitted
                return(proposal)
              }))
}




##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    Metropolis-Hastings MCMC Sampler                                                            #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

mcmcSampler <- function(initParams, ## initial parameter guess
                        randInit = T, ## if T then randomly sample initial parameters instead of above value
                        obsDat = myDat, ## data
                        proposer = sequential.proposer(sdProps=sdProps), ## proposal distribution
                        niter = 100, ## MCMC iterations
                        particles =100,##number of particles for particle filter
                        nburn = 0, ## iterations to automatically burn
                        monitoring=0, ## if >2 browses, if >1 prints progress
                        tell = 100) { ## how often to print progress
  if(monitoring>2) browser()
  if(randInit) initParams <- initRand(initParams)
  currentParams <- initParams
  nfitted <- length(currentParams) ## number fitted parameters
  iter <- 2 ## mcmc iteration (started at 1 so we're already on 2
  accept <- 0 ## initialize proportion of iterations accepted
  ## Calculate log(likelihood X prior) for first value
  curVal <- llfnc(currentParams, particles=particles)
  print(curVal)
  ## Initialize matrix to store MCMC chain
  out <- matrix(NA, nr = niter, nc=length(currentParams)+1)
  out[1,] <- c(currentParams, ll = -curVal) ## add first value
  colnames(out) <- c(names(currentParams), 'll') ## name columns
  ## Store original covariance matrix
  while(iter <= niter) {
    if ((monitoring > 1) || (monitoring && (iter%%tell == 0))) print(paste("on iteration",iter,"of", niter + 1))
    proposal <- proposer$fxn(currentParams)
    print(proposal)
    propVal <- llfnc(proposal, particles=particles)
    lmh <- propVal - curVal ## likelihood ratio = log likelihood difference
    if (is.na(lmh)) { ## if NA, do not accept
    } else { ## if it's not NA then do acception/rejection algorithm
      if (monitoring > 1) print( c(lmh=lmh, propVal=propVal, curVal=curVal))
      if ( (lmh >= 0) | (runif(1,0,1) <= exp(lmh)) ) {
        currentParams <- proposal
        if (iter>nburn) accept <- accept + 1 ## only track acceptance after burn-in
        curVal <- propVal
      }
    }
    out[iter, ] <- c(currentParams, ll=curVal)
    iter <- iter+1
    aratio <- accept/((iter-nburn))
  }
  colnames(out) <- c(names(currentParams), 'll')
  results <- as.mcmc(out[1:nrow(out)>(nburn+1),])
  return(list(initParams = initParams
              , aratio = aratio
              , results = results
  ))
}



