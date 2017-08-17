
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
  globalParms<-fitParams[c(1:10)]
  globalParms[10]<-fitParams[c(i+9)]#i+6 to fit the scaling factor specific for each village
garkDat<-garkiObsX[,c(1,i+1)]#i+1 as first column is "time"
garkDat<-na.omit(garkDat)
rFc<-as.numeric(substring(names(garkDat)[2],1,1))#which rainfall cluster to use
#run particle filter for village
p1<-pFilt(particles,iState,modStep3,dataLikFunc,garkDat,pr=globalParms,rFclust=rFc)+lprior(globalParms)

pX<-rbind(pX,p1)
}
return(sum(pX))
}




## Log-Prior function
lprior <- function(parms) {
  uoEprior<-dnorm(parms[1],mean=0.035,sd=0.007,log=T)
  uoLprior<-dnorm(parms[2],mean=0.035,sd=0.007,log=T)
  uPprior<-dnorm(parms[3],mean=0.25,sd=0.0457,log=T)
  Yprior<-dnorm(parms[4],mean=13.06,sd=4.53,log=T)
  p0prior<-dunif(parms[6],min=0,max=1,log=T)
  FpPrior<-dunif(parms[8],min=0,max=1,log=T)
  oprior<-dunif(parms[7],min=0.1,max=1000)
  priorSum<-as.vector(uoEprior+uoLprior+uPprior+Yprior+p0prior+FpPrior+oprior)
  return(priorSum)
}


# set bounds on initial parameter guesses for if using random start values
initBounds <- data.frame(rbind( ## for initial conditions
  c(0.03,0.04), ## uoE
  c(0.03,0.04), ## uoL
  c(0.2,0.3), ## uP
  c(8,15), ## Y
  c(25,25),##n
  c(0.4,0.6),## p0
  c(0.01,0.9),##o
  c(1,5),#fracPop
  c(2,6),
  c(2,6),
  c(2,6),
  c(2,6)))##sf


colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('uoE','uoL','uP','Y','n','p0','logo','logFp','sf1','sf2','sf3','sf4')
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
                else  propVal <- proposal[on + 2] + rnorm(1, mean = 0, sd = sdProps[on + 2])
                proposal[on + 2]<- if(propVal<0) 0 else propVal
                on <<- (on+2) %% nfitted
                return(proposal)
              }))
}

##proposal sd tuning function
tuner <- function(curSd, acptR,curAcptR,minSddProps){
  print(paste0("aratio = ",curAcptR))
  
  if(curAcptR ==1) curAcptR <- 0.99
  if(curAcptR == 0) curAcptR <- 0.01
  curSd = (curSd*qnorm(acptR/2))/qnorm(curAcptR/2)
  #curSd[curSd > 1] <- 1
  #curSd[curSd < sdProps] <- 0.01
  curSd[c(which(curSd<minSddProps))]<-minSddProps[c(which(curSd<minSddProps))]
  
  return(curSd)
}


## function for adaptive blocked proposals based on var-covar matrix
multiv.proposer <- function(covar,blockLS = list(rownames(covar))) {
  return(list(type = 'block',
              fxn = function(current, sdTune) {
                proposal <- current + (rmnorm(1, mean = 0, varcov = covar)*sdTune)
                propsosal <- as.vector(proposal)
                proposal[proposal<=0]<-1e-3
                names(proposal) <- names(current)
                proposal
              }))
}

##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    Metropolis-Hastings MCMC Sampler                                                            #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################
library(coda)
mcmcSampler <- function(initParams, ## initial parameter guess
                        randInit = T, ## if T then randomly sample initial parameters instead of above value
                        obsDat = myDat, ## data
                        proposer = multiv.proposer(covar),## proposal distribution
                        sdProps=c(0.1,0.1,0.1,1,0,0.1,0.1,0.1,0.1,2,2,2,2),##starting proposal dists
                        minSddProps=c(0.01,0.01,0.01,0.1,0,0.01,0.01,0.01,0.01,0.3,0.3,0.3,0.3),#min values for tuner for each param
                        niter = 100, ## MCMC iterations
                        particles =100,##number of particles for particle filter
                        nburn = 0, ## iterations to automatically burn
                        monitoring=0, ## if >2 browses, if >1 prints progress
                        adaptiveMCMC = F, ## adapt proposal distribution?
                        startAdapt = 150, ## start adapting at what iteration?
                        adptBurn = 200, ## ignore first x number of iterations for adapting posterior
                        acceptanceRate=0.9,##acceptance rate for adaptive mcmc
                        tell = 100) { ## how often to print progress
  
  aratio<-0#starting acceptance ratio
  sdp<-sdProps
  
  if(monitoring>2) browser()
  if(randInit) initParams <- initRand(initParams)
  currentParams <- initParams
  nfitted <- length(currentParams) ## number fitted parameters
  iter <- 2 ## mcmc iteration (started at 1 so we're already on 2
  accept <- 0 ## initialize proportion of iterations accepted
  acceptR<- 0 #number of accepts for aratio
  ## Calculate log(likelihood X prior) for first value
  curVal <- llfnc(currentParams, particles=particles)
  print(curVal)
  ## Initialize matrix to store MCMC chain
  out <- matrix(NA, nr = niter, nc=length(currentParams)+1)
  out[1,] <- c(currentParams, ll = -curVal) ## add first value
  colnames(out) <- c(names(currentParams), 'll') ## name columns
  ## Store original covariance matrix
  if(proposer$type=='block') originalCovar <- get('covar', envir = environment(proposer$fxn)) 
  while(iter <= niter) {
    if ((monitoring > 1) || (monitoring && (iter%%tell == 0))) print(paste("on iteration",iter,"of", niter + 1))
    
    ##var covar matrix update - currently every 50 iterations
    if(adaptiveMCMC & proposer$type=='block' & iter > startAdapt & iter %% 50 == 0) { ##modulur division of 50, update covar every 50 iterations
      adptBurn <- min((startAdapt-50), adptBurn)
      adaptedCovar <- 2.38^2 / nfitted-1 * cov(log(out[adptBurn:(iter-1),c(1:4,6:13)]))
      adaptedCovar <- adaptedCovar*.95 + originalCovar*.05 ## 95% adapted & 5% original
      rownames(adaptedCovar) <- colnames(adaptedCovar) <- names(currentParams)
      assign('covar', adaptedCovar, envir = environment(proposer$fxn))
    }
    
    if(proposer$type=='block'){sdp<-tuner(sdp,acceptanceRate,aratio,minSddProps)
    print(sdp)
    proposal <- proposer$fxn(currentParams,sdTune=sdp)}
    else proposal <- proposer$fxn(currentParams)
    print(proposal)
    propVal <- llfnc(proposal, particles=particles)
    lmh <- propVal - curVal ## likelihood ratio = log likelihood difference
    if (is.na(lmh)) { ## if NA, do not accept
    } else { ## if it's not NA then do acception/rejection algorithm
      if (monitoring > 1) print( c(lmh=lmh, propVal=propVal, curVal=curVal))
      if ( (lmh >= 0) | (runif(1,0,1) <= exp(lmh)) ) {
        currentParams <- proposal
        acceptR<-acceptR + 1
        if (iter>nburn) accept <- accept + 1 ## only track acceptance after burn-in
        curVal <- propVal
        
      }
    }
  
    out[iter, ] <- c(currentParams, ll=curVal)
    iter <- iter+1
    aratio <- acceptR/(iter)
  }
  colnames(out) <- c(names(currentParams), 'll')
  results <- as.mcmc(out[1:nrow(out)>(nburn+1),])
  return(list(initParams = initParams
              , aratio = aratio
              , results = results
  ))
}

#Error in out[adptBurn:(iter - 1), 1:nfitted] : 
# only 0's may be mixed with negative subscripts 
