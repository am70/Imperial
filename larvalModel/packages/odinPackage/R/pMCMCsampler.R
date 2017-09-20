

##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    pMCMC functions                                                                             #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

##function that sums log-likelihood & log-prior inside MCMC sampler - runs for both local (sf) parameters on each targeted village and global parameters
# @param fitParams new model parameters
# @param particles number of particles

# @param fxdParams fixed parameters
# @param cluster T/F if true, using dide cluster, if false running locally
# @return log likelihood for global and local parameters
llfnc <- function(fitParams, particles, fxdParams, cluster) {
  pX <- NULL
  for (i in 1:4) {
    globalParms <-
      fitParams[c(1:8)]
    globalParms[8] <-
      fitParams[c(i + 7)]#i+8 to fit the scaling factor specific for each village
    garkDat <-
      garkiObsX[, c(1, i + 1)]#i+1 as first column is "time"
    garkDat <- na.omit(garkDat)
    rFc <-
      as.numeric(substring(names(garkDat)[2], 1, 1))#which rainfall cluster to use
    #run particle filter for village
    p1 <- 
      pFilt(particles,iState,modStep3,dataLikFunc,garkDat,pr = globalParms,rFclust = rFc,fxdParams = fxdParams,resM = F, cluster = cluster) + lprior(globalParms)
    pX <- 
      rbind(pX, p1)
  }
  return(sum(pX))
}




## Log-Prior function - specific to this model
# @parms current parameters
# @return sum of log priors
lprior <- function(parms) {
  
 priorSum<-(dnorm(parms[1], mean = 0.035, sd = 0.06, log = T)
  +dunif(parms[1], min = 0.001, max = 0.99, log = T)
  +dnorm(parms[2], mean = 0.035, sd = 0.06, log = T)
  +dunif(parms[2], min = 0.001, max = 0.99, log = T)
  +dnorm(parms[3], mean = 0.25, sd = 0.11, log = T)
  +dunif(parms[3], min = 0.001, max = 0.99, log = T)
  +dnorm(parms[4], mean = 13.06, sd = 4.53, log = T)
  +dunif(parms[5], min = 1, max = 1e+5, log = T)
  +dunif(parms[6], min = 1e-06, max = 1, log = T)
  +dunif(parms[5], min = 0.01, max = 1000, log = T)
  +dunif(parms[7], min = 1, max = 10, log = T)
  +dunif(parms[8], min = 1, max = 1e+20, log = T)
  )
return(as.vector(priorSum))
}


# set bounds on initial parameter guesses if using random start values
initBounds <- data.frame(rbind(
  ## for initial conditions
  c(0.03, 0.04), ## uoE
  c(0.03, 0.04), ## uoL
  c(0.2, 0.3), ## uP
  c(8, 15), ## Y
  c(0.4, 0.6), ## p0
  c(0.01, 0.9), ##o
  c(1, 5),#fracPop
  c(2, 6),##sf1:4
  c(2, 6),
  c(2, 6),
  c(2, 6)
))
colnames(initBounds) <- c('lower', 'upper')
rownames(initBounds) <-
  c('uoE','uoL','uP','Y','p0','logo','logFp','sf1','sf2','sf3','sf4')
class(initBounds[, 2]) <- class(initBounds[, 1]) <- 'numeric'
initBounds


##  randomly select a value that is uniformly distributed between these bounds
# @param fitParams list of model parameters
# @return random parameters from within set boundaries 
initRand <- function(fitParams) {
  fitParams <- logParms(fitParams)
  tempnm <- names(fitParams)
  for (nm in tempnm)
    fitParams[nm] <-
    runif(1, min = initBounds[rownames(initBounds) == nm, 'lower'],
          max =  initBounds[row.names(initBounds) ==
                              nm, 'upper'])
  return(fitParams)
}


## Sequential proposal function
# @param current current parameters
# @param sdTune s.d. tuning from tuner function
# @param prmNum which parameter to update
# @param sdProps 
# @param tune
# @return
sequential.proposer <- function(current, prmNum, sdProps) {
  proposal <- current
  propVal <-
    proposal[prmNum] + (rnorm(1, mean = 0, sd = sdProps[prmNum]))
  proposal[prmNum] <- propVal
  return(proposal)
}



## function for adaptive blocked proposals based on var-covar matrix
# @param current current parameters
# @param sdTune tuned s.d.
# @return proposal parameters
multiv.proposer <-  function(current, sdProps) {
  proposal <- 
    current + (rmnorm(1, mean = 0, varcov = covar) * sdTune)
  propsosal <- 
    as.vector(proposal)
  names(proposal) <- 
    names(current)
  proposal
}


##proposal sd tuning function
# @param current s.d.
# @param target acceptance ratio
# @param current acceptance ratio
# @param maximum proposal s.d.
# @return tuned s.d.
tuner <- function(curSd, acptR, curAcptR, maxSddProps) {
  if (curAcptR == 1)
    curAcptR <- 0.99
  if (curAcptR == 0)
    curAcptR <- 0.01
  curSd = (curSd * qnorm(acptR / 2)) / qnorm(curAcptR / 2)
  curSd[c(which(curSd < maxSddProps))] <-
    maxSddProps[c(which(curSd < maxSddProps))]
  return(curSd)
}


##proposal sd tuning function for sequential  adaptive tuning
# @param current s.d.hh
# @param target acceptance ratio
# @param current acceptance ratio
# @param maximum proposal s.d.
# @param i which parameter to work on
# @return tuned s.d.
tunerSeq <- function(curSd, acptR, curAcptR, maxSddProps, i) {
  if (curAcptR[i] == 1)
    curAcptR[i] <- 0.99
  if (curAcptR[i] == 0)
    curAcptR[i] <- 0.01
  curSd[i] = (curSd[i] * qnorm(acptR[i] / 2)) / qnorm(curAcptR[i] / 2)
  curSd[i][curSd[i] > maxSddProps[i]] <- maxSddProps[i]
  return(curSd[i])
}

##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    Metropolis-Hastings MCMC Sampler                                                            #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################
library(coda)
# @param initial parameter guess
# @param if T then randomly sample initial parameters instead of initParams value
# @param fixedParam specific fixed parameters in model, currently just n - needs updating to include more than one value for future flexibility
# @param proposer proposal function, multivariate block (adaptiveMCMC must = T) or sequential can have adaptive or not adaptive tuning
# @param sdProps standard deviation for proposal distributions - in adaptive this is the starting sd
# @param maxSddProps maximum values for sd proposals if using adaptive MCMC
# @param niter number of iterations to run the MCMC for
# @param particles number of particles for particle filter
# @param nburn number of mcmc iterations to burn
# @param monitoring 0 = no monitoring, > 0 prints more progress information
# @param adaptiveMCMC T/F if true uses tuning for s.d. of parmater proposal distributions based on acceptance ratios
# @param proposerType "seq" or "block", seq for sequential proposing, block for blocked proposing, blocks based on var-covar matrix - block proposing only available when using adaptiveMCMC
# @param startAdapt starting iteration for adapting
# @param adaptBurn burn n number of iterations for defining var-covar matrix in block proposing
# @param acceptanceRate acceptance rates, for adaptive sequential use string, one rate for each param, for block use single value
# @param tell print monitoring information every x number of iterations
# @param cluster T/F if true, using dide cluster, if false running locally
mcmcSampler <- function(initParams,
                        randInit = T,
                        fixedParam = 40,
                        proposer = sequential.proposer,
                        sdProps = c(0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 0.1, 2),
                        maxSddProps = c(0.01, 0.01, 0.01, 0.1, 0.01, 0.01, 0.01, 0.01, 0.3),
                        niter = 100,
                        particles = 100,
                        nburn = 0,
                        monitoring = 0,
                        adaptiveMCMC = F,
                        proposerType = 'seq',
                        startAdapt = 150,
                        adptBurn = 200,
                        acceptanceRate = 0.9,
                        tell = 100,
                        cluster = F) {
  aratio <- rep(0, length(sdProps))#starting acceptance ratio
  acceptRseq <- rep(0, length(sdProps))
  iterR <- rep(0, length(sdProps))
  sdp <- sdProps
  prmNum <- 1 # which paramter is currently being fitted
  if (randInit) initParams <- initRand(initParams)
  currentParams <- initParams
  nfitted <- length(currentParams) ## number fitted parameters
  iter <- 2 ## mcmc iteration started at 1 so already on 2
  accept <- 0 ## initialize proportion of iterations accepted
  acceptR <- 0 #number of accepts for aratio
  ## Calculate log likelihood for first value
  curVal <-
    llfnc(currentParams,
          particles = particles,
          fxdParams = fixedParam,
          cluster)
  ## Initialize matrix to store MCMC chain
  out <- matrix(NA, nr = niter, nc = length(currentParams) + 1)
  out[1,] <- c(currentParams, ll = -curVal) ## add first value
  colnames(out) <- c(names(currentParams), 'll') ## name columns
  if (adaptiveMCMC == T & proposerType == 'block')
    originalCovar <-
    get('covar', envir = environment(proposer))## Store original covariance matrix
  while (iter <= niter) {
 
    
    ##var covar matrix update - currently every 50 iterations
    if (adaptiveMCMC == T &
        proposerType == 'block' & iter > startAdapt & iter %% 50 == 0) {
      ##modulur division of 50, update covar every 50 iterations
      adptBurn <- min((startAdapt - 50), adptBurn)
      adaptedCovar <-
        (2.38 ^ 2 / nfitted) * cov(log(out[adptBurn:(iter - 1), 1:nfitted]))
      adaptedCovar <-
        adaptedCovar * .95 + originalCovar * .05 ## 95% adapted & 5% original
      rownames(adaptedCovar) <-
        colnames(adaptedCovar) <- names(currentParams)
      assign('covar', adaptedCovar, envir = environment(proposer))
    }
    
    if (adaptiveMCMC == T & proposerType == 'block') {
      sdp <- tuner(sdp, acceptanceRate, aratio, maxSddProps)
      proposal <- 
        proposer(currentParams, sdProps = sdp)
    }
    
    if (adaptiveMCMC == T & proposerType == 'seq') {
      if (iter >= startAdapt)
        sdp[prmNum] <-
          tunerSeq(sdp, acceptanceRate, aratio, maxSddProps, prmNum)
      proposal <- 
        proposer(currentParams, prmNum, sdp)
    
    }
    
    if (adaptiveMCMC == F & proposerType == 'seq') {
      proposal <- 
        proposer(currentParams, prmNum, sdp)
    }
    
    propVal <-
      llfnc(proposal,
            particles = particles,
            fxdParams = fixedParam,
            cluster)
    lmh <-
      propVal - curVal ## likelihood ratio = log likelihood difference
    
    if ((monitoring > 1 && iter %% tell == 0)){
      print(paste("on iteration", iter, "of", niter + 1))
      print(curVal)
      print(sdp)
      print(iterR)
      print(paste0("a ratio ", aratio))
      print(proposal)
    print(c(
      lmh = lmh,
      propVal = propVal,
      curVal = curVal
    ))}
    
    if (is.na(lmh)) {
      ## if NA, do not accept
    } else {
      ## if it's not NA then do acception/rejection algorithm
      if ((lmh >= 0) | (runif(1, 0, 1) <= exp(lmh))) {
        currentParams <- proposal
        acceptRseq[prmNum] <- acceptRseq[prmNum] + 1
        if (iter > nburn)
          accept <- accept + 1 ## track acceptance after burn-in
        curVal <- propVal
      }
    }
    out[iter, ] <- c(currentParams, ll = curVal)
    iter <- iter + 1
    iterR[prmNum] <- iterR[prmNum] + 1

    aratio[prmNum] <-
      acceptRseq[prmNum] / (iterR[prmNum])#acceptrance ratio change for specific parameter number
    prmNum <- prmNum + 1#progress parameter number
    if (prmNum > length(sdProps))
      prmNum <-1#if parameter number reaches end of parameters, switch back to start
    
  }
  colnames(out) <- c(names(currentParams), 'll')
  results <- as.mcmc(out[1:nrow(out) > (nburn + 1),])
  return(list(
    initParams = initParams
    ,
    aratio = aratio
    ,
    results = results
  ))
}
