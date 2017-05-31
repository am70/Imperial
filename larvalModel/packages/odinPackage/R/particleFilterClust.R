
##for use when running parralel within particlefilterMCMC
runFilt<-function(runs, particles,thetaList,parameter){
  
  th<-thetaList[runs]
  parms <- theta()
  parms[parameter] <- th 
  
  x<-particleFilter(larvalModP,
                    thetaX=parms,
                    init.state,
                    data = garkiObs,
                    nParticles = particles)
  return(x)
}

#for use when running parallel on single node with lots of cores - need to figure out rrq function?
runFiltAll<-function(runs, particles,theta){

  x<-particleFilter(larvalModP,
                    thetaX=theta,
                    init.state,
                    data = garkiObs,
                    nParticles = particles)
  return(x)
}



parRunS<-function(x,y,z,p){parallel::parLapply(NULL,x,function(x)runFilt(x,y,z,p))}
parRunSAll<-function(runs,particles,theta){parallel::parLapply(NULL,runs,function(runs)runFiltAll(runs,particles,theta))}
parRunWorker<-function(x,y,z,p){lapply(x,function(x)runFilt(x,y,z,p))}#for use with single node workers




##model step function - runs model in steps taking the liklihood of the subsequent observed data point
modStep<-function(pr, initialState, times){
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],sf=pr[6])
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(times[1]*10, times[2]*10, length.out=times[2]-times[1])))
  res1<-(simDat[simDat$step/10 == times[1],])
  res2<-(simDat[simDat$step/10 == times[2],])
  return(as.data.frame(rbind(res1,res2)))
}




# set up data likelihood function 
dataLik <- function(simPoint, obsDat,phi)
{
  ll = dzipois(obsDat[2], simPoint[4], pstr0=0, log = T)
  return (exp(ll))
}

#return log likelihood subtract off the maximum value at each point, exp then add back in
#write down log likelihood of poisson, rather than rely on function calls

particleFilter <-

  function(fitmodel,
           thetaX,
           init.state,
           data,
           nParticles) {
    
  
    margLogLike <- 0
    wpM<-data.frame(1:nParticles)
    colnames(wpM)<-c("particleNumber")
    
    
    #create initial particle states
    stateParticles  <- rep(list(init.state), nParticles)
    
    #weight particles evenly
    weightParticles <- rep(1 / nParticles, length = nParticles)
    
    #initiate start time
    currentTime <- 0
    
    #loop to run the model in steps using starting values drawn at random 
    #but weighted in favour of the most likely values from the previous step. 
   
     for (i in seq_len(nrow(data))) {
      dataPoint <- unlist(data[i,])
      dataPointAnc <- if(i>1)unlist(data[i-1,]) else unlist(data[i,])
      
      nextTime <- dataPoint["time"]
      # Resample particles according to their weights.
      weightParticles=as.numeric(weightParticles)

      weightParticles<-log(weightParticles)-max(log(weightParticles))
      weightParticles<-(exp(weightParticles))
     
      swP=sum(weightParticles)
      norm.weightedParticles=weightParticles/swP

      indexResampled <- sample(
        x = nParticles,
        size = nParticles,
        replace = TRUE,
        prob = norm.weightedParticles 
      )
      stateParticles <- stateParticles[indexResampled]


      #main function for weighting particles
      weightP<-function(weightInput){
        
        dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
        stateParticles<-as.numeric(dataInput[,c(1:4)])
        currentTime<-as.numeric(dataInput[,c(5)])
        nextTime<-as.numeric(dataInput[,c(6)])
        current.state.particle <- stateParticles
        # Propagate the particle from current observation time to the next one 
        traj <-  modStep(pr = thetaX,
                         initialState = current.state.particle,
                         times =c(currentTime,nextTime))
     
        # Extract state of the model at next observation time
        modelPoint <- unlist(traj[2, c("E","L","P","M")])

        # Weight the particle with the likelihood of the observed
        weightX <-dataLik(obsDat = dataPoint,
                          simPoint = modelPoint
                          )
        
        
        res<-(c(as.numeric(weightX),as.numeric(modelPoint)))
        return(res)
      }
      
      
      dat<-data.frame(t(sapply(stateParticles,c)))
      colnames(dat)<-c("E","L","P","M")
      
      partWeights<-lapply(paste(dat$E,dat$L,dat$P, dat$M,currentTime,nextTime,sep=","),weightP)     
      
      stateParticles <- lapply(partWeights, '[', c(2:5))
      weightParticles<-lapply(partWeights, '[', 1)
      print((stateParticles))
      print(log(unlist(weightParticles)))
      wpMtemp<-as.data.frame(as.vector(log(unlist(weightParticles))))
      colnames(wpMtemp) <- paste("Obs", i, sep = "_")
      wpM<-cbind(wpM,wpMtemp)

      # Increment time
      currentTime <- nextTime
      
      ## Increment the marginal log-likelihood
      # Add the log of the mean of the particles weights
      margLogLike <- margLogLike + log(mean(as.numeric(weightParticles)))
    }
    
    wpMRes<-sample(
      x = nParticles,
      size = 1,
      replace = TRUE,
      prob = exp(wpM$Obs_12) 
    )
    
   
    resLike<-( wpM[wpM$particleNumber == wpMRes,])
    
    wpM$sumLikes<-rowSums(wpM[,-1])
    ## Return marginal log-likelihood
    #return(sum(resLike[,-1]))
    return(wpM)
    
    #NEED TO KEEP RUNNING LIKLIHOOD TOTAL FOR EACH PARTICLE THEN SAMPLE A FINAL SINGLE PARTICLE BASED ON FINAL WEIGHTS
  }


theta <-
  function(
    uoE = 0.1736114,uoL = 0.4792658,uP = 0.2163627,y = 0.1301954,n =  2.932363,sf = 13.74987)
    return(c(uoE=uoE,uoL=uoL,uP=uP,y=y,n=n,sf=sf))

thetaGarki <-
  function(
    uoE = 0.2732377,uoL = 0.4218158,uP = 0.5580947,y = 1.168431,n =  4.397768,sf = 24.79025)
    return(c(uoE=uoE,uoL=uoL,uP=uP,y=y,n=n,sf=sf))


#test params for simDat2
# uoE        uoL        uP        Y  n sf  
# 0.5496511 0.03253728 0.5188584 13 15 17.63029 

init.state <-
  c(E = 177,L = 8,P = 1,M = 7)



####second method#########################

t0=0


# now define a sampler for the prior on the initial state
simx0 <- function(N,t0)
{
  mat=cbind(rpois(N,177),rpois(N,8),rpois(N,1),rpois(N,7))
  colnames(mat)=c("E","L","P","M")
  mat
}

##model step function - runs model in steps taking the liklihood of the subsequent observed data point
modStep2<-function(weightInput){
  
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  initialState<-as.numeric(dataInput[,c(1:4)])
  currentTime<-as.numeric(dataInput[,c(5)])
  nextTime<-as.numeric(dataInput[,c(6)])
  pr<-as.numeric(dataInput[,c(7:12)])
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],sf=pr[6])
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(currentTime, nextTime, length.out=nextTime-currentTime)))
  res2<-(simDat[simDat$step == nextTime,])
  return(as.data.frame(rbind(res2[,c(8:11)])))
}


dataLik2 <- function(input)
{
  dataInput<-read.table(text = input, sep = ",", colClasses = "numeric")
  ll = dzipois(dataInput[2], dataInput[1], pstr0=0, log = T)
  return (exp(ll))
}

#n=particle number
#simx0 =

pfMLLik <- function (n, simx0, t0, stepFun, dataLik, data,pr) 
{
  times = c(as.numeric(data$time))
  xmat = simx0(n, t0) #initial state
  ll = 0
  for (i in 1:length(times[-length(times)])) {
    wp<-paste(xmat[,1],xmat[,2],xmat[,3], xmat[,4],times[i],
              times[i + 1],pr[1],pr[2],pr[3],pr[4],pr[5],pr[6],sep=",")
    
    xmatTemp = parLapply(cl,wp,stepFun)
    xmat<-data.frame(t(sapply(xmatTemp, `[`)))
    
    likeDat<-paste(xmat$M,times[i],sep=",")
    w = lapply(likeDat,dataLik)
    w<-as.vector(unlist(w))
    
    swP=sum(w)
    w=w/swP
    
    w<-log(w)-max(log(w))
    w<-(exp(w))
    
    ll = ll + log(mean(w))
    w[is.na(w)] <- 1e-200
    rows = sample(1:n, n, replace = TRUE, prob = w)
    xmat = xmat[rows, ]
  }
  
  ll
  
}


