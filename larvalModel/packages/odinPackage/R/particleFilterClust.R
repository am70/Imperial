
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
  print(pr[1])
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],sf=pr[6])
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(times[1]*10, times[2]*10, length.out=times[2]-times[1])))
  res1<-(simDat[simDat$step/10 == times[1],])
  res2<-(simDat[simDat$step/10 == times[2],])
  return(as.data.frame(rbind(res1,res2)))
}




# set up data likelihood function 
dataLik <- function(simPoint, obsDat, theta, log = TRUE)
{
  ll = dzipois(obsDat[2]*25, simPoint[4]*25, pstr0=0, log = TRUE)
  ell=exp(ifelse(ll<-700,-700,ll))
  return (ell)
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
      nextTime <- dataPoint["time"]
      
      # Resample particles according to their weights.
      weightParticles=as.numeric(weightParticles)
      weightParticles[is.na(weightParticles)] <- 0
      swP=sum(weightParticles)
      norm.weightedParticles=weightParticles/swP
      indexResampled <- sample(
        x = nParticles,
        size = nParticles,
        replace = TRUE,
        prob = norm.weightedParticles # incase of NA's add a baseline probability
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
                          simPoint = modelPoint,
                          theta = thetaX)
        
        res<-(c(as.numeric(weightX),as.numeric(modelPoint)))
        return(res)
      }
      
      
      dat<-data.frame(t(sapply(stateParticles,c)))
      colnames(dat)<-c("E","L","P","M")
      
      partWeights<-parallel::parLapply(cl,paste(dat$E,dat$L,dat$P, dat$M,currentTime,nextTime,sep=","),weightP)     
      
      stateParticles <- lapply(partWeights, '[', c(2:5))
      weightParticles<-lapply(partWeights, '[', 1)
      
      # Increment time
      currentTime <- nextTime
      
      ## Increment the marginal log-likelihood
      # Add the log of the mean of the particles weights
      margLogLike <- margLogLike + log(mean(as.numeric(weightParticles)))
    }
    
    ## Return marginal log-likelihood
    return(margLogLike)
    
  }


theta <-
  function(
    uoE = 0.1736114,uoL = 0.4792658,uP = 0.2163627,y = 0.1301954,n =  2.932363,sf = 13.74987)
    return(c(uoE=uoE,uoL=uoL,uP=uP,y=y,n=n,sf=sf))

thetaGarki <-
  function(
    uoE = 0.1392127,uoL = 0.03115276,uP = 0.7523195,y = 3.150944,n =  3.911705,sf = 21.89244)
    return(c(uoE=uoE,uoL=uoL,uP=uP,y=y,n=n,sf=sf))


init.state <-
  c(E = 177,L = 8,P = 1,M = 7)


