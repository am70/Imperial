
##for use when running parralel within particlefilterMCMC
runFilt<-function(runs, particles,thetaList,parameter){
  
  th<-thetaList[runs]
  parms <- theta()
  parms[parameter] <- th 
  
  x<-particleFilter(larvalModP,
                    theta=parms,
                    init.state,
                    data = garkiObs,
                    n.particles = particles)
  return(x)
}

#for use when running parallel on single node with lots of cores
runFiltAll<-function(runs, particles,theta){

  x<-particleFilter(larvalModP,
                    theta=theta,
                    init.state,
                    data = garkiObs,
                    n.particles = particles)
  return(x)
}



parRunS<-function(x,y,z,p){parallel::parLapply(NULL,x,function(x)runFilt(x,y,z,p))}
parRunSAll<-function(runs,particles,theta){parallel::parLapply(NULL,runs,function(runs)runFiltAll(runs,particles,theta))}


modStep<-function(pr, initialState, times){
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],n=pr[4],sf=pr[5],Y=pr[6])
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(times[1]*10, times[2]*10, length.out=times[2]-times[1])))
  res1<-(simDat[simDat$step/10 == times[1],])
  res2<-(simDat[simDat$step/10 == times[2],])
  return(as.data.frame(rbind(res1,res2)))
}




# set up data likelihood
dataLik <- function(simPoint, obsDat, theta, log = TRUE)
{
  ll = sum(dzipois(obsDat[2]*25, simPoint[4]*25, log = TRUE))
  return(exp(ll/25))
}



particleFilter <-
  function(fitmodel,
           theta,
           init.state,
           data,
           n.particles) {
    
    margLogLike <- 0
    
    state.particles  <- rep(list(init.state), n.particles)
    
    weight.particles <- rep(1 / n.particles, length = n.particles)
    
    current.time <- 0
    
    for (i in seq_len(nrow(data))) {
      data.point <- unlist(data[i,])
      next.time <- data.point["time"]
      
      # Resample particles according to their weights.
      weight.particles[is.na(weight.particles)] <- 0
      index.resampled <- sample(
        x = n.particles,
        size = n.particles,
        replace = TRUE,
        prob = 0.001+as.numeric(weight.particles)
      )
      
      state.particles <- state.particles[index.resampled]
      
      weightP<-function(weightInput){
        
        dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
        colnames(dataInput)<-c("E","L","P","M","startTime","nextTime")
        state.particles<-as.numeric(dataInput[,c(1:4)])
        current.time<-as.numeric(dataInput[,c(5)])
        next.time<-as.numeric(dataInput[,c(6)])
        current.state.particle <- state.particles
        # Propagate the particle from current observation time to the next one 
        traj <-  modStep(pr = theta,
                         initialState = current.state.particle,
                         times =c(current.time,next.time))
        # Extract state of the model at next observation time
        model.point <- unlist(traj[2, c("E","L","P","M")])
        
        # Weight the particle with the likelihood of the observed
        weightX <-dataLik(obsDat = data.point,
                          simPoint = model.point,
                          theta = theta)
        
        res<-(c(as.numeric(weightX),as.numeric(model.point)))
        return(res)
      }
      
      
      dat<-data.frame(t(sapply(state.particles,c)))
      colnames(dat)<-c("E","L","P","M")
      dat<- paste(dat$E,dat$L,dat$P, dat$M,current.time,next.time,sep=",")
      
      
      
      partWeights<-lapply(dat,weightP)
      
      state.particles <- lapply(partWeights, '[', c(2:5))
      weight.particles<-lapply(partWeights, '[', 1)
      
      # Increment time
      current.time <- next.time
      
      ## Increment the marginal log-likelihood
      # Add the log of the mean of the particles weights
      margLogLike <- margLogLike + log(mean(as.numeric(weight.particles)))
    }
    
    ## Return marginal log-likelihood
    return(margLogLike)
    
  }


theta <-
  function(
    uoE = 0.1343961,uoL = 0.01005042,uP = 0.06521256,n =  10,sf = 2.972,y = 12.62243)
    return(c(uoE=uoE,uoL=uoL,uP=uP,n=n,sf=sf,y=y))


init.state <-
  c(E = 177,L = 8,P = 1,M = 7)



particleFilterMCMC <-
  function(fitmodel,
           theta,
           init.state,
           data,
           n.particles) {
    
    margLogLike <- 0
    
    state.particles  <- rep(list(init.state), n.particles)
    
    weight.particles <- rep(1 / n.particles, length = n.particles)
    
    current.time <- 0
    
    for (i in seq_len(nrow(data))) {
      data.point <- unlist(data[i,])
      next.time <- data.point["time"]
      
      # Resample particles according to their weights.
      weight.particles[is.na(weight.particles)] <- 0
      index.resampled <- sample(
        x = n.particles,
        size = n.particles,
        replace = TRUE,
        prob = 0.001+as.numeric(weight.particles)
      )
      
      state.particles <- state.particles[index.resampled]
      
      
      
      dat<-data.frame(t(sapply(state.particles,c)))
      colnames(dat)<-c("E","L","P","M")
      
      dat<- paste(dat$E,dat$L,dat$P, dat$M,current.time,next.time,data.point[1],data.point[2],sep=",")
clVal<-length(dat)/10
      
      partWeights1<-obj$enqueue(pWeightApply(dat[1:(clVal*1)]),name="particle filter MCMC1")
      partWeights2<-obj$enqueue(pWeightApply(dat[(1+clVal*1):(clVal*2)]),name="particle filter MCMC2")
      partWeights3<-obj$enqueue(pWeightApply(dat[(1+clVal*2):(clVal*3)]),name="particle filter MCMC3")
      partWeights4<-obj$enqueue(pWeightApply(dat[(1+clVal*3):(clVal*4)]),name="particle filter MCMC4")
      partWeights5<-obj$enqueue(pWeightApply(dat[(1+clVal*4):(clVal*5)]),name="particle filter MCMC5")
      partWeights6<-obj$enqueue(pWeightApply(dat[(1+clVal*5):(clVal*6)]),name="particle filter MCMC6")
      partWeights7<-obj$enqueue(pWeightApply(dat[(1+clVal*6):(clVal*7)]),name="particle filter MCMC7")
      partWeights8<-obj$enqueue(pWeightApply(dat[(1+clVal*7):(clVal*8)]),name="particle filter MCMC8")
      partWeights9<-obj$enqueue(pWeightApply(dat[(1+clVal*8):(clVal*9)]),name="particle filter MCMC9")
      partWeights10<-obj$enqueue(pWeightApply(dat[(1+clVal*9):(clVal*10)]),name="particle filter MCMC10")
      
      p1<-partWeights1$wait(1e6)
      p2<-partWeights2$wait(1e6)
      p3<-partWeights3$wait(1e6)
      p4<-partWeights4$wait(1e6)
      p5<-partWeights5$wait(1e6)
      p6<-partWeights6$wait(1e6)
      p7<-partWeights7$wait(1e6)
      p8<-partWeights8$wait(1e6)
      p9<-partWeights9$wait(1e6)
      p10<-partWeights10$wait(1e6)
      
      
      partWeights<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)

      state.particles <- lapply(partWeights, '[', c(2:5))
      weight.particles<-lapply(partWeights, '[', 1)
      
      # Increment time
      current.time <- next.time
      
      ## Increment the marginal log-likelihood
      # Add the log of the mean of the particles weights
      margLogLike <- margLogLike + log(mean(as.numeric(weight.particles)))
    }
    
    ## Return marginal log-likelihood
    return(margLogLike)
    
  }

weightPmcmc<-function(weightInput){
  
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  colnames(dataInput)<-c("E","L","P","M","startTime","nextTime","dpTime","dpM")
  data.point<-head(dataInput[,c("dpTime","dpM")],1)
  state.particles<-as.numeric(dataInput[,c(1:4)])
  current.time<-as.numeric(dataInput[,c(5)])
  next.time<-as.numeric(dataInput[,c(6)])
  current.state.particle <- state.particles
  # Propagate the particle from current observation time to the next one 
  traj <-  modStep(pr = theta(),
                   initialState = current.state.particle,
                   times =c(current.time,next.time))
  # Extract state of the model at next observation time
  model.point <- unlist(traj[2, c("E","L","P","M")])
  
  # Weight the particle with the likelihood of the observed
  weightX <-dataLik(obsDat = data.point,
                    simPoint = model.point,
                    theta = theta())
  
  res<-(c(as.numeric(weightX),as.numeric(model.point)))
  return(res)
}

pWeightApply<-function(x){parallel::parLapply(NULL,x,weightPmcmc)}

