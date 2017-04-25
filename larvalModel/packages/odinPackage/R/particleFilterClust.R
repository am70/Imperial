

runFilt<-function(runs, particles){
  
  x<-particleFilter(larvalModP,
                    theta(),
                    init.state,
                    data = garkiObs,
                    n.particles = particles)
  return(x)
}



parRunS<-function(x,y){parallel::parLapply(NULL,x,function(x)runFilt(x,y))}


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
  ll = sum(dzipois(obsDat*25, simPoint[4]*25, log = TRUE))
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

