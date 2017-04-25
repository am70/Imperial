#########################################################################################################################################
#                                                                                                                                       #
#                                                     particle MCMC                                                                     #
#                                                                                                                                       #
#########################################################################################################################################


load_all()
document()
devtools::install(args = "--no-multiarch")

cl <- makeSOCKcluster(c(7),outfile="")
clusterEvalQ(cl, library("odinPackage"))
clusterEvalQ(cl, library("lhs"))
clusterEvalQ(cl, library("VGAM"))
clusterEvalQ(cl, library("deSolve"))
res<-NULL
clusterExport(cl, c("res","rFx","delta","garkiObs"), envir=environment())




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
      
      partWeights<-(parLapply(cl,dat,weightP))
      
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



#theta <-
 # c(uoE = 0.01,uoL = 0.01,uP = 0.25,n = 5,sf = 2,y = 1)

theta <-
  function(
  uoE = 0.1343961,uoL = 0.01005042,uP = 0.06521256,n =  10,sf = 2.972,y = 12.62243)
  return(c(uoE=uoE,uoL=uoL,uP=uP,n=n,sf=sf,y=y))


init.state <-
  c(E = 177,L = 8,P = 1,M = 7)


clusterExport(cl, c("modStep","theta","init.state","dataLik"), envir=environment())

set.seed(10)
# run the particle filter
for(i in c(1:10)){

  x<-particleFilter(larvalModP,
                    theta(),
                    init.state,
                    data = garkiObs,
                    n.particles = 500)
  print(x)
}





####################################################################################
####################################################################################
####################################################################################
####################################################################################





test.n.particles <- c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)

n.replicates <- 100

sample.log.like <- vector("numeric", length = n.replicates)
res <- data.frame()

for(n.particles in test.n.particles){
  
  start.time  <- Sys.time()
  print("number of particles")
  print(n.particles)
  print("marg-logLike's:")
  
  for(i in 1:n.replicates){
    sample.log.like[i] <- particleFilter(larvalModP,
                                         theta(),
                                         init.state,
                                         data = garkiObs,
                                         n.particles)
    print(sample.log.like[i])
  }
  end.time  <- Sys.time()
  
  sample.finite.log.like <- sample.log.like[is.finite(sample.log.like)]
  
  ans <- c(mean = mean(sample.finite.log.like), 
           sd = sd(sample.finite.log.like), 
           prop.depletion = 1-length(sample.finite.log.like)/length(sample.log.like), 
           time = end.time - start.time)
  
  res <- rbind(res, t(ans))
}

write.table(res,"C:\\Users\\ALM210\\Dropbox (SPH Imperial College)\\pFilterResults2.csv")
write.table(test.n.particles,"C:\\Users\\ALM210\\Dropbox (SPH Imperial College)\\pFilterCounts2.csv")



sd<-ggplot()+
  xlab("Number of Particles")+
  ylab("S.D. of Log Liklihood")+
  geom_point(aes(x=test.n.particles, y=res$sd), color='red',alpha=0.9)+
  geom_smooth(aes(x=test.n.particles, y=res$sd), color='blue',alpha=0.3)+
  theme_bw()


mean<-ggplot()+
  xlab("Number of Particles")+
  ylab("Mean of Log Liklihood")+
  ggtitle("100 Iterations of particle filter with fixed parameter values")+
  geom_point(aes(x=test.n.particles, y=res$mean), color='red',alpha=0.9)+
  geom_smooth(aes(x=test.n.particles, y=res$mean), color='blue',alpha=0.3)+
  theme_bw()


res2$time2<-res$time*10000
res2$time2<-((res2$time2/60)/60)/24
res2$time2[13]<-7.9

time<-ggplot()+
  xlab("Number of Particles")+
  ylab("Time in Days for 10,000 iter. MCMC")+
  geom_point(aes(x=test.n.particles, y=res2$time2), color='red',alpha=0.9)+
  geom_smooth(aes(x=test.n.particles, y=res2$time2), color='blue',alpha=0.3)+
  theme_bw()


grid.arrange(mean,sd)


mLLik1=function(x){particleFilter(larvalModP,theta(uoE=x),
                                                 init.state=init.state,
                                                 data = garkiObs,
                                                 n.particles=2000)}
x=seq(0.0001,0.15,length.out=1000)
y=sapply(x,mLLik1)

uoLplot<-ggplot()+
  xlab("uoL")+
  ylab("Marginal log Like")+
  geom_point(aes(x=x, y=y), color='red',alpha=0.9)+
  theme_bw()


mLLik1=function(x){particleFilter(larvalModP,theta(uoE=x),
                                  init.state=init.state,
                                  data = garkiObs,
                                  n.particles=2000)}
x=seq(0.01,0.1,length.out=500)
y=sapply(x,mLLik1)

uPplot<-ggplot()+
  xlab("uP")+
  ylab("Marginal log Like")+
  geom_point(aes(x=x, y=y), color='red',alpha=0.9)+
  theme_bw()