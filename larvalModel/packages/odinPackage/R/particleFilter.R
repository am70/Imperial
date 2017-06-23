
##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter functions                                                                   #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

# initial state sampler, samples random initial states
iState <- function(N,t0)
{
  mat=cbind(rpois(N,177),rpois(N,8),rpois(N,1),rpois(N,7))
  colnames(mat)=c("E","L","P","M")
  mat
}

#likelihood function
dataLikFunc <- function(input)
{
  dataInput<-read.table(text = input, sep = ",", colClasses = "numeric")
  ll = dnbinom(as.numeric(dataInput[2]), as.numeric(dataInput[1]), prob=as.numeric(dataInput[3]), log = T)
  return (ll)
}


##model step function - runs model in steps taking the likelihood of the subsequent observed data point
modStep3<-function(weightInput){
  #parse input data
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  initialState<-as.numeric(dataInput[,c(1:4)])
  currentTime<-as.numeric(dataInput[,c(5)]*10)#*10 as model is running in discrete time steps
  nextTime<-as.numeric(dataInput[,c(6)]*10)
  pr<-as.numeric(dataInput[,c(7:12)])
  #initialise model
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],sf=pr[6])
  #run model between two discrete time periods and return results
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(currentTime, nextTime, length.out=nextTime-currentTime)))
  res2<-(simDat[simDat$step == nextTime,])
  return(as.data.frame(rbind(res2[,c(8:11)])))
}



##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter                                                                             #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################


pFilt <- function (n, iState, t0, stepFun, dataLik, obsData,prms) 
{
  times = c(obsData$time)
  particles = iState(n, t0) #initial state
  ll = 0
  for (i in 1:length(times[-length(times)])) {
    wp<-paste(particles[,1],particles[,2],particles[,3], particles[,4],times[i],
              times[i + 1],prms[1],prms[2],prms[3],prms[4],prms[5],prms[6],sep=",")
    
    particlesTemp = parLapply(cl,wp,stepFun) #use NULL for dide cluster, cl for local
    particles<-data.frame(t(sapply(particlesTemp, `[`)))
    
    likeDat<-paste(particles$M,obsData[i+1,2],prms[7],sep=",")
    weights = lapply(likeDat,dataLik)
    weights<-as.vector(unlist(weights))

    ll = ll + mean(weights)
  
    #normalise weights
    swP=sum(weights)
    weights=weights/swP
    weights<-(weights)-max((weights))
    weights<-exp(weights)
    
    weights[is.na(weights)] <- 1e-4##only keep in if needed
    
    rows = sample(1:n, n, replace = TRUE, prob = weights)
    particles = particles[rows, ]
  }
  
  ll
  
}




