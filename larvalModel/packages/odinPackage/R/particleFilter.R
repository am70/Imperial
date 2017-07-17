
##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter functions                                                                   #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

# initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M, should fitted value be random?
iState <- function(N,t0)
{
  mat=cbind(rpois(N,177),rpois(N,8),rpois(N,1),rpois(N,7))
  colnames(mat)=c("E","L","P","M")
  mat
}

##negative binomial - gamma poisson (mixture) distribution
nBgP<-function(x,n,p,log=F){
  if(x==0 && n==0) 1
  else if (n==0) 0
  else if(log==T)(lgamma(x+n)-(lgamma(n)+lfactorial(x)))+(log(p^n)+log((1-p)^x))
  else  exp(lgamma(x+n)-(lgamma(n)+lfactorial(x))+log(p^n)+log((1-p)^x))
}



#likelihood function
dataLikFunc <- function(input)
{
  dataInput<-read.table(text = input, sep = ",", colClasses = "numeric")
  #dataInput[1]=Simulated data point
  #dataInput[2]=Observed data point
  ll = nBgP(as.numeric(round(dataInput[1],0)), as.numeric(dataInput[2]+1), p=as.numeric(dataInput[3]), log = T)
  return (ll)
}



##model step function - runs model in steps using Odin, returning the model state at a subsequent time point
modStep3<-function(weightInput){
  #parse input data
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  initialState<-as.numeric(dataInput[,c(1:4)])#state of E, L, P and M
  currentTime<-as.numeric(dataInput[,c(5)])#time in model to run from
  nextTime<-as.numeric(dataInput[,c(6)])#time in model to run to
  pr<-as.numeric(dataInput[,c(7:12)])#parameters
  #initialise model
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],sf=pr[6])
  #run model between two discrete time periods and return results
  modR <- larvalModP(user=params)#run model
  simDat <- as.data.frame(modR$run(seq(currentTime, nextTime, length.out=nextTime-currentTime)))
  res2<-(simDat[simDat$step == nextTime,])#return model state
  return(as.data.frame(rbind(res2[,c(8:11)])))
}



##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter                                                                             #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################



pFilt <- function (n, iState, t0, stepFun, likeFunc, obsData,prms,resM=F) 
{
  times = c(obsData$time/delta) #/delta as model is running in discrete time steps
  particles = iState(n, t0) #initial state
  ll = 0
  for (i in 1:length(times[-length(times)])) {
    wp<-paste(particles[,1],particles[,2],particles[,3], particles[,4],times[i],
              times[i + 1],prms[1],prms[2],prms[3],prms[4],prms[5],prms[7],sep=",")
    particlesTemp = parLapply(cl,wp,stepFun) #use NULL for dide cluster, cl for local
    particles<-data.frame(t(sapply(particlesTemp, `[`)))
    likeDat<-paste(particles$M,obsData[i+1,2],prms[6],sep=",")
    weights = parLapply(cl,likeDat,likeFunc)
    
    weights<-as.vector(unlist(weights))
    if(resM==T)ll=rbind(ll,mean(as.numeric(particles$M)))
    else ll = ll + mean(weights)
    #normalise weights
    swP=sum(weights)
    weights=weights/swP
    weights<-(weights)-max(weights)
    weights<-(exp(weights)*0.9)
    rows = sample(1:n, n, replace = TRUE, prob = weights)
    particles = particles[rows, ]
  }
  
  ll
  
}

# llx<-ll[,sample(1:length(tail(ll,1)),1,prob=exp(head(tail(ll,2),1)))]
# sum(llx[-1])

