
##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter functions                                                                   #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

# initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M

iState<-function(N,t,prms){
parms<-mosParamsP(uoE=prms[1],uoL = prms[2],uP=prms[3],Y=prms[4],n=prms[5],sf=prms[10],o=prms[9])
  o<-prms[7]
  dE<-parms$dE
  dL<-parms$dL
  dP<-parms$dP
  y<-parms$Y
  S<-parms$S
  sf<-parms$sf
  n<-parms$n
  UoE<-parms$uoE
  UoL<-parms$uoL
  rF<-parms$rF
  Um<-parms$uM
  Up<-parms$uP
  tr<-parms$tr/delta
  

  
  K<-(1+(sf*((1/tr)*(sum(rF[(t-tr):t-1])))))
  
  a=((n/S)*dP*dL)/((2*Um)*(Up*dP))
  b=(UoE/(y*UoL))*(dL+UoL)-dE-UoE
  c=-(UoE*dE)/(UoL*y)
  x=(-b+sqrt(b^2*-4*a*c))/(2*a)
  
  L<-o#((dE*x-dL-UoL)/(UoL*y))*(1/o)*(K/(x+1))#L
  E<-L/x
  P<-(dL*L)/(Up+dP)
  M<-(dP*P)/(2*Um)
  conds<-cbind(E,L,P,M)
  conds<-conds[rep(seq_len(nrow(conds)), each=N),]
  return(round(conds,0))
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
  #dataInput[3]=negBin prob
  #dataInput[4]=population scaling factor
  #dataInput[5]=time in model for rainfall calc
  #rain<-if (dataInput[5]<140) sum(rFx[1:(as.numeric(dataInput[5])-1)]) else sum(rFx[(as.numeric(dataInput[5])-140):(as.numeric(dataInput[5])-1)])
  fracPop<-dataInput[2]/(dataInput[4]+1e-5)#*(rain/(dataInput[2]+1))
  ll = nBgP(as.numeric(round(dataInput[1],0)), as.numeric(1+fracPop), p=as.numeric(dataInput[3])+1e-10, log = T)##+1e-10 to stop INF if MH proposes a 0
  return (ll)
}



##model step function - runs model in steps using Odin, returning the model state at a subsequent time point
modStep3<-function(weightInput){
  #parse input data
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  initialState<-as.numeric(dataInput[,c(1:4)])#state of E, L, P and M
  currentTime<-as.numeric(dataInput[,c(5)])#time in model to run from
  nextTime<-as.numeric(dataInput[,c(6)])#time in model to run to
  pr<-as.numeric(dataInput[,c(7:13)])#parameters
  #initialise model
  if(dataInput[14]==1) rFc<-rFx  else rFc<-rFx2
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],o=pr[6],sf=pr[7],rF=rFc,time1=currentTime)
  #run model between two discrete time periods and return results
  modR <- larvalModP(user=params)#run model
  simDat <- as.data.frame(modR$run(seq(currentTime, nextTime, length.out=nextTime-currentTime)))
  res2<-(simDat[simDat$step == nextTime,])#return model state
  return(as.data.frame(rbind(res2[,c(8:11)])))
}


##model step function - runs model in steps using Odin, returning the model state at a subsequent time point
modStep4<-function(weightInput){
  #parse input data
  dataInput<-read.table(text = weightInput, sep = ",", colClasses = "numeric")
  initialState<-as.numeric(dataInput[,c(1:4)])#state of E, L, P and M
  currentTime<-as.numeric(dataInput[,c(5)])#time in model to run from
  nextTime<-as.numeric(dataInput[,c(6)])#time in model to run to
  pr<-as.numeric(dataInput[,c(7:13)])#parameters
  #initialise model
  if(dataInput[14]==1) rFc<-rFx  else rFc<-rFx2
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],n=pr[5],o=pr[6],sf=pr[7],rF=rFc,time1=currentTime)
  #run model between two discrete time periods and return results
  modR <- larvalModP(user=params)#run model
  simDat <- as.data.frame(modR$run(seq(currentTime, nextTime)))

    return(as.data.frame(simDat))
}


##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter                                                                             #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################



pFilt <- function (n, iState, stepFun, likeFunc, obsData,prms,resM=F,rFclust) 
{
  times = c(obsData$time/delta) #/delta as model is running in discrete time steps

  particles = iState(n, t = times[1],prms = prms) #initial state
  ll = 0
  rMeans<-NULL
  for (i in 1:length(times[-length(times)])) {
    wp<-paste(particles[,1],particles[,2],particles[,3], particles[,4],times[i],
              times[i + 1],prms[1],prms[2],prms[3],prms[4],prms[5],prms[9],prms[10],rFclust,sep=",")
    
    #if statement for if outputting model results for plotting rather than LL

    particlesTemp = parLapply(cl,wp,stepFun)  #use NULL for dide cluster, cl for local
    particles<-data.frame(t(sapply(particlesTemp, `[`)))
    
    likeDat<-paste(particles$M,obsData[i+1,2],prms[6],prms[8],times[i],sep=",")
    weights = parLapply(cl,likeDat,likeFunc)
    
    weights<-as.vector(unlist(weights))
     ll = ll + mean(weights)
    #normalise weights
    swP=sum(weights)
    weights=weights/swP
    weights<-(weights)-max(weights)
    weights<-(exp(weights)*0.9)
    weights[is.na(weights)]<-1e-50
    
    #used for outputting full runs from pF for model visualisations. 
    if(resM==T){
      particlesTemp1 = parLapply(cl,wp,modStep4)
      pt<-NULL
      for(j in 1:n){
        M<-(particlesTemp1[[j]]$M)
        pt<-cbind(pt,M)
      }
      pt<-as.data.frame(rowMeans(pt))
      rMeans<-rbind(rMeans,pt)
    }
    
    
    rows = sample(1:n, n, replace = TRUE, prob = weights)
    particles = particles[rows, ]
  }
  
  if(resM==T) return (rMeans) 
  else return (ll)
  
}



