
##Zero inflated poisson density function
ziP<-function(x,lambda,p0,log=FALSE){
  if(log==FALSE){
    ifelse(x==0,p0+(1-p0)*dpois(x,lambda),(1-p0)*dpois(x,lambda))
  }
  else{
    ifelse(x==0,log(p0+(1-p0)*dpois(x,lambda)),log((1-p0)*dpois(x,lambda)))
  }
}


# define a sampler for the prior on the initial state
iState <- function(N,t0)
{
  mat=cbind(rpois(N,177),rpois(N,8),rpois(N,1),rpois(N,7))
  colnames(mat)=c("E","L","P","M")
  mat
}


dataLik2 <- function(input)
{
  dataInput<-read.table(text = input, sep = ",", colClasses = "numeric")
  ll = dzipois(dataInput[2], dataInput[1], pstr0=0, log = T)
  return (exp(ll))
}


###Particle Filter
pFilt <- function (n, iState, t0, stepFun, dataLik, obsData,prms) 
{
  times = c(as.numeric(obsData$time))
  particles = iState(n, t0) #initial state
  ll = 0
  for (i in 1:length(times[-length(times)])) {
    wp<-paste(particles[,1],particles[,2],particles[,3], particles[,4],times[i],
              times[i + 1],prms[1],prms[2],prms[3],prms[4],prms[5],prms[6],sep=",")
    
    particlesTemp = parLapply(cl,wp,stepFun) #use NULL for dide cluster, cl for local
    particles<-data.frame(t(sapply(particlesTemp, `[`)))
    
    likeDat<-paste(particles$M,obsData[i,2],sep=",")
    weights = lapply(likeDat,dataLik)
    weights<-as.vector(unlist(weights))
    
    swP=sum(weights)
    weights=weights/swP
    
    weights<-log(weights)-max(log(weights))
    weights<-(exp(weights))
    
    ll = ll + log(mean(weights))
    weights[is.na(weights)] <- 1e-200##only keep in if needed
    rows = sample(1:n, n, replace = TRUE, prob = weights)
    particles = particles[rows, ]
  }
  
  ll
  
}



theta <-
  function(
    uoE = 0.1736114,uoL = 0.4792658,uP = 0.2163627,y = 0.1301954,n =  2.932363,sf = 13.74987)
    return(c(uoE=uoE,uoL=uoL,uP=uP,y=y,n=n,sf=sf))



init.state <-
  c(E = 177,L = 8,P = 1,M = 7)

###test parms
#uoE        uoL        uP  Y  n       sf   
# 0.05985359 0.01074641 0.3381503 13 20 4.207474 

###

## Log-Prior 
lprior <- function(parms) with(parms, {
  uoEprior<-dnorm(parms$uoE,mean=0.035,sd=0.00485,log=T)
  uoLprior<-dnorm(parms$uoL,mean=0.035,sd=0.00485,log=T)
  uPprior<-dnorm(parms$uP,mean=0.25,sd=0.0357,log=T)
  Yprior<-dnorm(parms$Y,mean=13.06,sd=3.53,log=T)
  sfprior<-dunif(parms$sf,min=0,max=100,log=T)
  priorSum<-(uoEprior+uoLprior+uPprior+Yprior+sfprior)
  return(priorSum)
})