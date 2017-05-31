


# define a sampler for the prior on the initial state
simx0 <- function(N,t0)
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
pFilt <- function (n, simx0, t0, stepFun, dataLik, data,pr) 
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
    w[is.na(w)] <- 1e-200##only keep in if needed
    rows = sample(1:n, n, replace = TRUE, prob = w)
    xmat = xmat[rows, ]
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
#uoE         uoL        uP          Y  n      sf   
# 0.03188414 0.005738873 0.2690644 13 20 0.9429963 

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