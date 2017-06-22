

# initial state sampler
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


## Log-Prior 
lprior <- function(parms) {
  uoEprior<-dnorm(parms[1],mean=0.035,sd=0.00485,log=T)
  uoLprior<-dnorm(parms[2],mean=0.035,sd=0.00485,log=T)
  uPprior<-dnorm(parms[3],mean=0.25,sd=0.0357,log=T)
  Yprior<-dnorm(parms[4],mean=13.06,sd=3.53,log=T)
  sfprior<-dunif(parms[6],min=0,max=100,log=T)
  p0prior<-dnorm(parms[7],mean=0.5,sd=0.03,log=T)
  priorSum<-as.vector(uoEprior+uoLprior+uPprior+Yprior+sfprior+p0prior)
  return(priorSum)
}


###Particle Filter
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


init.state <-
  c(E = 177,L = 8,P = 1,M = 7)


