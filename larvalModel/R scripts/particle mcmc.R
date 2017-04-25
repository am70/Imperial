#########################################################################################################################################
#                                                                                                                                       #
#                                                     particle MCMC                                                                     #
#                                                                                                                                       #
#########################################################################################################################################
require(smfsb)

#x0,t0,delt,th
modStep<-function(dat){
  dat<-read.table(text = dat, sep = ",", colClasses = "numeric")
  colnames(dat)<-c("xmat","times","deltas","uoE","uoL","uP","sf","Y","n")
  params<-mosParamsP(M0=dat$xmat,uoE=dat$uoE,uoL=dat$uoL,uP=dat$uP,n=dat$n,sf=dat$sf,Y=dat$Y)
  modR <- larvalModP(user=params)
  simDat <- as.data.frame(modR$run(seq(0, 2000, length.out=201)))
  result<-(simDat[simDat$step/10 == dat$times+dat$deltas,])
  return(unlist(result[,c(11)]))
}

pfMLLik <- function (n, simx0, t0, stepFun, dataLik, obs) 
{
  times = c(t0, as.numeric(rownames(obs)))
  deltas = diff(times)
  return(function(...) {
    xmat = simx0(n, t0, ...)
    ll = 0
    for (i in 1:length(deltas)) { #for each time step, run each particle
      dat<-cbind(xmat,times[i],deltas[i],th[1],th[2],th[3],th[4],th[5],th[6])
      dat<-paste(dat[,1],dat[,2],dat[,3],dat[,4],dat[,5],dat[,6],dat[,7],dat[,8],dat[,9],sep=",")
      xmat<-data.frame(as.vector(t(sapply(parLapply(cl,dat,modStep), `[`))))
      
            w = apply(xmat, 1, dataLik, t = times[i + 1], y = obs[i,], log = FALSE, ...)
     
      ll = ll + log(mean(w))
      rows = sample(1:n, n, replace = TRUE, prob = w)
      xmat = as.data.frame(as.vector(xmat[rows, ]))
      colnames(xmat)<-c("M")
    }
    
    ll
  })
}

# set up data likelihood
dataLik <- function(x,t,y,log=TRUE,...)
{
  ll=sum(dzipois(y[2],x,log=TRUE))
  if (log)
    return(ll)
  else
    return(exp(ll))
}

##initial conditions
simx0 <- function(N,t0,...)
{
  mat=cbind(rpois(N,10))
  colnames(mat)=c("M")
  mat
}


garkiObsMat<-as.matrix(garkiObs)#convert observed data to matrix
# create marginal log-likelihood functions, based on particle filter
mLLik=pfMLLik(100,simx0,0,modStep,dataLik,garkiObsMat)
iters=10000
tune=0.01
thin=1
th=c(uoE = 0.104366, uoL = 0.015119, uP = 0.171963, sf=6.630103, Y=13.76813, n=8.059338)
p=length(th)
ll=-1e99
thmat=matrix(0,nrow=iters,ncol=p)
colnames(thmat)=names(th)
# Main pMCMC loop
for (i in 1:iters) {
  message(paste(i,""),appendLF=FALSE)
  for (j in 1:thin) {
    thprop=th*exp(rnorm(p,0,tune))
    llprop=mLLik(thprop)
    if (log(runif(1)) < llprop - ll) {
      th=thprop
      ll=llprop
    }
  }
  thmat[i,]=th
}

#plot summaries
mcmcSummary(thmat)

