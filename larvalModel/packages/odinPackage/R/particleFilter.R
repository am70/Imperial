
##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                    particle Filter functions                                                                   #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

# initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M
# @param N number of particles for particle filter
# @param t starting time period for rainfall
# @param prms model parameters
# @fxdParams fixed parameter (currently just for n in mosParamsP) - maybe update for multiple fixed parameters
# @return conditions for E, L, P and M
iState<-function(N,t,prms,fxdParams){
parms<-mosParamsP(uoE=prms[1],uoL = prms[2],uP=prms[3],Y=prms[4],n=as.numeric(fxdParams),sf=prms[8],o=prms[7])
Lx <- prms[5]
dE <- parms$dE
dL <- parms$dL
dP <- parms$dP
y <- parms$Y
S <- parms$S
sf <- parms$sf
n <- parms$n
UoE <- parms$uoE
UoL <- parms$uoL
rF <- parms$rF
Um <- parms$uM
Up <- parms$uP
tr <- parms$tr / delta

K <- (1 + (sf * ((1 / tr) * (sum(
  rF[(t - tr):t - 1]
)))))

a = ((n / S) * dP * dL) / ((2 * Um) * (Up * dP))
b = (UoE / (y * UoL)) * (dL + UoL) - dE - UoE
c = -(UoE * dE) / (UoL * y)
x = (-b + sqrt(b ^ 2 * -4 * a * c)) / (2 * a)

L <- Lx#((dE*x-dL-UoL)/(UoL*y))*(1/o)*(K/(x+1))#L
E <- L / x
P <- (dL * L) / (Up + dP)
M <- (dP * P) / (2 * Um)
conds <- cbind(E, L, P, M)
conds <- conds[rep(seq_len(nrow(conds)), each = N), ]
return(round(conds, 0))
}



##negative binomial - gamma poisson (mixture) likelihood function
# @param x vector of quantiles
# @param n target number of succesful trials
# @param p probability of success in each tria
# @param log True or False return results in log
# @return likelihood value using negative binomial - gamma poisson (mix)
nBgP <- function(x, n, p, log = F) {
  p[p > 1] <- 1
  if (x == 0 && n == 0)
    1
  else if (n == 0)
    0
  else if (log == T)
    (lgamma(x + n) - (lgamma(n) + lfactorial(x)))+(log(p^n)+log((1-p)^x))
  else  exp(lgamma(x+n)-(lgamma(n)+lfactorial(x))+log(p^n)+log((1-p)^x))
}

# Beta binomial likelihood function
# @param k Observed data point
# @param n Simulated data point
# @param p fraction of population
# @param w overdisperion parameter
# @return log likelihood
betaBinom <- function(k, n, p, w) {
  a <- p * ((1 / w) - 1)
  b <- (1 - p) * ((1 / w) - 1)
  lbeta(k + a, n - k + b) - lbeta(a, b) + lchoose(n, k)
}



#likelihood function - obsolete?
# @param input string of values for likelihood function - in order: sim data point, Obs data point, negBin prob, population scaling factor, time in model for rainfall 
# @return returns likelihood value
dataLikFunc <- function(input)
{
  dataInput <-
    read.table(text = input,
               sep = ",",
               colClasses = "numeric")
  #dataInput[1]=Simulated data point
  #dataInput[2]=Observed data point
  #dataInput[3]=negBin prob
  #dataInput[4]=population scaling factor
  #dataInput[5]=time in model for rainfall calc
  #rain<-if (dataInput[5]<140) sum(rFx[1:(as.numeric(dataInput[5])-1)]) else sum(rFx[(as.numeric(dataInput[5])-140):(as.numeric(dataInput[5])-1)])
  fracPop <- dataInput[2] / (dataInput[4] + 1e-5)#*(rain/(dataInput[2]+1))
  ll = nBgP(
    as.numeric(round(dataInput[1], 0)),
    as.numeric(1 + fracPop),
    p = as.numeric(dataInput[3]) + 1e-10,
    log = T
  )##+1e-10 to stop INF if MH proposes a 0
  return (ll)
}



##model step function - runs model in steps using Odin, returning the model state at a subsequent time point
# @ param weightInput string of values for model state, time to run model from/to and model parameters
modStep3 <- function(weightInput) {
  #parse input data
  dataInput <-
    read.table(text = weightInput,
               sep = ",",
               colClasses = "numeric")
  initialState <-
    as.numeric(dataInput[, c(1:4)])#state of E, L, P and M
  currentTime <-
    as.numeric(dataInput[, c(5)])#time in model to run from
  nextTime <- as.numeric(dataInput[, c(6)])#time in model to run to
  pr <- as.numeric(dataInput[, c(7:12)])#parameters
  fixed <- dataInput[14]#fixed parameters, currently just n
  #initialise model
  if (dataInput[13] == 1)
    rFc <- rFx
  else
    rFc <- rFx2
  params<-mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],o=pr[5],sf=pr[6],rF=rFc,time1=currentTime,n=as.numeric(fixed))
  #run model between two discrete time periods and return results
  modR <- larvalModP(user = params)#run model
  simDat <-
    as.data.frame(modR$run(seq(
      currentTime, nextTime, length.out = nextTime - currentTime
    )))
  res2 <- (simDat[simDat$step == nextTime, ])#return model state
  return(as.data.frame(rbind(res2[, c(8:11)])))
}


##model step function - for use when returning p.filter full results
# @ param weightInput string of values for model state, time to run model from/to and model parameters
modStep4 <- function(weightInput) {
  #parse input data
  dataInput <-
    read.table(text = weightInput,
               sep = ",",
               colClasses = "numeric")
  initialState <-
    as.numeric(dataInput[, c(1:4)])#state of E, L, P and M
  currentTime <-
    as.numeric(dataInput[, c(5)])#time in model to run from
  nextTime <- as.numeric(dataInput[, c(6)])#time in model to run to
  pr <- as.numeric(dataInput[, c(7:12)])#parameters
  fixed <- dataInput[14]#fixed parameters, currently just n
  #initialise model
  if (dataInput[13] == 1)
    rFc <- rFx
  else
    rFc <- rFx2
  params<-
    mosParamsP(E0=initialState[1],L0=initialState[2],P0=initialState[3],M0=initialState[4],
                     uoE=pr[1],uoL=pr[2],uP=pr[3],Y=pr[4],o=pr[5],sf=pr[6],rF=rFc,time1=currentTime,n=as.numeric(fixed))
  #run model between two discrete time periods and return results
  modR <- larvalModP(user = params)#run model
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


# Particle filter function - runs model in steps between each observed data point 
# before using weighted re-sampling of model state at each step to start next
# @param n number of particles 
# @param iState function for calculating initial state of particles
# @param stepFun function for taking parmeters and time point, running model in steps and outputting model state
# @param likeFunc likelihood function - obsolete
# @param obsData observed data to fit model to
# @param prms model parameters
# @param resM True/False whether to output likelihood values or results of simulation for plotting figures
# @param fxedParams fixed parameters - currently just accepting n but may expand to include more than one parameter
# @param cluster T/F if true, using the dide cluster if false, run locally
# @return if resM = F log likelihood value, if resM = T mean results of simulation
pFilt <- function (n, iState, stepFun, likeFunc, obsData,prms,resM=F,rFclust,fxdParams, cluster=F) 
{
  times = c(obsData$time/delta) #/delta as model is running in discrete time steps

  particles = iState(n, t = times[1],prms = prms,fxdParams) #initial state
 ll= mean(betaBinom(obsData[1,2],particles[,4],0.01, prms[6]))#starts ll with betaBinom of starting values, else losing data from first data point
  rMeans<-NULL
  for (i in 1:length(times[-length(times)])) {
    wp<-paste(particles[,1],particles[,2],particles[,3], particles[,4],times[i],
              times[i + 1],prms[1],prms[2],prms[3],prms[4],prms[7],prms[8],rFclust,fxdParams,sep=",")
    
    
    
    if (cluster == F)
      particlesTemp = parLapply(cl, wp, stepFun) #use NULL for dide cluster, cl for local
    else
      particlesTemp = parLapply(NULL, wp, stepFun)
    particles <- data.frame(t(sapply(particlesTemp, `[`)))
    
    # fracPop<-obsData[i+1,2]/prms[7]
    weights <-
      betaBinom(obsData[i + 1, 2], (1 + round(as.vector(
        unlist(particles$M, 0)
      ))), 0.01, prms[6])#dnbinom(round(as.vector(unlist(particles$M,0))),round(1+fracPop,0),prms[5],log=T)#as.vector(unlist(weights))
    ll = ll + mean(weights)
    #normalise weights
    swP = sum(weights)
    weights = weights / swP
    weights <- (weights) - max(weights)
    weights <- (exp(weights) * 0.9)
    weights[is.na(weights)] <- 1e-50
    
    #used for outputting full runs from pF for model visualisations.
    if (resM == T) {
      if (cluster == F)
        particlesTemp1 = parLapply(cl, wp, modStep4)
      else
        particlesTemp1 = parLapply(NULL, wp, modStep4)
      pt <- NULL
      for (j in 1:n) {
        Mt <- (particlesTemp1[[j]]$M)
        pt <- cbind(pt, Mt)
      }
      pt <- as.data.frame(rowMeans(pt))
      rMeans <- rbind(rMeans, pt)
    }
    
    
    rows = sample(1:n, n, replace = TRUE, prob = weights)
    particles = particles[rows,]
  }
  
  if (resM == T)
    return (rMeans)
  else
    return (ll)
  
}
