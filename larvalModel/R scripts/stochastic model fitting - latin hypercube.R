library(lhs)
library(lubridate)
library(ggplot2)
library(VGAM)
library(odin)
library(parallel)
library(devtools)




#########################################################################################################################################
#                                                                                                                                       #
#                                                     stochastic larval model                                                           #
#                                                                                                                                       #
#########################################################################################################################################

# 
# mos_params <- function(
#   #parameter estimates from white et al (2011)
#   rF=rFx,
#   dE = 0.150, #development time of early larval instars
#   dL = 0.269, #development time of late larval instars
#   dP = 1.563, #development time of pupae
#   uoE = 0.034, #per capita daily mortality rate of early instars (low density)
#   uoL = 0.035, #per capita daily mortality rate of late instars (low density)
#   uP = 0.25, #per capita daily mortality rate of pupae
#   uM = 0.096, # per capita daily mortality rate of adult An.gambiae - KEEP
#   B = 21.19, #No. of eggs laid per day per mosquito
#   Y = 13.25, #effect of density dependence on late instars relative to early instars
#   S = 3, #duration of gonotrophic cycle - KEEP
#   Emax = 93.6, #max number of eggs per oviposition per mosquito - KEEP
#   tr = 4, #days of rainfall contributing to carrying capacity
#   sf = 10, #scaling factor
#   dt=delta,
#   O=1,
#   n=10,
#   E0=177,
#   L0=8,
#   P0=1,
#   M0=7
# )
# return(as.list(environment()))

##model in odin
larvalR <- odin::odin({
  #param
  dE <- user()
  dL <- user()
  dP <- user()
  uoE <- user()
  uoL <- user()
  uP <- user()
  uM <- user()
  Y <- user()
  S <- user()
  tr <- user()
  sf <- user()
  dt<-user()
  n<-user()
  rF[] <- user()
  dim(rF) <- user()
  trx<-(tr/dt)
  Emax<-user()
  E0<-user()
  L0<-user()
  P0<-user()
  M0<-user()
  time1<-user()
  o<-user()
  Mg<-user()
  
  #initial values
  
  initial(E) <- E0
  initial(L) <- L0
  initial(P) <- P0
  initial(M) <- M0
  initial(timeX)<-time1
 
  tt<-round(trx)

  K <-if (timeX<=tt) ((sf*((1/tt)*(sum(rF[0:(timeX-1)]))))) 
  else ((sf*((1/tt)*(sum(rF[(timeX-tt):timeX-1])))))
  
  uE = uoE*(1+((E + L) / (K)))
  uL = uoL*(1+(Y*((E + L) / (K))))
  
  update(timeX)<-timeX+1
  
  Be<-if((dE+uE)*dt<1) rbinom(E,(dE+uE)*dt) else rbinom(E,1)
  Bl<-if ((dL+uL)*dt<1) rbinom(L, (dL+uL)*dt) else rbinom(L,1)
  Bp<-rbinom(P,(dP+uP)*dt)
  Bm<-rbinom(M,uM*dt)
  nt<-rbinom(M,(dt/S))

  update(E)<-round(E-Be+rpois(nt*n)) #else rpois(nt*n)
  update(L)<-round(L-Bl+rbinom(Be,(dE/(uE+dE)))) #else rbinom(Be,(dE/(uE+dE)))
  update(P)<-round(P-Bp+rbinom(Bl,(dL/(uL+dL)))) #else rbinom(Bl,(dL/(uL+dL)))
  mTemp =M+round(0.5*(rbinom(Bp,(dP/(uP+dP)))))-Bm+rpois(round(Mg))
  update(M)<-mTemp

})




#########################################################################################################################################
#                                                                                                                                       #
#                                                     parameterisation with LHC's                                                       #
#                                                                                                                                       #
#########################################################################################################################################

parms<-mosParamsP(uoE=0.09007434,uoL=0.03156214,uP=0.2499843,Y=11.501,
           o=7.827497,sf=57110,n=50,time1=1244,E0=4311, L0=35, P0=5, M0=42)
modR <- larvalR(user=parms)
simDat <- as.data.frame(modR$run(1244:2000))
summary(is.na(simDat$M))
plot(simDat$M)

##model simulation function - starts simulation from same random seed each time
modSim <- function(parms, obsDat) {
  modR <- larvalR(user=parms)
  simDat <- as.data.frame(modR$run(1244:3000))
  simDat<-simDat[seq(1, NROW(simDat), by = 1/delta),]
  simDat$step<-simDat$step*delta
  matchedTimes <- simDat$step %in% garkiObs$time
  simDat$M[simDat$M<=0|is.na(simDat$M)] = 1e-3
  return(simDat)
}
##goodness of fit function - log liklihood
GOF<-function(pr){ #add parameters
  parameters<-read.table(text = pr, sep = ",", colClasses = "numeric")
  print(parameters)
  modSim <- modSim(parms=mosParamsP(uoE=parameters$V1,uoL=parameters$V2,uP=parameters$V3,
                                    Y=parameters$V4,n=parameters$V5,sf=parameters$V6,
                                    E0=177,L0=9,P0=1,M0=7),obsDat = garkiObs)
  matchedTimes <- modSim$step %in% garkiObs$time
  log_like <- (sum(dzipois(garkiObs$M,modSim$M[matchedTimes], log = TRUE)))+ lprior(as.numeric(parameters))
  print(log_like)
  if (!is.na(log_like)){
    res<-c(parameters$V1,parameters$V2,parameters$V3,parameters$V4,parameters$V5,parameters$V6,log_like)
  }
  return(res)
}

pGOF<-function(pr){
  pr<-read.table(text = pr, sep = ",", colClasses = "numeric")
  res1<-pFilt(96,iState,modStep3,dataLikFunc,garkiObs101,pr=pr,rFclust=1,fxdParams=50)+lprior(as.numeric(pr))
  pr[9]<-pr[10]
  res2<-pFilt(96,iState,modStep3,dataLikFunc,garkiObs104,pr=pr,rFclust=1,fxdParams=50)+lprior(as.numeric(pr))
  pr[9]<-pr[11]
  res3<-pFilt(96,iState,modStep3,dataLikFunc,garkiObs219,pr=pr,rFclust=2,fxdParams=50)+lprior(as.numeric(pr))
  pr[9]<-pr[12]
  res4<-pFilt(96,iState,modStep3,dataLikFunc,garkiObs220,pr=pr,rFclust=2,fxdParams=50)+lprior(as.numeric(pr))
 # print(paste0("iteration ",pr[13]))
  print(sum(res1,res2,res3,res4))
  return(sum(res1,res2,res3,res4))
}

for(i in 1:500){
  pGOF(paste(parms[1],parms[2],parms[3],parms[4],parms[5],parms[6],parms[7],parms[8],parms[9],parms[10],parms[11],parms[12] ,sep=","))
}

HC<-NULL
##initialise latin hypercube
HCtemp<- randomLHS(10000, 12)
HC$count<-c(1:10000)
HC$uoE <-qnorm(HCtemp[,1],mean=0.035,sd=0.007)
HC$uoL <-qnorm(HCtemp[,2],mean=0.035,sd=0.007)
HC$uP <- qnorm(HCtemp[3],mean=0.25,sd=0.0457)
HC$Y <- qnorm(HCtemp[4],mean=13.06,sd=4.53)
HC$p0<-qnorm(HCtemp[5],mean=0.5,sd=0.05)
HC$Lx<-(HCtemp[,6])
HC$Fp<-qunif(HCtemp[,7],min=0.8,max=1)
HC$o<-qnorm(HCtemp[,8],mean=1,sd=0.05)
HC$sf1<-(15*HCtemp[,9])
HC$sf2<-(15*HCtemp[,10])
HC$sf3<-(15*HCtemp[,11])
HC$sf4<-(15*HCtemp[,12])


HC<-as.data.frame(HC)



HC<- paste(HC$uoE,HC$uoL,HC$uP, HC$Y,HC$p0,HC$Lx,HC$Fp,HC$o,HC$sf1,HC$sf2,HC$sf3, HC$sf4,HC$count, sep=",")#concatonate into single vector 

clusterExport(cl, c("pGOF","garkiObs101","garkiObs104","garkiObs219","garkiObs220"), envir=environment())


system.time(lhcResults2<-lapply(HC,pGOF))#convert results to dataframe

colnames(lhcResults)<-c("uoE","uoL","uP","Y","n","sf","logLike")

mins<-lhcResults[which(lhcResults$logLike == max(lhcResults$logLike)), ]#find minimum values


#"0.0275392720335688,0.0391684551658979,0.317601558489952,21.8705305912987,0.474725347165231,0.162143399122823,0.921936883796244,1.07328942136458,2.62381167667732,7.33830458177731,4.35943666368211,4.29330090065394,4392"

HC[which(grepl(max(g), g))]

#pr<-paste(mins$uoE,mins$uoL,mins$uP, mins$sf,mins$Y, mins$n, sep=",")#concatonate into single vector 


##run and plot model with set parameter values
set.seed(10) 
mod <- odinPackage::larvalModP(user=mos_params(uoE=mins$uoE,uoL=mins$uoL,uP=mins$uP,
                               Y=mins$Y,sf=mins$sf,n=mins$n, E0=177,L0=9,P0=1,M0=7)) #parameters estimated from LHC sampling
sim <- as.data.frame(mod$run(0:2000))

df<-sim[seq(1, NROW(sim), by = 1/delta),]
df$time<-df$step*delta

ggplot(data=df, aes(time))+ 
  geom_line(data=df,aes(x=time, y=M), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=garkiObs,(aes(x=time,y=M)),col="red")+
  theme_bw()


###make simulated data
resSim<-c(0:3000)
for (i in 1:1000){
  
  mod <- larvalR(user=mos_params(uoE=0.03437649,uoL=0.03540515,uP=0.24770876,
                                                 Y=7.43381088,n=2,sf=2, E0=177,L0=9,P0=1,M0=7)) #parameters estimated from LHC sampling
  sim <- as.data.frame(mod$run(0:1650))
  plot(sim$M,col="white")
  lines(sim$M)
  points(garkiObs101$time*10,garkiObs101$`101`,col="red")
  lines(rFx)
  resSim<- cbind(resSim,sim$M)
  
}

simDat2<-as.data.frame(c(1:20*10))
colnames(simDat2)<-c("time")
resSim<-as.data.frame(resSim)
resSim<-resSim[ , colSums(is.na(resSim)) == 0]
resSimMean<-rowMeans(resSim[,-1])
resSimMean[1:20*100]

simDat2$M<-resSimMean[1:20*100]

simDat2<-rbind(data.frame(time = 0, M = 0), simDat2)

#


garkiObs<-rbind(garkiObs,data.frame(time = 202, M = 7))

