library(lhs)
library(lubridate)
library(ggplot2)
library(VGAM)
library(odin)
library(parallel)
library(devtools)

setwd("C:\\Users\\Aaron\\Documents\\Imperial\\larval model\\data")




#########################################################################################################################################
#                                                                                                                                       #
#                                                     load in data                                                                      #
#                                                                                                                                       #
#########################################################################################################################################

##load in Garki rainfall data
rainfall<-read.csv("garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
rainfall$time<-(1:nrow(rainfall))

#load in mosquito data
garki<-read.table("spraycollect.csv",sep=",",head=T)
garki$ag_sum<-rowSums(garki[,c(8:16)])#create sum of A.gambiae spray samples
garki$Date<-as.Date(garki$Date)
#garki$Date<-dmy(as.character(garki$Date))
garki72<-subset(garki,Date >= as.Date("1973-05-27") & Date <= as.Date("1974-06-03"))
garki72_101<-subset(garki72,E_Station==101) ##subset by village - needs checking with michael specific villagse
garki72_101<-merge(garki72_101,rainfall,by.x="Date",by.y="date",all=T) #merge rainfall and mosquito data
garkiObs<-garki72_101[,c(39,37)]
colnames(garkiObs)<-c("time","M")
garkiObs<-subset(garkiObs,M>=0)

#########################################################################################################################################
#                                                                                                                                       #
#                                                     stochastic larval model                                                           #
#                                                                                                                                       #
#########################################################################################################################################
delta<-0.01#discrete time period
rF<-as.data.frame(rainfall$rainfall)
rFx<-rF[rep(seq_len(nrow(rF)), each=1/delta),] #split rainfall data into discreet time periods


mos_params <- function(
  #parameter estimates from white et al (2011)
  rF=rFx,
  dE = 0.150, #development time of early larval instars
  dL = 0.269, #development time of late larval instars
  dP = 1.563, #development time of pupae
  uoE = 0.034, #per capita daily mortality rate of early instars (low density)
  uoL = 0.035, #per capita daily mortality rate of late instars (low density)
  uP = 0.25, #per capita daily mortality rate of pupae
  uM = 0.096, # per capita daily mortality rate of adult An.gambiae
  B = 21.19, #No. of eggs laid per day per mosquito
  Y = 10, #effect of density dependence on late instars relative to early instars
  S = 3, #duration of gonotrophic cycle
  Emax = 93.6, #max number of eggs per oviposition per mosquito
  tr = 4, #days of rainfall contributing to carrying capacity
  sf = 64, #scaling factor
  dt=delta,
  O=1,
  n=10 #clutch size
)
return(as.list(environment()))

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
  
  
  #initial values
  initial(Be)<-0
  initial(Bl)<-0
  initial(Bp)<-0
  initial(Bm)<-0
  initial(nt)<-0
  
  initial(E) <- E0
  initial(L) <- L0
  initial(P) <- P0
  initial(M) <- M0
  
  initial(Reff)<-0
  
  K <-if (step<=trx) 1+(sf*(1/(trx*(1-exp(-step/trx)))))else 1+(sf*(1/(trx*(1-exp(-step/trx)))))*(sum(rF[(step-trx):step]))

  uE<-uoE*dt*(1+((E+L)/(K)))
  uL<-uoL*dt*(1+(Y*(E+L)/(K)))
  
  update(Be)<-rbinom(E,(dE+uE)*dt)
  update(Bl)<-rbinom(L, (dL+uL)*dt)
  update(Bp)<-rbinom(P,(dP+uP)*dt)
  update(Bm)<-rbinom(M,uM*dt)
  update(nt)<-rbinom(M,(dt/S))
  
  update(Reff)<-0.5*(Emax/(exp(uM*S)-1))*(1/(1+uE/dE))*(1/(1+uL/dL))*(1/(1+(uP*dt)/dP))
  
  update(E)<-if(E-Be>0)E-Be+rpois(nt*n) else rpois(nt*n)
  update(L)<-if(L-Bl>0)L-Bl+rbinom(Be,(dE/(uE+dE))) else rbinom(Be,(dE/(uE+dE)))
  update(P)<-if(P-Bp>0)P-Bp+rbinom(Bl,(dL/(uL+dL))) else rbinom(Bl,(dL/(uL+dL)))
  update(M)<-if(M-Bm>0)M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))-Bm else M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))
  
})




#########################################################################################################################################
#                                                                                                                                       #
#                                                     parameterisation with LHC's                                                       #
#                                                                                                                                       #
#########################################################################################################################################

##model simulation function - starts simulation from same random seed each time
modSim <- function(parms, obsDat) {
  modR <- larvalR(user=parms)
  simDat <- as.data.frame(modR$run(0:2000))
  simDat<-simDat[seq(1, NROW(simDat), by = 1/delta),]
  simDat$step<-simDat$step*delta
  matchedTimes <- simDat$step %in% garkiObs$time
  simDat$M[simDat$M<=0|is.na(simDat$M)] = 1e-3
  return(simDat)
}
##goodness of fit function - log liklihood
GOF<-function(pr){ #add parameters
  parameters<-read.table(text = pr, sep = ",", colClasses = "numeric")
  modSim <- modSim(parms=mosParamsP(uoE=parameters$V1,uoL=parameters$V2,uP=parameters$V3,
                                    sf=parameters$V4,Y=parameters$V5,n=parameters$V6,
                                    E0=177,L0=9,P0=1,M0=7),obsDat = garkiObs)
  matchedTimes <- modSim$step %in% garkiObs$time
  log_like <- sum(-dzipois(garkiObs$M,modSim$M[matchedTimes], log = TRUE))
  if (!is.na(log_like)){
    res<-c(parameters$V1,parameters$V2,parameters$V3,parameters$V4,parameters$V5,parameters$V6,log_like)
  }
  return(res)
}



HC<-NULL
##initialise latin hypercube
HCtemp<- randomLHS(5000, 6)
HC$uoE <- (HCtemp[,1])
HC$uoL <- (HCtemp[,2])
HC$uP <- (HCtemp[,3])
HC$Y <- (100*HCtemp[,4])
HC$n <- (25*HCtemp[,5])
HC$sf<-(25*HCtemp[,6])

HC<-as.data.frame(HC)


HC<- paste(HC$uoE,HC$uoL,HC$uP, HC$sf,HC$Y, HC$n, sep=",")#concatonate into single vector 

clusterExport(cl, c("rFx","delta","garkiObs","modSim"), envir=environment())


system.time(lhcResults<-data.frame(t(sapply(lapply(HC,GOF), `[`))))#convert results to dataframe

colnames(lhcResults)<-c("uoE","uoL","uP","Y","sf","n","logLike")

mins<-lhcResults[which(lhcResults$logLike == min(lhcResults$logLike)), ]#find minimum values



##run and plot model with set parameter values
set.seed(10) 
mod <- larvalR(user=mosParamsP(
                               uoE=mins$uoE,uoL=mins$uoL,uP=mins$uP,
                               Y=mins$Y,sf=mins$sf,n=mins$n, E0=177,L0=9,P0=1,M0=7)) #parameters estimated from LHC sampling
sim <- as.data.frame(mod$run(0:2000))

df<-sim[seq(1, NROW(sim), by = 1/delta),]
df$time<-df$step*delta

ggplot(data=df, aes(time))+ 
  geom_line(data=df,aes(x=time, y=M), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=garkiObs,(aes(x=time,y=M)),col="red")+
  theme_bw()


