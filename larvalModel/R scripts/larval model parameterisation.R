library(FME)
library(lubridate)
setwd("C:\\ImperialMalaria\\larval model\\Data")

##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
rainfall$time<-(1:nrow(rainfall))

#load in mosquito data
garki<-read.table("spraycollect.txt",sep=",",head=T)
garki$ag_sum<-rowSums(garki[,c(8:16)])#create sum of A.gambiae spray samples
garki$date<-gsub(" 00:00:00", "", garki$date)
garki$date<-dmy(garki$date)
garki72<-subset(garki,date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
garki72_101<-subset(garki72,e_station==101) ##subset by village - needs checking with michael specific villages
garki72_101<-merge(garki72_101,rainfall,by.x="date",by.y="date",all=T) #merge rainfall and mosquito data



##Parameters from White et al (2011)
parsX <- list(
              dE = 6.64, #development time of early larval instars
              dL = 3.72, #development time of late larval instars
              dP = 0.64, #development time of pupae
              uoE = 0.034, #per capita daily mortality rate of early instars (low density)
              uoL = 0.035, #per capita daily mortality rate of late instars (low density)
              uP = 0.25, #per capita daily mortality rate of pupae
              uM = 0.096, # per capita daily mortality rate of adult An.gambiae
              B = 21.19, #No. of eggs laid per day per mosquito
              Y = 13.25, #effect of density dependence on late instars relative to early instars
              S = 3, #duration of gonotrophic cycle
              Emax = 93.6, #max number of eggs per oviposition per mosquito
              tr = 4, #days of rainfall contributing to carrying capacity 
              Si= 1, #number of sites
              sf = 50 #scaling factor
)

##Parameters from White et al (2011) upper 95%
parsH <- list(rF=rainfall$rainfall,
              dE = 4.82, #development time of early larval instars
              dL = 2.03, #development time of late larval instars
              dP = 0.07, #development time of pupae
              uoE = 0.024, #per capita daily mortality rate of early instars (low density)
              uoL = 0.025, #per capita daily mortality rate of late instars (low density)
              uP = 0.18, #per capita daily mortality rate of pupae
              uM = 0.087, # per capita daily mortality rate of adult An.gambiae
              B = 25.31, #No. of eggs laid per day per mosquito
              Y = 9.82, #effect of density dependence on late instars relative to early instars
              S = 3, #duration of gonotrophic cycle
              Emax = 93.6, #max number of eggs per oviposition per mosquito
              tr = 7, #days of rainfall contributing to carrying capacity 
              Si= 1, #number of sites
              sf = 66.97094 #scaling factor
)

##Parameters from White et al (2011) lower 95%
parsL <- list(rF=rainfall$rainfall,
              dE = 8.53, #development time of early larval instars
              dL = 5.61, #development time of late larval instars
              dP = 1.47, #development time of pupae
              uoE = 0.044, #per capita daily mortality rate of early instars (low density)
              uoL = 0.044, #per capita daily mortality rate of late instars (low density)
              uP = 0.32, #per capita daily mortality rate of pupae
              uM = 0.1, # per capita daily mortality rate of adult An.gambiae
              B = 11.57, #No. of eggs laid per day per mosquito
              Y = 17.51, #effect of density dependence on late instars relative to early instars
              S = 3, #duration of gonotrophic cycle
              Emax = 93.6, #max number of eggs per oviposition per mosquito
              tr = 3, #days of rainfall contributing to carrying capacity 
              Si= 1, #number of sites
              sf = 66.97094 #scaling factor
)



##white et al larval model
larvalRs <- odin::odin({
  ##Parameters
  dE <- user()
  dL <- user()
  dP <- user()
  uoE <- user()
  uoL <- user()
  uP <- user()
  uM <- user()
  B <- user()
  Y <- user()
  S <- user()
  Emax <- user()
  tr <- user()
  Si <- user()
  sf <- user()
  rF[] <- user()
  
  dim(rF) <- user()
  
  ## Derivatives
  deriv(E) <- (B*M)-(E*dE)-((uoE*((1+(E+L)/(K/Si))))*E)
  deriv(L) <- (E*dE)-(L*dL)-((uoL*(1+(Y*(E+L)/(K/Si))))*L)
  deriv(P) <- (L*dL)-(P*dP)-(uP*P)
  deriv(M) <- 0.5*(P/dP)-(uM*M)
  
  uE<- (uoE*((1+(E+L)/(K/Si))))
  uL <- (uoL*(1+(Y*(E+L)/(K/Si))))
  
  K <-if (step<=trx) (1+(sf*((1/trx)))) else (1+(sf*((1/trx)*(sum(rF[(step-trx):step])))))
  deriv(time) <- +1
  
  
  # deriv(R0)<-0.5*(Emax/(exp(uM*S-1)))*(1/(1+uoE*E))*(1/(1+uoL*L))*(1/(1+uP*dP))-R0
  deriv(Reff)<-0.5*(Emax/(exp(uM*S)-1))*(1/(1+uE*dE))*(1/(1+uL*dL))*(1/(1+uP*dP))-Reff
  
  ## Initial conditions
  initial(E) <- 0
  initial(L) <- 0
  initial(P) <- 0
  initial(M) <- 1
  initial(time) <-0
  #initial(R0)<-0
  initial(Reff)<-0
  
  config(base) <- "larvalModelR"
})


Objective <- function (x) {         
  parsX[] <- x
  garkiObs<-garki72_101[,c(39,37)]
  colnames(garkiObs)<-c("time","M")
  garkiObs<-subset(garkiObs,M>=0)
  
  t <- seq(0, 150)
  modR <- larvalRs(user=parsX)
  yR <- modR$run(t)
  yR<-as.data.frame(yR)
  yRx<-yR[,c(1,5)]
  colnames(yRx)<-c("time","M")
  return(modCost(obs = garkiObs, model = yRx))
}



Coll <- collin(sF <- sensFun(func = Objective, parms = parsX))
plot(Coll, log = "y")
abline(h = 20, col = "red")


print(system.time(Fit <- modFit(p = c( sf = 50 #development time of early larval instars
                                       ), f = Objective, lower = c(0.0))))




init <- larvalRs(user=parsX)
init<- modR$run(t)
parsX[c("sf")] <- Fit$par
modR <- larvalRs(user=parsX)
yR <- modR$run(t)
yR<-as.data.frame(yR)
yRx<-yR[,c(1,5)]
colnames(yRx)<-c("time","M")
Cost <- modCost(obs = garkiObs, model = yRx)
Cost

plot(yRx, init, xlab = "time, hour", ylab = "molC/m3", lwd = 2,
      obs = garkiObs, obspar = list(cex = 2, pch = 18))


