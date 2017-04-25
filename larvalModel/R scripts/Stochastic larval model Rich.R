library(lubridate)
library(ggplot2)
library(VGAM)
library(odin)


setwd("C:\\ImperialMalaria\\larval model\\Data")

##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
rainfall$time<-(1:nrow(rainfall))

delta<-0.01##discrete time period
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
Si= 1, #number of sites
sf = 64, #scaling factor
dt=delta,
n=10
)
return(as.list(environment()))


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
  Si <- user()
  sf <- user()
  dt<-user()
  n<-user()
  rF[] <- user()
  dim(rF) <- user()
  trx<-(tr/dt)
  
  #initial values
  initial(Be)<-0
  initial(Bl)<-0
  initial(Bp)<-0

  initial(E) <- 0
  initial(L) <- 0
  initial(P) <- 0
  initial(M) <- 1
  initial(Bm)<-0
  
  initial(nt)<-0

  K <-if (step<=trx) (1+(sf*((1/trx)*(0)))) else (1+(sf*((1/trx)*(sum(rF[(step-trx):step])))))

  uE<-uoE*(1+(E+L)/(K/Si))
  uL<-uoL*(1+(Y*(E+L)/(K/Si)))
  
  
  update(Be)<-rbinom(E,(dE+uE)*dt)
  update(Bl)<-rbinom(L, (dL+uL)*dt)
  update(Bp)<-rbinom(P,(dP+uP)*dt)
  update(Bm)<-rbinom(M,uM*dt)
  
  update(nt)<-rbinom(M,(dt/S))
  
  update(E)<-if(E-Be>0)E-Be+rpois(nt*n) else rpois(nt*n) 
  update(L)<-if(L-Bl>0)L-Bl+rbinom(Be,(dE/(uE+dE))) else rbinom(Be,(dE/(uE+dE))) 
  update(P)<-if(P-Bp>0)P-Bp+rbinom(Bl,(dL/(uL+dL))) else rbinom(Bl,(dL/(uL+dL)))
  update(M)<-if(M-Bm>0)M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))-Bm else M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))
  
 
})


set.seed(10)
mod <- larvalR(user=mos_params(sf=111.9714,Y=6.793488)) #parameters estimated from LHC sampling
sim <- as.data.frame(mod$run(0:20000))

df = data.frame(time=sim$step, adMos=sim$M)
df2<-df[seq(1, NROW(df), by = 1/delta),]
df2$time<-df2$time*delta
ggplot(data=df2, aes(time))+ 
  #  geom_ribbon(aes(x=time, ymax=yH$M,ymin=yL$M), color='black',alpha=0.2) +
  geom_line(data=df2,aes(x=time, y=adMos), color='blue',lwd=1,alpha=0.5)+
#  geom_point(data=garkiObs,(aes(x=time,y=M)),col="red")+
  # geom_line(data=yR,aes(x=time, y=Reff), color='black',lwd=1)+
  # geom_line(data=rainfall,aes(x=date, y=rainfall), color='grey')+
  theme_bw()


