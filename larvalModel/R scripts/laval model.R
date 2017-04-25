#install_github("richfitz/odin")
library(odin)
library(stringr)
library(ggplot2)
library(lattice)
library(growthrates)
setwd("C:\\ImperialMalaria\\larval model\\Data")

##load in Garki rainfall data
rainfall<-read.csv("garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-(1:nrow(rainfall))
#limit data to one year
rainfall<-subset(rainfall, date >730 & date<1104)

##white et al larval model
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
  
  #initial values
  initial(Be)<-0
  initial(Bl)<-0
  initial(Bp)<-0
  
  initial(E) <- 177
  initial(L) <- 8
  initial(P) <- 1
  initial(M) <- 7
  initial(Bm)<-0
  initial(nt)<-0
  initial(Reff)<-0
  
  K <-if (step<=trx) (1+(sf*((1/trx)*(0)))) else (1+(sf*((1/trx)*(sum(rF[(step-trx):step])))))
  
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
  
  
  
  config(base) <- "larvalModelR"
})


##Parameters from White et al (2011)
pars <- function(rF=rainfall$rainfall,
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
                 sf =  64 #scaling factor calculated using mcmc for current garki data
)
return(as.list(environment()))



##fit and plot model##
t <- seq(0, 200)


modR <- larvalR(user=pars())
yR <- modR$run(t)
yR<-as.data.frame(yR)


par(mfrow=c(1,1))
ggplot(data=yR, aes(time)) + 
  #  geom_ribbon(aes(x=time, ymax=yH$M,ymin=yL$M), color='black',alpha=0.2) +
  geom_line(data=yR,aes(x=time, y=M), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=garkiObs,(aes(x=time,y=M)),col="red")+
  # geom_line(data=yR,aes(x=time, y=Reff), color='black',lwd=1)+
  # geom_line(data=rainfall,aes(x=date, y=rainfall), color='grey')+
  theme_bw()


##calculate change in Reff and plot changes

runMod<-function(sites){
  modR <- larvalR(user=pars(Si=sites))
  yR <- modR$run(t)
  yR<-as.data.frame(yR)
  max(yR$Reff)
}

yReff<-as.data.frame(sapply((1:50),FUN=runMod))
yReff$si<-c(1:nrow(yReff))
colnames(yReff)<-c("Reff","si")

ggplot(data=yReff, aes(si)) + 
  geom_line(data=yReff,aes(x=si, y=Reff), color='blue',lwd=1,alpha=0.5) +
  geom_point(data=yReff,(aes(x=si,y=Reff)),col="red")+
  xlab("N")+ylab("Effective reproduction number (Reff)")+ggtitle("Methods from White et al for effective reproduction number")+
  theme_bw()




##########################################################Logistic growth models############################################################
df = data.frame(time=yR$time, adMos=yR[,5])
df2<-df[c(1:90),]##data for parametric fit

#fit parametric non least squares model
population=df2$adMos
time=df2$time
#guesstimate intial starting parameters
SS<-getInitial(adMos~SSlogis(time,alpha,xmid,scale),data=df2)
K_start<-SS["alpha"]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
#fit the nonlinear least squares regression
m<-nls(population~K*N0*exp(R*time)/(K+N0*(exp(R*time)-N0)),start=list(K=K_start,R=R_start,N0=N0_start))
#estimated parameters
summary(m)
res<-as.data.frame(cbind(time,predict(m)))

ggplot(data=res, aes(V2)) + 
  geom_line(data=res,aes(x=time, y=V2), color='blue',lwd=1,alpha=0.5) +
  geom_point((aes(x=time,y=population)),col="red",alpha=0.5)+
  xlab("Time")+ylab("population")+ggtitle("logistic growth model")+
  theme_bw()


logGrowthMod <- odin::odin({
  ##Parameters
  K <- user()
  r <- user()
  No <- user()
  
  ## Derivatives
  deriv(N) <- (K*No*exp(r*time))/(K+No*(exp(r*time)-No))
  deriv(time)<- +1
  
  ## Initial conditions
  initial(N) <- 1
  initial(time) <-1
  
  config(base) <- "logGrowthMod"
})

parsG <- list(
  K=(coef(m)[1]), #carrying capacity
  r=(coef(m)[2]), #growth rate
  No=(coef(m)[3]) #starting pop
)

t <- seq(0, 300)
modG <- logGrowthMod(user=parsG)
yG <- modG$run(t)
yG<-as.data.frame(yG)
yG[,3]<-c(1:301)

#ggplot(col=Legend,data=yR, aes(time)) +
 # ylim(0,max(garkiObs$M))+ 
  # geom_ribbon(aes(x=time, ymax=yH$M,ymin=yL$M), color='black',alpha=0.2)+
 # geom_point(data=garkiObs,(aes(x=time,y=M,col="Spray Catches")))+
 # geom_line(data=yR,aes(x=time, y=M, col='Model'),lwd=1,alpha=0.5) +
 # scale_color_manual(values=c("Spray Catches"="red", "Model"="blue"))+
 # geom_line(aes(x=yG[,3], y=yG[,2]), color='black',lty="dashed",lwd=1,alpha=0.5)+
 # theme_bw()+ theme(legend.title = element_blank())+ylab("Adult Mosquitoes")



##put growth models into function to run for different numbers of sites

growthMods<-function(sites){
  
  modR <- larvalR(user=pars(Si=sites))
  yR <- modR$run(t)
  yR<-as.data.frame(yR)
  
  df = data.frame(time=t, adMos=yR[,5])
  df2<-df[c(1:90),]##data for parametric fit
  
  #fit parametric non least squares model
  population=df2$adMos
  time=df2$time
  #guesstimate intial starting parameters
  SS<-getInitial(adMos~SSlogis(time,alpha,xmid,scale),data=df2)
  K_start<-SS["alpha"]
  R_start<-1/SS["scale"]
  N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
  #fit the nonlinear least squares model
  m<-nls(population~K*N0*exp(R*time)/(K+N0*(exp(R*time)-N0)),start=list(K=K_start,R=R_start,N0=N0_start))
  
 # plot(population~time)
 #  lines(predict(m)~time)
  
  parsG <- list(
    K=(coef(m)[1]), #carrying capacity
    r=(coef(m)[2]), #growth rate
    No=(coef(m)[3]) #starting pop
  )
  
  #modG <- logGrowthMod(user=parsG)
  #yG <- modG$run(t)
 #yG<-as.data.frame(yG)
  #yG[22,]$N
   parsG$r
}


##run growth model for differing numbers of sites
yGr<-as.data.frame(sapply((1:40),FUN=growthMods))
yGr$time<-c(1:nrow(yGr))
colnames(yGr)<-c("growth","si")

##plot growth model~sites results
ggplot(data=yGr, aes(si))+ 
  geom_line(data=yGr,aes(x=si, y=growth), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=yGr,(aes(x=si,y=growth)),col="red")+
  xlab("N")+ylab("r")+ggtitle("Effect of N on r in logistic growth model")+
  theme_bw()

