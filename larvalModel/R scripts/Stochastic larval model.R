delta<-0.01##discrete time period
rF<-as.data.frame(rainfall$rainfall)
rFx<-rF[rep(seq_len(nrow(rF)), each=1/delta),] #split rainfall data into discreet time periods

mos_params <- function(
#parameter estimates from white et al (2011)
rF=rFx,
dE = 1/6.64, #development time of early larval instars
dL = 1/3.72, #development time of late larval instars
dP = 1/0.64, #development time of pupae
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
sf = 8, #scaling factor
dt=delta,
n=20
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
  

  update(E)<-if(E-Be>0)(E-Be+rpois(nt*n)) else rpois(nt*n)
  update(L)<-if(L-Bl>0)L-Bl+rbinom(Be,(dE/(uE+dE))) else rbinom(Be,(dE/(uE+dE)))
  update(P)<-if(P-Bp>0)P-Bp+rbinom(Bl,(dL/(uL+dL))) else rbinom(Bl,(dL/(uL+dL)))
  update(M)<-if(M-Bm>0)M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))-Bm else M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))
  

 
})


set.seed(10)
mod <- larvalR(user=mos_params(Si=1,sf=252.0122,Y=11.39073,n=10))
yy <- as.data.frame(mod$run(0:20000))

df = data.frame(time=yy$step, adMos=yy$M)
df2<-df[seq(1, NROW(df), by = 1/delta),]
df2$time<-df2$time*delta
ggplot(data=df2, aes(time))+ 
  #  geom_ribbon(aes(x=time, ymax=yH$M,ymin=yL$M), color='black',alpha=0.2) +
  geom_line(data=df2,aes(x=time, y=adMos), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=garkiObs,(aes(x=time,y=M)),col="red")+
  # geom_line(data=yR,aes(x=time, y=Reff), color='black',lwd=1)+
  # geom_line(data=rainfall,aes(x=date, y=rainfall), color='grey')+
  theme_bw()

#sum(-dzipois(garkiObs$M,df2$adMos[matchedTimes], log = TRUE))




###growth models


##put growth models into function to run for different numbers of sites

growthMods<-function(sites){
  set.seed(1)
  repeat{
  mod <- larvalR(user=mos_params(Si=sites,sf=8, n=1))
  yy <- as.data.frame(mod$run(0:100000))
  df = data.frame(time=yy$step, adMos=yy$M)
  df2<-df[seq(1, NROW(df), by = 1/delta),]
  df2$time<-df2$time*delta
  df2<-df2[c(20:80),]##data for parametric fit
  
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
  
   plot(population~time)
    lines(predict(m)~time)
  
  parsG <- list(
    K=(coef(m)[1]), #carrying capacity
    r=(coef(m)[2]), #growth rate
    No=(coef(m)[3]) #starting pop
  )
 break(sum(is.na(df2$adMos))==0) 
}
  #t <- seq(0, 300)
#  modG <- logGrowthMod(user=parsG)
 # yG <- modG$run(t)
 # yG<-as.data.frame(yG)
 # yG[,3]<-c(1:301)
  
  #yG[22,2]

  parsG$r
}


##run growth model for differing numbers of sites
yGr<-as.data.frame(sapply((1:2),FUN=growthMods))
yGr$time<-c(1:nrow(yGr))
colnames(yGr)<-c("growth","si")

##plot growth model~sites results
ggplot(data=yGr, aes(si))+ 
  geom_line(data=yGr,aes(x=si, y=growth), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=yGr,(aes(x=si,y=growth)),col="red")+
  xlab("N")+ylab("r")+ggtitle("variation of site number")+
  theme_bw()




