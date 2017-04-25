
##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
rainfall$time<-(1:nrow(rainfall))

larvMod <- function (Time, State, Pars) {
  

  with(as.list(c(State, Pars)), {

    K <- if (Time<=tr) (1+(sf*((1/tr)*(0))))  else (1+(sf*((1/tr)*(sum(rF[(Time-tr):Time])))))
    #print(E)
    
    #bb<-(rpois(1,lambda=(B*M)/Si)*Si)
    
    ## Derivatives
    dE <- (B*M)/Si-(E/e)-((uoE*((1+(E+L)/(K))))*E)
    dL <- (E/e)-(L/l)-((uoL*(1+(Y*(E+L)/(K))))*L)
    dP <- (L/l)-(P/p)-(uP*P)
    dM <- 0.5*(P/p)-(uM*M)
    
    uE<- (uoE*((1+(E+L)/(K))))
    uL <- (uoL*(1+(Y*(E+L)/(K))))
    

    # deriv(R0)<-0.5*(Emax/(exp(uM*S-1)))*(1/(1+uoE*dE))*(1/(1+uoL*dL))*(1/(1+uP*xP))-R0
   # Reff<-0.5*(Emax/(exp(uM*S)-1))*(1/(1+uE*xE))*(1/(1+uL*xL))*(1/(1+uP*xP))
    list(c(dE,dL,dP,dM))
  })
}

Pars <- c(mos_params2())
State <- c(E=0,L=0,P=0,M=2)
Time <- seq(1, 200, by = 1)

clusterExport(cl, c("State","Pars","Time","larvMod"), envir=environment())


out <- as.data.frame(ode(func = larvMod, y = State, parms = Pars, times = Time))

plot(out$M~out$time)



##Parameters from White et al (2011)
mos_params2 <- function(rF=rainfall$rainfall,
                       e = 6.64, #development time of early larval instars
                       l = 3.72, #development time of late larval instars
                       p = 0.64, #development time of pupae
                       uoE = 0.034, #per capita daily mortality rate of early instars (low density)
                       uoL = 0.035, #per capita daily mortality rate of late instars (low density)
                       uP = 0.25, #per capita daily mortality rate of pupae
                       uM = 0.096, # per capita daily mortality rate of adult An.gambiae
                       B = 21.19, #No. of eggs laid per day per mosquito
                       Y = 13.25, #effect of density dependence on late instars relative to early instars
                       S = 3, #duration of gonotrophic cycle
                       Emax = 93.6, #max number of eggs per oviposition per mosquito
                       tr = 4, #days of rainfall contributing to carrying capacity 
                       Si= 15 #number of sites
                       ,sf = 61.143070 #scaling factor
)
return(as.list(environment()))


g<-NULL
x<-1
B<-1
while(x<1000){
  B<-B+1
  M<-sample(100)
  g<-rbind(g,(rpois(1,lambda=(B*M)/Si)*Si))
  x<-x+1
  print(x)
  
}

