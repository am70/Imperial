mHyperC <-
function(sites){
HC<-NULL
##initialise latin hypercube
HCtemp<- randomLHS(1000, 9)
HC$dE <- 0.151# (HCtemp[,1])
HC$dL <- 0.269#(HCtemp[,2])
HC$dP <- 1.563#round((10*HCtemp[,3]),2)
HC$uoE <- (HCtemp[,4])
HC$uoL <- (HCtemp[,5])
HC$uP <- (HCtemp[,6])
HC$Y <- round((100*HCtemp[,7]),2)
#HC$tr <- (10*HCtemp[,3])
HC$sf <- round((50*HCtemp[,8]),2)
HC$n<-80#round((HCtemp[,9])+13,2)
HC2<-HC
HC<-as.data.frame(HC)


HC<- paste(HC$dE, HC$dL, HC$dP,HC$uoE,HC$uoL,HC$uP, HC$Y,HC$sf, HC$n, sep=",")#concatonate into single vector 
#system.time(lhcResults<-data.frame(t(sapply(lapply(HC,GOF), `[`))))#convert results to dataframe

system.time(lhcResults<-data.frame(t(sapply(snow::parLapply(cl,HC,likeFunc), `[`))))#convert results to dataframe

colnames(lhcResults)<-c("dE","dL","dP","uoE","uoL","uP","Y","sf","n","logLike")

mins<-lhcResults[which(lhcResults$logLike == min(lhcResults$logLike)), ]#find minimum values
plot(lhcResults$uoE,log(lhcResults$logLike))

modSim <- modSim(parms=mos_params(dE=mins$dE,dL=mins$dL,dP=mins$dP,
                                  uoE=mins$uoE,uoL=mins$uoL,uP=mins$uP,
                                  Y=mins$Y,sf=mins$sf,n=mins$n),obsDat = garkiObs)

par(mfrow=c(1,1))
ggplot(data=modSim, aes(step))+ 
  #  geom_ribbon(aes(x=time, ymax=yH$M,ymin=yL$M), color='black',alpha=0.2)+
  geom_line(data=modSim,aes(x=step, y=M), color='blue',lwd=1,alpha=0.5)+
  geom_point(data=garkiObs,(aes(x=time,y=M)),col="red")+
  # geom_line(data=yR,aes(x=time, y=Reff), color='black',lwd=1)+
  # geom_line(data=rainfall,aes(x=date, y=rainfall), color='grey')+
  theme_bw()


sens<-function(No){
  modSim <- modSim(parms=mos_params(dE=mins$dE,dL=mins$dL,dP=mins$dP,
                                    uoE=mins$uoE,uoL=mins$uoL,uP=mins$uP,
                                    Y=No,sf=mins$sf,n=mins$n),obsDat = garkiObs)
  matchedTimes <- modSim$step %in% garkiObs$time
  log_like <-max(modSim$Reff)#sum(-dzipois(garkiObs$M,modSim$M[matchedTimes], log = TRUE))#
  log_like<-cbind(log_like,No)

return(log_like)
  
}



clusterExport(cl, c("mins"), envir=environment())

like<-as.data.frame(t(sapply(snow::clusterApply(cl,c(1:200),sens), `[`)))
#like<-as.data.frame(t(sapply(lapply(c(1:100),sens), `[`)))


plot(like$V2,like$V1)


prof.lower = like$V1[1:which.min(like$V1)]
prof.avec = like$V2[1:which.min(like$V1)]

a<-approx(prof.lower, prof.avec, xout = mins$logLike+
         qchisq(0.95, 1)/2)



prof.upper = like$V1[nrow(like):which.min(like$V1)]
prof.avecU = like$V2[nrow(like):which.min(like$V1)]

a1<-approx(prof.upper, prof.avecU, xout = mins$logLike +
         qchisq(0.95, 1)/2)


ss1<-smooth.spline(prof.avec,prof.lower)
ss1$x[ss$y[mins$logLike + qchisq(0.95, 1)/2]]

ss2<-smooth.spline(prof.avecU, prof.upper)
ss2$x[ss$y[ mins$logLike + qchisq(0.95, 1)/2]]


##plot results
#par(mfrow=c(1,2))
#plot(lhcResults$X1,log(lhcResults$X3))
#points(mins$X1,log(mins$X3),col="red",pch=16)
#plot(lhcResults$X2,log(lhcResults$X3))
#points(mins$X2,log(mins$X3),col="red",pch=16)
##run and plot model with set parameter values
#set.seed(10)
mod <- larvalMod(user=mos_params(sf=mins$X1,Y=mins$X2,O=mins$X3,n=sites)) #parameters estimated from LHC sampling
sim <- as.data.frame(mod$run(0:20000))
#plot(sim$M,col="white")
plot(sim$M,col="blue")
##run logistic regression to find max growth rate
R<-logisticReg(mos_params(sf=mins$X1,Y=mins$X2,n=sites))
print(R)
return(R)
}
