
#########################################################################################################################################
#                                                                                                                                       #
#                                                     visualise results and extract C.I.                                                #
#                                                                                                                                       #
#########################################################################################################################################
library(miscTools)

cIfunc<-function(x){
  quantRes<-as.data.frame(c(1:2))
  for(i in 1:(ncol(x)-1)){
    resD<-density(x[,i],adjust=1)################CHECK DENSITY ADJUSTMENTS
    samp <- sample(resD$x, 1e6, replace = TRUE, prob = resD$y)
    quants<-as.data.frame(quantile(samp, c(0.05, 0.95)))
    colnames(quants)<-c(names(x[i]))
    quantRes<-cbind(quantRes,quants)
  }
  quantRes[,-1]
}

iqrFunc<-function(x){
  iqrRes<-as.data.frame(c(1:2))
  
  for(i in 1:(ncol(x)-1)){
    resD<-density(x[,i])
    resD<-IQR(x[,i])
    l<-median(x[,i])-resD
    h<-median(x[,i])+resD
    h<-rbind(l,h)
    iqrRes<-cbind(iqrRes,h)
  }
  iqrRes[,-1]
}

<<<<<<< Updated upstream
plotMod<-function(n,mcmcRes,obsDat,rf,sf,rFx,z){
=======


plotMod<-function(n,mcmcRes,obsDat,rf,sf,rFxs,z){
>>>>>>> Stashed changes
  parms<-colMedians(mcmcRes)
  #cI<-cIfunc(mcmcRes)
  #iqr<-iqrFunc(mcmcRes)
  startTime<-min(obsDat$time)/delta
  endTime<-max(obsDat$time)/delta
<<<<<<< Updated upstream
  
  
  #  crosscorr.plot(g$results[,-13])
  
  ####95% CI####
  par(mfrow=c(1,1))
  mid<- pFilt(2,iState,modStep3,dataLikFunc,obsDat,pr=c(parms[1:6],10^parms[z],10^parms[sf],parms[15:22]),rFclust=3,resM=T,rFx=rFx3,cluster=F)
  colnames(mid)<-c("E","L","P","M")
  mid <- mid$M[seq(1, length(mid$M), 4)]
  plot(mid,col="white")
  lines(mid,col="blue")
  points(obsDat$`408`/parms[22]~obsDat$time,col="red")
  
  
  #for(i in 1:nrow(g)){
  #fs<-betaBinom(g[i, 3], g[i,2], 0.01, parms[6])
  
  # print(fs)
  # }
  
  
  low<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(cI[1,1],cI[1,2],cI[1,3],
                                                        cI[1,4],cI[2,z],cI[2,9],cI[2,sf],cI[1,14]),rFclust=rf,resM=T,fxdParams=fxd)
  colnames(low)<-c("E","L","P","M")
  
  high<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(cI[2,1],cI[2,2],cI[2,3],
                                                         cI[2,4],cI[2,z],cI[1,9],cI[1,sf],cI[2,14]),rFclust=rf,resM=T,fxdParams=fxd)
  colnames(high)<-c("E","L","P","M")
  
  ####IQR####
  
  lowIqr<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(iqr[1,1],iqr[1,2],iqr[1,3],
                                                           iqr[1,4],iqr[2,5],iqr[2,9],iqr[2,sf],iqr[1,14]),rFclust=rf,resM=T,fxdParams=fxd)
  colnames(lowIqr)<-c("E","L","P","M")
  
  highIqr<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(iqr[2,1],iqr[2,2],iqr[2,3],
                                                            iqr[2,4],cI[2,z],iqr[1,9],iqr[1,sf],iqr[2,14]),rFclust=rf,resM=T,fxdParams=fxd)
  colnames(highIqr)<-c("E","L","P","M")
  
  mid$time<-1:nrow(mid)*0.25
  
  timeX<-as.data.frame(1:nrow(mid))*delta
  colnames(timeX)<-c("time")
  obsDat$time<-obsDat$time-min(obsDat$time)
  timeX<-merge(obsDat,timeX,all=T)
  colnames(timeX)<-c("time","M")
  timeX<-timeX[-1,]
  #timeX[is.na(timeX)]<-0
  print(as.numeric(max(timeX$M,na.rm=T))/0.01)
  
  ggplot(data = mid,aes(x=mid$time,y=timeX$M/0.01))+
    expand_limits(y=c(0,(as.numeric(max(timeX$M,na.rm=T))/0.01)))+
    # geom_ribbon(aes(x=mid$time, ymax=highIqr$M, ymin=lowIqr$M), fill="dark grey", alpha=.5)+
    #geom_ribbon(aes(x=mid$time, ymax=high$M, ymin=low$M), fill="grey", alpha=.5)+
    #geom_line(aes(x=mid$time,y = low$M), colour = 'dark grey')+
    # geom_line(aes(x=mid$time,y = high$M), colour = 'dark grey')+
    # geom_line(aes(x=mid$time,y = lowIqr$M), colour = 'dark grey')+
    # geom_line(aes(x=mid$time,y = highIqr$M), colour = 'dark grey')+
    geom_point(x=timeX$time,y=timeX$M/0.01,col="red")+
    geom_line(aes(x=mid$time,mid$M))+
    xlab("time")+
    ylab("M")+
    theme_bw()
  
}

hh1<-plotMod(44,cTest,garkiObs154,1,10,4,5)
=======

####95% CI####

  par(mfrow=c(1,1))
  mid<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(parms[1:6],10^parms[z],10^parms[sf],parms[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F)
  colnames(mid)<-c("E","L","P","M","rEff")
  plot(mid$rEff,col="white")
  lines(mid$rEff)
  mid <- mid$M[seq(1, length(mid$M), 4)]
 plot(mid*parms[22],col="white")
  lines(mid*parms[22],col="blue")
  points(obsDat$`408`~obsDat$time,col="red")

 
 #low<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(cI[2,1],cI[2,2],cI[2,3],
 #                                                      cI[2,4],cI[1,5],cI[1,6],10^cI[1,z],10^cI[1,sf],cI[2,15],cI[2,16],cI[2,17],cI[1,18],cI[1,19],cI[2,20],cI[2,21],cI[1,22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F)
# colnames(low)<-c("E","L","P","M")
# low <- low$M[seq(1, length(low$M), 4)]
 #lines(low*cI[2,22])
# lines(mid*parms[22],col="blue")
 
 
# high<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(cI[2,1],cI[2,2],cI[2,3],
#                                                        cI[2,4],cI[2,5],cI[2,6],10^cI[1,z],10^cI[1,sf],cI[1,15],cI[1,16],cI[1,17],cI[2,18],cI[2,19],cI[1,20],cI[2,21],cI[2,22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F)
 #colnames(high)<-c("E","L","P","M")
 #high <- high$M[seq(1, length(high$M), 4)]
 #plot(high*cI[2,22])
 
 
 ####IQR####
 
 #lowIqr<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(iqr[1,1],iqr[1,2],iqr[1,3],
  #                                                             iqr[1,4],iqr[2,5],iqr[2,9],iqr[2,sf],iqr[1,14]),rFclust=rf,resM=T,fxdParams=fxd)
 #colnames(lowIqr)<-c("E","L","P","M")
 
# highIqr<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(iqr[2,1],iqr[2,2],iqr[2,3],
 #                                                               iqr[2,4],cI[2,z],iqr[1,9],iqr[1,sf],iqr[2,14]),rFclust=rf,resM=T,fxdParams=fxd)
 #colnames(highIqr)<-c("E","L","P","M")
 
#mid$time<-1:nrow(mid)*0.25


  
 # ggplot(data = mid,aes(x=mid$time,y=timeX$M/0.01))+
  # expand_limits(y=c(0,(as.numeric(max(timeX$M,na.rm=T))/0.01)))+
   # geom_ribbon(aes(x=mid$time, ymax=highIqr$M, ymin=lowIqr$M), fill="dark grey", alpha=.5)+
    #geom_ribbon(aes(x=mid$time, ymax=high$M, ymin=low$M), fill="grey", alpha=.5)+
    #geom_line(aes(x=mid$time,y = low$M), colour = 'dark grey')+
   # geom_line(aes(x=mid$time,y = high$M), colour = 'dark grey')+
   # geom_line(aes(x=mid$time,y = lowIqr$M), colour = 'dark grey')+
   # geom_line(aes(x=mid$time,y = highIqr$M), colour = 'dark grey')+
   # geom_point(x=timeX$time,y=timeX$M/0.01,col="red")+
    #geom_line(aes(x=mid$time,mid$M))+
  #  xlab("time")+
   # ylab("M")+
    #theme_bw()
  
}

hh1<-plotMod(44,cTest,garkiObs408,3,11,rFx3,7)
>>>>>>> Stashed changes
hh2<-plotMod(44,pp$results,garkiObs104,1,11,7,6)
hh3<-plotMod(44,pp$results,garkiObs219,2,12,7,7)
hh4<-plotMod(44,pp$results,garkiObs220,2,13,7,8)

grid.arrange(hh1,hh2,hh3,hh4)


<<<<<<< Updated upstream
simX<-c(1:964)
for (i in 1:100){
  
  mod <- odinPackage::larvalModP(user=mosParamsP(uoE=parms[1],uoL=parms[2],uP=parms[3],
                                                 Y=parms[4],n=10,sf=parms[13],o=parms[9],tr=14,time1=1244,E0=f[1],L0=f[2],P0=f[3],M0=f[4]))
  sim <- as.data.frame(mod$run(1:964))
  simX<-cbind(simX,sim$M)
  
}
# 101  ,   104   ,  108  ,   113
simMean<-rowMeans(simX[,-1])
plot(simMean~sim$timeX,ylim=c(0,50),col="white")
lines(simMean~sim$timeX)
points(garkiObs220$time/delta,(garkiObs220$`220`/parms[8]),col="red")
=======

>>>>>>> Stashed changes
#########################################################################################################################################
#                                                                                                                                       #
#                                                     density plots                                                                     #
#                                                                                                                                       #
#########################################################################################################################################

library("gridExtra")
library("grid")
library("miscTools")

plotDens<-function(mcmc,mcmc2,title){
  uoEprior<-rnorm( nrow(mcmc),mean=0.035,sd=0.0056)
  uoLprior<-rnorm( nrow(mcmc),mean=0.035,sd=0.0056)
  uPprior<-rnorm( nrow(mcmc),mean=0.25,sd=0.05)
  Yprior<-rnorm( nrow(mcmc),mean=13.06,sd=3)
  
dEprior<-rnorm( nrow(mcmc),mean=0.150602,sd=0.04)
dLprior<-rnorm( nrow(mcmc),mean=0.268812,sd=0.06)
dPprior<-rnorm( nrow(mcmc),mean=1,sd=0.1)

<<<<<<< Updated upstream
plotDens<-function(mcmc){
  uoEprior<-rnorm( nrow(mcmc),mean=0.035,sd=0.009)
  uoLprior<-rnorm( nrow(mcmc),mean=0.035,sd=0.009)
  uPprior<-rnorm( nrow(mcmc),mean=0.25,sd=0.0557)
  Yprior<-rnorm( nrow(mcmc),mean=13.06,sd=2)
  dEprior<-rnorm( nrow(mcmc),mean=0.150602,sd=0.01)
  dLprior<-rnorm( nrow(mcmc),mean=0.268812,sd=0.055)
  dPprior<-rnorm( nrow(mcmc),mean=1.563,sd=0.05)
  #p0prior<-rnorm(94999,mean=0.5,sd=0.015)
=======
tauprior<-rnorm( nrow(mcmc),mean=7,sd=1.5)
uMprior<-rnorm( nrow(mcmc),mean=0.091,sd=0.005)
    #p0prior<-rnorm(94999,mean=0.5,sd=0.015)
>>>>>>> Stashed changes
  
  
  p1<-ggplot(data = mcmc,aes(uoE))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(uoEprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p2<-ggplot(data = mcmc,aes(uoL))+
    geom_density(kernel = "gaussian", adjust = 3, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(uoLprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p3<-ggplot(data = mcmc,aes(uP))+
    geom_density(kernel = "gaussian", adjust = 5, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(uPprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p4<-ggplot(data = mcmc,aes(Y))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(Yprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p5<-ggplot(data = mcmc,aes(dE))+
<<<<<<< Updated upstream
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(dEprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
=======
   geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(dEprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
   scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
>>>>>>> Stashed changes
    theme_bw()
  p6<-ggplot(data = mcmc,aes(dL))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(dLprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p7<-ggplot(data = mcmc,aes(dP))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(dPprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
<<<<<<< Updated upstream
  
  
  grid.arrange(p1, p2, p3, p4,p5,p6,p7)
}





#########################################################################################################################################
#                                                                                                                                       #
#                                                     likelihood test                                                                   #
#                                                                                                                                       #
#########################################################################################################################################


lk<-function(fitParams,f){
  ll<-0
  pX<-NULL
  for (i in 1:4){
    globalParms<-fitParams[c(1:10)]
    globalParms[10]<-fitParams[c(i+9)]#i+6 to fit the scaling factor specific for each village
    garkDat<-garkiObsX[,c(1,i+1)]#i+1 as first column is "time"
    
    simX<-c(1:2440)
    for (j in 1:100){
      
      mod <- odinPackage::larvalModP(user=mosParamsP(uoE=globalParms[1],uoL=globalParms[2],uP=globalParms[3],
                                                     Y=globalParms[4],n=globalParms[5],sf=globalParms[10],o=globalParms[9],tr=14,time1=3110,E0=f[1],L0=f[2],P0=f[3],M0=f[4]))
      sim <- as.data.frame(mod$run(1:2440))
      simX<-cbind(simX,sim$M)
      
    }
    # 101  ,   104   ,  108  ,   113
    simMean<-rowMeans(simX[,-1])
    simMean
    simMeanT<-cbind(simMean,sim$timeX/10)
    colnames(simMeanT)<-c("M","time")
    
    lDat<-merge(garkiObsX,simMeanT,by.x="time",by.y="time")
    print(lDat[i+1])
    
    ll<-ll+sum(dnbinom(as.numeric(unlist(lDat[i+1])),lDat$M,p=globalParms[6],log=T))+lprior(globalParms)
    
=======
  p8<-ggplot(data = mcmc,aes(tau))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(tauprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p9<-ggplot(data = mcmc,aes(uM))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(uMprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()

  p10<-ggplot(data = mcmc,aes(w))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p11<-ggplot(data = mcmc,aes(n))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p12<-ggplot(data = mcmc,aes(z1))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p13<-ggplot(data = mcmc,aes(z4))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p14<-ggplot(data = mcmc,aes(z5))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p15<-ggplot(data = mcmc,aes(z6))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p12<-ggplot(data = mcmc,aes(sf1))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p13<-ggplot(data = mcmc,aes(sf4))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p14<-ggplot(data = mcmc,aes(sf5))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p15<-ggplot(data = mcmc,aes(sf6))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p16<-ggplot(data = mcmc,aes(Mg))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  p17<-ggplot(data = mcmc,aes(p))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
    theme_bw()
  if(title =="power"){
    p18<-ggplot(data = mcmc,aes(o))+
      geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
      scale_colour_manual(values=c("Posterior"="red"), name="Densities")+
      theme_bw()
    grid.arrange(p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,main=textGrob(title,gp=gpar(fontsize=20,font=3)))
                 
>>>>>>> Stashed changes
  }
else
  grid.arrange(p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,main=textGrob(title,gp=gpar(fontsize=20,font=3))
)
}


plotDensLL<-function(mcmc,mcmc2){
  ggplot(data = mcmc,aes(ll))+
    geom_density(kernel = "gaussian", adjust = 1, aes(colour="Clumped"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(mcmc2$ll,colour="Non-Clumped"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Clumped"="red", "Non-Clumped"="blue"), name="Densities")+
    geom_vline(xintercept = median(mcmc$ll),lty="dashed",col="red")+
    geom_vline(xintercept = median(mcmc2$ll),lty="dashed",col="blue")+
    theme_bw()
}




##########plot and record

<<<<<<< Updated upstream
plotRec<-function(clustRes,name,tr){
=======
plotRec<-function(clustRes,clustRes2,name,tr){
>>>>>>> Stashed changes
  res<-clustRes
  w=1800
  h=1200
  
  fileName<-paste0("C:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",name)
  dir.create(fileName)
  
  #write.table(res,paste0(fileName,"\\mcmcRes.csv"))

  
  png(filename=paste0(fileName,"\\density.png"),width = w, height = h)
  plotDens(as.data.frame(res),as.data.frame(clustRes2),name)
  dev.off()


  png(filename=paste0(fileName,"\\fig1.png"),width = w, height = h,pointsize=30)
  plot(res[,c(1:3)])
  dev.off()
  png(filename=paste0(fileName,"\\fig2.png"),width = w, height = h,pointsize=30)
  plot(res[,c(4:6)])
  dev.off()
  png(filename=paste0(fileName,"\\fig3.png"),width = w, height = h,pointsize=30)
  plot(res[,c(7:9)])
  dev.off()
  png(filename=paste0(fileName,"\\fig4.png"),width = w, height = h,pointsize=30)
  plot(res[,c(10:12)])
  dev.off()
  png(filename=paste0(fileName,"\\fig5.png"),width = w, height = h,pointsize=30)
  plot(res[,c(13:15)])
  dev.off()
<<<<<<< Updated upstream
  png(filename=paste0(fileName,"\\fig5.png"),width = w, height = h,pointsize=30)
  plot(res[,c(15:18)])
  dev.off()
  png(filename=paste0(fileName,"\\fig6.png"),width = w, height = h,pointsize=30)
  plot(res[,c(19:21)])
  dev.off()
  png(filename=paste0(fileName,"\\fig7.png"),width = w, height = h,pointsize=30)
  plot(res[,c(22:23)])
  dev.off()
  
  hh1<-plotMod(44,res,garkiObs101,1,10,tr,5)
  hh2<-plotMod(44,res,garkiObs104,1,11,tr,6)
  hh3<-plotMod(44,res,garkiObs219,2,12,tr,7)
  hh4<-plotMod(44,res,garkiObs220,2,13,tr,8)
  
  png(filename=paste0(fileName,"\\fits.png"),width = w, height = h)
  plot(grid.arrange(hh1,hh2,hh3,hh4))
=======
  png(filename=paste0(fileName,"\\fig6.png"),width = w, height = h,pointsize=30)
 plot(res[,c(16:18)])
  dev.off()
  png(filename=paste0(fileName,"\\fig7.png"),width = w, height = h,pointsize=30)
  plot(res[,c(19:21)])
>>>>>>> Stashed changes
  dev.off()
  png(filename=paste0(fileName,"\\fig8.png"),width = w, height = h,pointsize=30)
  plot(res[,c(22)])
  dev.off()

  
  
  png(filename=paste0(fileName,"\\crossCor.png"),width = w, height = h)
  crosscorr.plot(res)
  dev.off()
  
<<<<<<< Updated upstream
  png(filename=paste0(fileName,"\\density.png"),width = w, height = h)
  plotDens(as.data.frame(res))
  dev.off()
  dev.off()
  
  
=======
  #hh1<-plotMod(44,res,garkiObs101,1,10,tr,5)
  #hh2<-plotMod(44,res,garkiObs104,1,11,tr,6)
  #hh3<-plotMod(44,res,garkiObs219,2,12,tr,7)
  #hh4<-plotMod(44,res,garkiObs220,2,13,tr,8)
  
  #png(filename=paste0(fileName,"\\fits.png"),width = w, height = h)
  #plot(grid.arrange(hh1,hh2,hh3,hh4))

>>>>>>> Stashed changes
  
  print(colMedians(res))
  
  
}

<<<<<<< Updated upstream
cTest<-read.csv("C:\\Imperial\\particleTestLinearN\\35\\results.txt", sep = " ",head =F)
colnames(cTest)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                   ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
colMedians(cTest)
mcmcRes<-cTest


par(mfrow=c(2,2))
p=0.0201241

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki408.txt",sep=",")
test <- test[seq(1, length(test), 4)]
plot(test,col="white",xlab="time",ylab="M")
lines(test)
points(garkiObs408$time,garkiObs408$`408`/p,col="red")

test<-scan("C:\\Imperial\\particleTest\\results2\\4\\garki202.txt",sep=",")
test <- test[seq(1, length(test), 4)]
#plot(test,col="white",xlab="time",ylab="M")
#lines(test)
#points(garkiObs202$time,garkiObs202$`202`/0.01,col="red")

test<-scan("C:\\Imperial\\particleTest\\results2\\4\\garki218.txt",sep=",")
test <- test[seq(1, length(test), 4)]
#plot(test,col="white",xlab="time",ylab="M")
#lines(test)
#points(garkiObs218$time,garkiObs218$`218`/0.001,col="red")

test<-scan("C:\\Imperial\\particleTest\\results2\\4\\garki801.txt",sep=",")
test <- test[seq(1, length(test), 4)]
plot(test,col="white",xlab="time",ylab="M")
lines(test)
points(garkiObs801$time,garkiObs801$`801`/p,col="red")

test<-scan("C:\\Imperial\\particleTest\\results2\\4\\garki553.txt",sep=",")
test <- test[seq(1, length(test), 4)]
plot(test,col="white",xlab="time",ylab="M",ylim=c(0,70000))
lines(test)
points(garkiObs553$time,garkiObs553$`553`/p,col="red")

test<-scan("C:\\Imperial\\particleTest\\results2\\4\\garki802.txt",sep=",")
test <- test[seq(1, length(test), 4)]
plot(test,col="white",xlab="time",ylab="M")
lines(test)
points(garkiObs802$time,garkiObs802$`802`/p,col="red")



par(mfrow=c(2,2))
p=0.0121101

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki408Reff.txt",sep=",")
#test <- test[seq(1, length(test), 4)]
max(test)
plot(test,col="white",xlab="time",ylab="M")
lines(test)
#points(garkiObs408$time,garkiObs408$`408`/p,col="red")

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki202Reff.txt",sep=",")
test <- test[seq(1, length(test), 4)]
#plot(test,col="white",xlab="time",ylab="M")
#lines(test)
#points(garkiObs202$time,garkiObs202$`202`/0.01,col="red")

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki218Reff.txt",sep=",")
test <- test[seq(1, length(test), 4)]
#plot(test,col="white",xlab="time",ylab="M")
#lines(test)
#points(garkiObs218$time,garkiObs218$`218`/0.001,col="red")

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki801Reff.txt",sep=",")
#test <- test[seq(1, length(test), 4)]
max(test)
plot(test,col="white",xlab="time",ylab="M")
lines(test)
#points(garkiObs801$time,garkiObs801$`801`/p,col="red")

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki553Reff.txt",sep=",")
#test <- test[seq(1, length(test), 4)]
max(test)
plot(test,col="white",xlab="time",ylab="M")
lines(test)
#points(garkiObs553$time,garkiObs553$`553`/p,col="red")

test<-scan("C:\\Imperial\\particleTestLinearN\\35\\garki802Reff.txt",sep=",")
#test <- test[seq(1, length(test), 4)]
max(test)
plot(test,col="white",xlab="time",ylab="M")
lines(test)
#points(garkiObs802$time,garkiObs802$`802`/p,col="red")


=======
#test<-read.csv("Q:\\Imperial\\fitPlots\\sTest3.txt",head=F,sep=",")
 #plot(test$V5[26:nrow(test)],col="white")
 #lines(test$V5[26:nrow(test)])

#plotRec(as.mcmc(cTest),as.mcmc(cTest2),"betaBinom 1.5mil power NEW",7)

for (i in c("power")){
fN=i

cTest<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\results.txt"), sep = " ",head =T)
colnames(cTest)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                   ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
ff<-colMedians(cTest)

cTest2<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\results.txt"), sep = " ",head =T)
colnames(cTest2)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                   ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
ss<-colMedians(cTest2)


#png(filename=paste0("Q:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",fN,"\\densityLL.png"),width = 1800, height = 1200)
#plotDensLL(as.data.frame(cTest),as.data.frame(cTest2))
#dev.off()
#png(filename=paste0("Q:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",fN,"Clumped","\\densityLL.png"),width = 1800, height = 1200)
#plotDensLL(as.data.frame(cTest),as.data.frame(cTest2))
#dev.off()

#plotRec(as.mcmc(cTest),as.mcmc(cTest2),fN,7)
#plotRec(as.mcmc(cTest2),as.mcmc(cTest1),paste0(fN,"Clumped"),7)



par(mfrow=c(2,2))
dtt=4
w=    ff[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki408.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]#current discrete time period is 0.25
res<-0
plot(test*w,col="white",xlab="time (days)",ylab="M",bty="n",xlim=c(0,200),main="Village 1",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test*w, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
points(garkiObs408$time,garkiObs408$`408`,type = 'p', pch=16, col = "black")
if(i=="exp"){
legend(130,3900, legend=c("Non-clumped", "Clumped","Observed"),
       col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200),"black"), lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,19),cex=1,bty="n")
}

 if(i=="linear"){
  legend(130,2400, legend=c("Non-clumped", "Clumped","Observed"),
         col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200),"black"), lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,19),cex=1,bty="n")
}

 if(i=="power"){
  legend(130,2700, legend=c("Non-clumped", "Clumped","Observed"),
         col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200),"black"), lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,19),cex=1,bty="n")
}

w=    ss[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki408.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]#current discrete time period is 0.25
res<-0
lines(test*w,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)


w=    ff[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki801.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]
plot(test*w,col="white",xlab="time (days)",ylab="M",bty="n",main="Village 2",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test*w, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
points(garkiObs801$time,garkiObs801$`801`,col="black",type = 'p', pch=16)
#legend(300,365, legend=c("Clumped egg laying", "Non-clumped egg laying","Observed data"),
 #      col=c("red", "blue","black"), lty=c(1,1,NA),pch=c(NA,NA,19), cex=0.9,bty="n")


w=    ss[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki801.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]
lines(test*w,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)


w=    ff[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki553.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]
plot(test*w,col="white",xlab="time (days)",ylab="M",ylim=c(0,2100),bty="n",main="Village 3",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test*w, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
points(garkiObs553$time,garkiObs553$`553`,col="black",type = 'p', pch=16)
#legend(125,2200, legend=c("Clumped egg laying", "Non-clumped egg laying","Observed data"),
 #      col=c("red", "blue","black"), lty=c(1,1,NA),pch=c(NA,NA,19), cex=0.9,bty="n")


w=    ss[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki553.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]
lines(test*w,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)


w=    ff[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki802.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]
plot(test*w,col="white",xlab="time (days)",ylab="M",ylim=c(0,200),bty="n",main="Village 4",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test*w, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
points(garkiObs802$time,garkiObs802$`802`,col="black",type = 'p', pch=16)
#legend(300,200, legend=c("Clumped egg laying", "Non-clumped egg laying","Observed data"),
 #      col=c("red", "blue","black"), lty=c(1,1,NA),pch=c(NA,NA,19), cex=0.9,bty="n")


w=    ss[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki802.txt"),sep=",")
test <- test[seq(1, length(test), dtt)]
lines(test*w,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)


}

#####Rm test calcs #################################


xt<-(colMedians(cTest))
UoE = xt[1]
UoL =xt[2]
dE=xt[15]
dL=xt[16]
n=xt[6]
uP =xt[3]
dP = xt[17]
Y=xt[4]
uM = xt[20]

E=n*(exp(S*uM)-1)/uM

0.5*(n/ (exp(uM*S) - 1))*(1 / (1 + UoE / dE))*(1 / (1 + UoL / dL))*(1 / (1 + (uP) / dP))

##white et al
UoE = 0.034 
UoL =0.035 
dE=6.64 
dL=3.72 
n=21.19
uP =0.25 
dP = 0.64 
Y=13.25 
uM = 0.096 
0.5*((E) / (exp(uM*S) - 1))*(1 / (1 + UoE * dE)) * (1 / (1 + UoL * dL)) * (1 / (1 + uP * dP))


E=n*(exp(S*uM)-1)/uM

#W = -0.5*(y*(UoL / UoE) - (dE / dL) + (y - 1)*UoL*dE) + sqrt(0.25*((y*(UoL / UoE) - (dE / dL) + (y - 1)*UoL*dE))^2 + y*((n*UoL*dE) / (2 * UoE*uM*dL*(1 + dP*uP))))


#uE = UoE*(1+(((E + L) / (K))))
#uL = UoL*(1+(Y*(((E + L) / (K)))))


r<-NULL
for (i in seq(1,20,by=1)){
 g<- n/(exp(uM*i) - 1)
  r<-rbind(r,g)
}
plot(r)



####################plot Rm from c++#################################

for(i in c( "1.5milResults\\power")){
  fN=i
  
par(mfrow=c(2,2))
test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki408Reff.txt"),head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
plot(test,col="white",xlab="time (days)",ylab="Rm",bty="n",main="Village 1",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
if(i=="exp"){
  legend(400,130, legend=c("Non-clumped", "Clumped","Observed"),
         col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200),"black"), lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,19),cex=1,bty="n")
}

if(i=="linear"){
  legend("topright",legend=c("Non-clumped", "Clumped","Observed"),
         col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200),"black"), lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,19),cex=1,bty="n")
}

if(i=="power"){
  legend("topright", legend=c("Non-clumped", "Clumped","Observed"),
         col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200),"black"), lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,19),cex=1,bty="n")
}

test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki408Reff.txt"),head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
res<-0
max(test)
lines(test,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)
abline(h=c(1),lty="dashed")


test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki801Reff.txt"),sep=",",head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
plot(test,col="white",xlab="time (days)",ylab="Rm",bty="n",main="Village 2",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki801Reff.txt"),sep=",",head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
lines(test,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)
abline(h=c(1),lty="dashed")


test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki553Reff.txt"),sep=",",head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
plot(test,col="white",xlab="time (days)",ylab="Rm",bty="n",main="Village 3",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)

test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki553Reff.txt"),sep=",",head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
lines(test,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)
abline(h=c(1),lty="dashed")


test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"\\garki802Reff.txt"),sep=",",head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
plot(test,col="white",xlab="time (days)",ylab="Rm",bty="n",main="Village 3",cex.lab=1.6,cex.main=1.5,cex.axis=1.5)
lines(test, col=rgb(123,50,148,maxColorValue=255, alpha=200),lwd=2)
test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki802Reff.txt"),sep=",",head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
max(test)
res<-0
lines(test,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)
abline(h=c(1),lty="dashed",lwd=1)
}
>>>>>>> Stashed changes




<<<<<<< Updated upstream






res<-NULL
for (i in 1:100){
  print(i)
  o=30
  E=i
  L=10
  K=1000
  uoE=0.035
  uE = uoE*(1+(E + L) / (K))^o;
  
  res<-rbind(res,uE)
}
=======
params <-
  mosParamsP(
    E0 = 24361,
    L0 = 244,
    P0 = 47,
    M0 = 252,
    uoE = 0.03496090 ,
    uoL = 0.03529590 ,
    uP = 0.25343600,
    Y = 13.39410000,
    sf = 10^4.52356000,
    n = 23,
    dE=0.14850950,dL=0.24706200,dP=0.99235600,tr=7.61548000,uM= 0.08081660,Mg=2.65532000,o=0.45936300,
    rF = rFx3,
    time1 = 5/0.25,
    time2 = 186 / 0.25
  )

modR <- larvalR(user = params)#run model
simDat <-
  as.data.frame(modR$run(seq(
    5/0.25, 186 / 0.25
  )))
simDat <- simDat$M[seq(1, length(simDat$M), 4)]
plot(simDat*ff[22])
lines(simDat*ff[22])

test<-scan("Q:\\Imperial\\fitPlots\\sTest3.txt",sep=",")
test <- test[seq(1, length(test), 4)]#current discrete time period is 0.25
plot(test,col="white",xlab="time",ylab="M")
lines(test,col="green")


simDat <- as.data.frame(larvalWhile(params,10,1,0.1))#run model
#simDatH <- simDat$H[seq(1, length(simDat$H), 4)]
#simDatM <- simDat$M[seq(1, length(simDat$M), 4)]
plot(simDat$M,col="white")
lines(simDat$M,col="red")
lines(simDat$H,col="blue")



ggplot(data = cTest2,aes(ll))+
  geom_density(kernel = "gaussian", adjust = 7, aes(colour="Clumped Egg Laying"),fill="red", size=1,alpha=0.1)+
  geom_density(kernel = "gaussian", adjust = 3, aes(cTest$ll,colour="Non-clumped"), size=1,fill="blue",alpha=0.1)+
  scale_colour_manual(values=c("Clumped Egg Laying"="red", "Non-clumped"="blue"), name="Densities")+
  xlab("Log Likelihood")+
  scale_x_continuous(limits=c(-320,-265))+
  theme_bw()


>>>>>>> Stashed changes

