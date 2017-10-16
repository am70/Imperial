

#########################################################################################################################################
#                                                                                                                                       #
#                                                     visualise results and extract C.I.                                                #
#                                                                                                                                       #
#########################################################################################################################################
library(miscTools)

cIfunc<-function(x){
  quantRes<-as.data.frame(c(1:2))
  for(i in 1:(ncol(x)-1)){
    resD<-density(x[,i],adjust=1)
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

plotMod<-function(n,mcmcRes,obsDat,rf,sf,fxd,z){
  parms<-colMedians(mcmcRes)
  cI<-cIfunc(mcmcRes)
  iqr<-iqrFunc(mcmcRes)
  startTime<-min(obsDat$time)/delta
  endTime<-max(obsDat$time)/delta
  
  
#  crosscorr.plot(g$results[,-13])

####95% CI####
  par(mfrow=c(1,1))
 mid<- pFilt(n,iState,modStep3,dataLikFunc,obsDat,pr=c(parms[1],parms[2],parms[3],
                                                                  parms[4],parms[z],parms[9],parms[sf],parms[14]),rFclust=rf,resM=T,fxdParams=fxd)
 colnames(mid)<-c("E","L","P","M")
 #plot(mid$M~timeX$time,col="white")
 #lines(mid$M~timeX$time,col="blue")
#points(obsDat$`101`/0.01~obsDat$time,col="red")


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

hh1<-plotMod(44,pp$results,garkiObs101,1,10,7,5)
hh2<-plotMod(44,pp$results,garkiObs104,1,11,7,6)
hh3<-plotMod(44,pp$results,garkiObs219,2,12,7,7)
hh4<-plotMod(44,pp$results,garkiObs220,2,13,7,8)

grid.arrange(hh1,hh2,hh3,hh4)

resSimMean<-rowMeans(resSimx[,-1])
plot(resSimMean,col="white",ylim=c(0,150))
points(garkiObs$time*10,garkiObs$M,col="red")
lines(resSimMean,col="blue",ylim=c(0,150))

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
#########################################################################################################################################
#                                                                                                                                       #
#                                                     density plots                                                                     #
#                                                                                                                                       #
#########################################################################################################################################

library("gridExtra")

plotDens<-function(mcmc){
  uoEprior<-rnorm( nrow(mcmc),mean=0.035,sd=0.009)
  uoLprior<-rnorm( nrow(mcmc),mean=0.035,sd=0.009)
  uPprior<-rnorm( nrow(mcmc),mean=0.25,sd=0.11)
  Yprior<-rnorm( nrow(mcmc),mean=13.06,sd=4.53)
  #p0prior<-rnorm(94999,mean=0.5,sd=0.015)
  
  
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
  #p5<-ggplot(data = mcmc,aes(o))+
   # geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
  #  geom_density(kernel = "gaussian", adjust = 1, aes(p0prior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
   # scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    #theme_bw()

  grid.arrange(p1, p2, p3, p4)
}





#########################################################################################################################################
#                                                                                                                                       #
#                                                     likelihood test                                                                   #
#                                                                                                                                       #
#########################################################################################################################################

parms<-c(0.03929046 , 0.04366398  ,0.33992035 ,18.40634427 ,10.00000000 , 0.51194118 , 6.67020976 , 0.67765564 , 1.17470130 , 5.36936347  ,2.80765721,  4.05669931 ,0)



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
    
  }
  
  ll
  
}

ps<-parms
parms[10]<-parms[12]
p1<-pFilt(particles,iState,modStep3,dataLikFunc,garkiObs101,pr=parms, resM=T)

dat<-cbind(p1[-1],garkiObsX[-18,])
plot(dat$`p1[-1]`~dat$time,col="white",ylim=c(0,40))
lines(dat$`p1[-1]`~dat$time)
points(dat$`111`/parms[8]~dat$time,col="red")





##########plot and record

plotRec<-function(clustRes,name,tr){
  res<-clustRes$results
  w=1800
  h=1200
  
  fileName<-paste0("Q:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",name)
  dir.create(fileName)
  
  write.table(res,paste0(fileName,"\\mcmcRes.csv"))
  
  
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
#  png(filename=paste0(fileName,"\\fig5.png"),width = w, height = h,pointsize=30)
# plot(res[,c(13)])
 # dev.off()
  
  hh1<-plotMod(44,res,garkiObs101,1,10,tr,5)
  hh2<-plotMod(44,res,garkiObs104,1,11,tr,6)
  hh3<-plotMod(44,res,garkiObs219,2,12,tr,7)
  hh4<-plotMod(44,res,garkiObs220,2,13,tr,8)
  
  png(filename=paste0(fileName,"\\fits.png"),width = w, height = h)
  plot(grid.arrange(hh1,hh2,hh3,hh4))
  dev.off()
  
  png(filename=paste0(fileName,"\\crossCor.png"),width = w, height = h)
  crosscorr.plot(res)
  dev.off()
 
  png(filename=paste0(fileName,"\\density.png"),width = w, height = h)
  plotDens(as.data.frame(res))
  dev.off()
  dev.off()


  
  print(colMedians(res))
  
  
}

f<-obj$task_get("c5ccf1e406561fe5309e0190b5ef0ecd")
plotRec(f$result(),"Z_Inf uoE0.015 50k tr7 short",7)

