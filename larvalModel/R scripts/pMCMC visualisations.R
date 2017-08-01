

#########################################################################################################################################
#                                                                                                                                       #
#                                                     visualise results and extract C.I.                                                #
#                                                                                                                                       #
#########################################################################################################################################
library(miscTools)

n100<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\n100.csv",sep=" ")
n80<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\n80.csv",sep=" ")
n60<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\n60.csv",sep=" ")
n40<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\n40.csv",sep=" ")
n20<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\n20.csv",sep=" ")
n10<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\n10.csv",sep=" ")


nx<-read.csv("Q:\\Imperial\\larvalModel\\Outputs\\pMCMC\\nx.csv",sep=" ")



medFunc<-function(x){
  colMedians(x)
}

cIfunc<-function(x){
  quantRes<-as.data.frame(c(1:2))
  res<-x
  for(i in 1:(ncol(res)-1)){
    resD<-density(res[,i])
    samp <- sample(resD$x, 1e6, replace = TRUE, prob = resD$y)
    quants<-as.data.frame(quantile(samp, c(0.05, 0.95)))
    colnames(quants)<-c(names(res[i]))
    quantRes<-cbind(quantRes,quants)
  }
  quantRes[,-1]
}

plotMod<-function(mcmcRes){
  parms<-medFunc(mcmcRes)
  cI<-cIfunc(mcmcRes)
  resSimx<-c(3110:5519)
  resSimxL<-c(3110:5519)
  resSimxH<-c(3110:5519)
  
  for (i in 1:150){
    mod <- odinPackage::larvalModP(user=mosParamsP(uoE=parms[1],uoL=parms[2],uP=parms[3],
                                                   Y=parms[4],n=parms[5],sf=parms[10],o=parms[9],tr=14,time1=3110,E0=f[1],L0=f[2],P0=f[3],M0=f[4])) 
    sim <- as.data.frame(mod$run(3110:5519))
   # plot(sim$M)
    #plot(rainfall$time*10,rainfall$rainfall,col="white")
   # lines(rainfall$time*10,rainfall$rainfall,col="blue")
    
    #lines(sim$M,ylim=c(0,50),col="white")
   # plot(garkiObs104$time[-1]*10,garkiObs104$`104`[-1],col="red",ylim=c(0,30))
   # lines(sim$M~sim$timeX)
    resSimx<- cbind(resSimx,sim$M)
    
    mod <- odinPackage::larvalModP(user=mosParamsP(uoE=cI[1,1],uoL=cI[1,2],uP=cI[1,3],
                                                   Y=cI[1,4],n=cI[1,5],sf=parms[10],tr=14,o=parms[9],E0=f[1],L0=f[2],P0=f[3],M0=f[4],time1=3110)) 
    sim <- as.data.frame(mod$run(3110:5519))
    resSimxL<- cbind(resSimxL,sim$M)
    
    mod <- odinPackage::larvalModP(user=mosParamsP(uoE=cI[2,1],uoL=cI[2,2],uP=cI[2,3],
                                                   Y=cI[2,4],n=cI[2,5],sf=parms[10],tr=14,o=parms[9],E0=f[1],L0=f[2],P0=f[3],M0=f[4],time1=3110)) 
    sim <- as.data.frame(mod$run(3110:5519))
    resSimxH<- cbind(resSimxH,sim$M)
  }
  resSimMean<-rowMeans(resSimx[,-1])
  resSimLMean<-rowMeans(resSimxL[,-1])
  resSimHMean<-rowMeans(resSimxH[,-1])
  
  time<-c(3110:5519)
  garkiObs101$timeX<-garkiObs101$time*10
  rMean<-as.data.frame(cbind(resSimMean,resSimLMean,resSimHMean,time))
  rMean2<-merge(garkiObs101,rMean,by.x="timeX",by.y="time", all=T)[-1,]
  
  ggplot(data = rMean,aes(rMean$resSimMean))+
   geom_point(x=rMean2$timeX,y=rMean2$`101`,col="red")+
    expand_limits(y=c(0,30))+
    geom_ribbon(aes(x=time, ymax=rMean$resSimLMean, ymin=rMean$resSimHMean), fill="grey", alpha=.5)+
    geom_line(aes(x=time,y = rMean$resSimLMean), colour = 'grey') +
    geom_line(aes(x=time,y = rMean$resSimHMean), colour = 'grey')+
    geom_line(aes(x=time,rMean$resSimMean))+
    xlab("time")+
    ylab("M")+
    theme_bw()
  
}

hh1<-plotMod(res)
grid.arrange(hh1,hh2,hh3,hh4)

resSimMean<-rowMeans(resSimx[,-1])
plot(resSimMean,col="white",ylim=c(0,150))
points(garkiObs$time*10,garkiObs$M,col="red")
lines(resSimMean,col="blue",ylim=c(0,150))

simX<-c(1:2440)
for (i in 1:100){
  
  mod <- odinPackage::larvalModP(user=mosParamsP(uoE=parms[1],uoL=parms[2],uP=parms[3],
                                                 Y=parms[4],n=parms[5],sf=parms[10],o=parms[9],tr=14,time1=3110,E0=f[1],L0=f[2],P0=f[3],M0=f[4]))
  sim <- as.data.frame(mod$run(1:2440))
  simX<-cbind(simX,sim$M)

}
# 101  ,   104   ,  108  ,   113
simMean<-rowMeans(simX[,-1])
plot(simMean~sim$timeX,ylim=c(0,30),col="white")
lines(simMean~sim$timeX)
points(garkiObs101$time*10,(garkiObs101$`101`/parms[8]),col="red")
#########################################################################################################################################
#                                                                                                                                       #
#                                                     density plots                                                                     #
#                                                                                                                                       #
#########################################################################################################################################



library("gridExtra")

plotDens<-function(mcmc){
  uoEprior<-rnorm(998999,mean=0.035,sd=0.00485)
  uoLprior<-rnorm(998999,mean=0.035,sd=0.00485)
  uPprior<-rnorm(998999,mean=0.25,sd=0.0357)
  Yprior<-rnorm(998999,mean=13.06,sd=3.53)
  
  p1<-ggplot(data = mcmc,aes(uoE))+
    geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 5, aes(uoEprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p2<-ggplot(data = mcmc,aes(uoL))+
    geom_density(kernel = "gaussian", adjust = 11, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 5, aes(uoLprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p3<-ggplot(data = mcmc,aes(uP))+
    geom_density(kernel = "gaussian", adjust = 5, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 5, aes(uPprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()
  p4<-ggplot(data = mcmc,aes(Y))+
    geom_density(kernel = "gaussian", adjust = 13, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 5, aes(Yprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
    scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
    theme_bw()

  grid.arrange(p1, p2, p3, p4)
}

plotDens(nx)
plotDens(n100)
plotDens(n80)
plotDens(n60)
plotDens(n40)
plotDens(n20)
plotDens(n10)





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
