

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



plotMod<-function(n,mcmcRes,obsDat,rf,sf,rFxs,z){
  parms<-colMedians(mcmcRes)
  #cI<-cIfunc(mcmcRes)
  #iqr<-iqrFunc(mcmcRes)
  startTime<-min(obsDat$time)/delta
  endTime<-max(obsDat$time)/delta

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

tauprior<-rnorm( nrow(mcmc),mean=7,sd=1.5)
uMprior<-rnorm( nrow(mcmc),mean=0.091,sd=0.005)
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
  p5<-ggplot(data = mcmc,aes(dE))+
   geom_density(kernel = "gaussian", adjust = 7, aes(colour="Posterior"),fill="red", size=1,alpha=0.1)+
    geom_density(kernel = "gaussian", adjust = 1, aes(dEprior,colour="Prior"), size=1,fill="blue",alpha=0.1)+
   scale_colour_manual(values=c("Posterior"="red", "Prior"="blue"), name="Densities")+
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

plotRec<-function(clustRes,clustRes2,name,tr){
  res<-clustRes
  w=1800
  h=1200
  
  fileName<-paste0("Q:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",name)
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
  png(filename=paste0(fileName,"\\fig6.png"),width = w, height = h,pointsize=30)
 plot(res[,c(16:18)])
  dev.off()
  png(filename=paste0(fileName,"\\fig7.png"),width = w, height = h,pointsize=30)
  plot(res[,c(19:21)])
  dev.off()
  png(filename=paste0(fileName,"\\fig8.png"),width = w, height = h,pointsize=30)
  plot(res[,c(22)])
  dev.off()

  
  
  png(filename=paste0(fileName,"\\crossCor.png"),width = w, height = h)
  crosscorr.plot(res)
  dev.off()
  
  #hh1<-plotMod(44,res,garkiObs101,1,10,tr,5)
  #hh2<-plotMod(44,res,garkiObs104,1,11,tr,6)
  #hh3<-plotMod(44,res,garkiObs219,2,12,tr,7)
  #hh4<-plotMod(44,res,garkiObs220,2,13,tr,8)
  
  #png(filename=paste0(fileName,"\\fits.png"),width = w, height = h)
  #plot(grid.arrange(hh1,hh2,hh3,hh4))

  
  print(colMedians(res))
  
  
}

#test<-read.csv("Q:\\Imperial\\fitPlots\\sTest3.txt",head=F,sep=",")
 #plot(test$V5[26:nrow(test)],col="white")
 #lines(test$V5[26:nrow(test)])

#plotRec(as.mcmc(cTest),as.mcmc(cTest2),"betaBinom 1.5mil power NEW",7)

for (i in c("2.5MilResults\\linear","2.5MilResults\\power","2.5MilResults\\exp")){
fN=i

cTest<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\results.txt"), sep = " ",head =T)
colnames(cTest)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                   ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
ff<-colMedians(cTest)

cTest2<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped\\results.txt"), sep = " ",head =T)
colnames(cTest2)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                   ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
ss<-colMedians(cTest2)


#png(filename=paste0("Q:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",fN,"\\densityLL.png"),width = 1800, height = 1200)
#plotDensLL(as.data.frame(cTest),as.data.frame(cTest2))
#dev.off()
#png(filename=paste0("Q:\\Imperial\\larvalModel\\Figures\\pMMH figures\\",fN,"Clumped","\\densityLL.png"),width = 1800, height = 1200)
#plotDensLL(as.data.frame(cTest),as.data.frame(cTest2))
#dev.off()

plotRec(as.mcmc(cTest),as.mcmc(cTest2),paste0(fN,"NoClumped"),7)
plotRec(as.mcmc(cTest2),as.mcmc(cTest1),paste0(fN,"Clumped"),7)



par(mfrow=c(2,2))
dtt=4
w=    ff[22]
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki408.txt"),sep=",")
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
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki801.txt"),sep=",")
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
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki553.txt"),sep=",")
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
test<-scan(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki802.txt"),sep=",")
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


xt<-(colMedians(cTest2))
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
E=n*(exp(S*uM)-1)/uM
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

for(i in c( "linear","power","exp")){
  fN=i
  
par(mfrow=c(2,2))
test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki408Reff.txt"),head=F)
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
  legend(400,30,  legend=c("Non-clumped", "Clumped"),
         col=c(rgb(123,50,148,maxColorValue=255, alpha=200), rgb(0,136,55,maxColorValue=255, alpha=200)), lty=c(1,1,NA),lwd=c(2,2,NA),cex=1,bty="n")
}

test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"Clumped","\\garki408Reff.txt"),head=F)
test[test == "-nan(ind)"] <- 0
test<-as.vector(as.factor(test$V1))
test<-as.numeric(test)
res<-0
max(test)
lines(test,col=rgb(0,136,55,maxColorValue=255, alpha=200),lwd=2)
abline(h=c(1),lty="dashed")


test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki801Reff.txt"),sep=",",head=F)
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


test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki553Reff.txt"),sep=",",head=F)
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


test<-read.csv(paste0("Q:\\Imperial\\particleTestLinearN\\",fN,"NoClumped","\\garki802Reff.txt"),sep=",",head=F)
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




cTest<-read.csv(paste0("C:\\Users\\ALM210\\Documents\\OverflowResults\\1milResultsFittedTau\\",fN,"NoClumped","\\results.txt"), sep = " ",head =F)
colnames(cTest)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                   ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
ff<-colMedians(cTest)

cTest2<-read.csv(paste0("C:\\Users\\ALM210\\Documents\\OverflowResults\\1milResultsFittedTau\\",fN,"Clumped\\results.txt"), sep = " ",head =F)
colnames(cTest2)<-c("uoE","uoL","uP","Y","w","n","z1","z4","z5","z6"
                    ,"sf1","sf4","sf5","sf6","dE","dL","dP","o","tau","uM","Mg","p","ll")
ss<-colMedians(cTest2)


initialP<-function(ff){
z <- ff[7]#fitted (E-L)
uoE = ff[1]
uoL = ff[2]
uP = ff[3]
y = ff[4]
sf = 10^ff[11]
n = ff[6]
dE=ff[15]
dL=ff[16]
dP=ff[17]
tr=ff[19]
uM= ff[20]
Mg=ff[21]
o=ff[18]
B=21
rF = rFx3
t = 5/0.25

if(tr>t){
  K <- (1 + (sf * ((1 / tr) * (sum(
    rF[0:(t - 1)]
  )))))}else{
  K <- (1 + (sf * ((1 / tr) * (sum(
    rF[(t - tr):(t - 1)]
  )))))
}
#power carrying cap
a = ((n / S) * dP * dL) / ((2 * uM) * (uP * dP))
b = (UoE / (y * UoL)) * (dL + UoL) - dE - UoE
c = -(UoE * dE) / (UoL * y)
x = (-b + sqrt((b^2) * -4 * a * c)) / (2 * a)

L = z
E =round(L / x)
P =round((dL * L) / (uP + dP))
M =round((dP * P) / (2 * uM))

if (M < 1)
  M = 1
return(round(c(E,L,P,M)))
}

startingClumped<-initialP(ss)
startingNoClumped<-initialP(ff)


paramsNoClump <-
  mosParamsP(
    E0 = startingNoClumped[1],
    L0 = startingNoClumped[2],
    P0 = startingNoClumped[3],
    M0 = startingNoClumped[4],
    uoE = ff[1] ,
    uoL = ff[2] ,
    uP = ff[3],
    Y = ff[4],
    sf = 10^ff[11],
    n = ff[6],
    dE=ff[15],dL=ff[16],dP=ff[17],tr=ff[19],uM= ff[20],Mg=ff[21],o=ff[18],
    rF = rFx3,
    time1 = 5/0.25,
    time2 = 250 / 0.25
  )



paramsClump <-
  mosParamsP(
    E0 = startingClumped[1],
    L0 = startingClumped[2],
    P0 = startingClumped[3],
    M0 = startingClumped[4],
    uoE = ss[1] ,
    uoL = ss[2] ,
    uP = ss[3],
    Y = ss[4],
    sf = 10^ss[11],
    n = ss[6],
    dE=ss[15],dL=ss[16],dP=ss[17],tr=ss[19],uM= ss[20],Mg=ss[21],o=ss[18],
    rF = rFx3,
    time1 = 5/0.25,
    time2 = 250 / 0.25
  )



inter<-250
interSizeClump = round(10/8.38098e-03)
interSizeNoClump = round(10/1.01411e-02)

par(mfrow=c(3,3))
for (i in seq(0.1,0.9,by=0.1)){
  
  yDrvPrms<-as.list(c(inter=inter,interSizeClump=0,efficacy=i,clumped=1))
  simDatClump<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(ss[1:6],10^ss[7],10^ss[11],ss[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F,yDrvPrms)
  colnames(simDatClump)<-c("E","L","P","M","G")
  
  yDrvPrms<-as.list(c(inter=inter,interSizeNoClump=0,efficacy=i,clumped=0))
  simDatNoClump<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(ff[1:6],10^ff[7],10^ff[11],ff[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F,yDrvPrms)
  colnames(simDatNoClump)<-c("E","L","P","M","G")
  
  simDatNoClumpM <- simDatNoClump$M[seq(1, length(simDatNoClump$M), 4)]
  simDatClumpM <- simDatClump$M[seq(1, length(simDatClump$M), 4)]
  plot(simDatClumpM,col="white",main=paste0("intervention efficacy =",i),ylab="M",xlab="time (days)")
  lines(simDatClumpM,col="cadetblue1")
  lines(simDatNoClumpM,col="salmon")
  
  yDrvPrms<-as.list(c(inter=inter,interSizeClump=interSizeClump,efficacy=i,clumped=1))
  simDatClump<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(ss[1:6],10^ss[7],10^ss[11],ss[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F,yDrvPrms)
  colnames(simDatClump)<-c("E","L","P","M","G")
  
  yDrvPrms<-as.list(c(inter=inter,interSizeNoClump=interSizeNoClump,efficacy=i,clumped=0))
  simDatNoClump<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(ff[1:6],10^ff[7],10^ff[11],ff[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F,yDrvPrms)
  colnames(simDatNoClump)<-c("E","L","P","M","G")
  
  simDatNoClumpM <- simDatNoClump$M[seq(1, length(simDatNoClump$M), 4)]
  simDatClumpM <- simDatClump$M[seq(1, length(simDatClump$M), 4)]
  
  lines(simDatClumpM,col="blue")
  lines(simDatNoClumpM,col="red")
  abline(v=c(inter*0.25-((5/0.25)/4)),lty="dashed")
  
  #simDatNoClumpM <- simDatNoClump$G[seq(1, length(simDatNoClump$G), 4)]
  #simDatClumpM <- simDatClump$G[seq(1, length(simDatClump$G), 4)]
  #lines(simDatClumpM,col="purple")
  #lines(simDatNoClumpM,col="orange")
  
}


clumpRes<-NULL
noClumpRes<-NULL
for (i in seq(0.01,0.9,by=0.001)){
clumpMG<-ss[21]/0.25
noClumpMg<-ff[21]/0.25
yDrvPrms<-as.list(c(inter=inter,interSizeClump=interSizeClump,efficacy=i,clumped=1))
simDatClump<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(ss[1:6],10^ss[7],10^ss[11],ss[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F,yDrvPrms)
colnames(simDatClump)<-c("E","L","P","M","G")

yDrvPrms<-as.list(c(inter=inter,interSizeNoClump=interSizeNoClump,efficacy=i,clumped=0))
simDatNoClump<- pFilt(25,iState,modStep3,dataLikFunc,obsDat,pr=c(ff[1:6],10^ff[7],10^ff[11],ff[15:22]),rFclust=rf,resM=T,rFx=rFxs,cluster=F,yDrvPrms)
colnames(simDatNoClump)<-c("E","L","P","M","G")

simDatNoClumpM <- simDatNoClump$M[seq(1, length(simDatNoClump$M), 4)]
simDatClumpM <- simDatClump$M[seq(1, length(simDatClump$M), 4)]

eradClumped<-min(which(simDatClumpM[round(inter*0.25):length(simDatClumpM)] < clumpMG))
eradNoClumped<-min(which(simDatNoClumpM[round(inter*0.25):length(simDatNoClumpM)] < noClumpMg))
if(eradClumped==Inf) eradClumped = 0
if(eradNoClumped==Inf) eradNoClumped = 0

clumpRes<-rbind(clumpRes,eradClumped)
print(eradClumped)
noClumpRes<-rbind(noClumpRes,eradNoClumped)
}
par(mfrow=c(1,1))
plot(noClumpRes,col="white")
lines(noClumpRes,col="red")
lines(clumpRes,col="blue")


simDatClump <- as.data.frame(larvalWhileYdrive(paramsClump,10/8.38098e-03,inter,0.1,1,0,TRUE))#run model
simDatNoClump <- as.data.frame(larvalWhileYdrive(paramsNoClump,10/1.01411e-02,inter,0.1,1,0,FALSE))#run model
par(mfrow=c(1,1))

lines(simDatNoClump$rEff,col="red")
abline(v=c(inter-5/0.25),lty="dashed")
plot(simDatClump$M,col="white")
lines(simDatClump$M,col="blue")
lines(simDatNoClump$M,col="red")
abline(v=c(inter-5/0.25),lty="dashed")



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







ggplot(data = cTest2,aes(ll))+
  geom_density(kernel = "gaussian", adjust = 7, aes(colour="Clumped Egg Laying"),fill="red", size=1,alpha=0.1)+
  geom_density(kernel = "gaussian", adjust = 3, aes(cTest$ll,colour="Non-clumped"), size=1,fill="blue",alpha=0.1)+
  scale_colour_manual(values=c("Clumped Egg Laying"="red", "Non-clumped"="blue"), name="Densities")+
  xlab("Log Likelihood")+
  scale_x_continuous(limits=c(-320,-265))+
  theme_bw()



