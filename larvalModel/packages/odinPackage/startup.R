library(lhs)
library(lubridate)
library(ggplot2)
library(VGAM)
library(odin)
library(parallel)
library(devtools)
library(snow)
library(pse)
library(provisionr)


#########################################################################################################################################
#                                                                                                                                       #
#                                                     load in data                                                                      #
#                                                                                                                                       #
#########################################################################################################################################
simDat<-read.csv("Q:\\simDat.csv",sep=" ")
##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("Q:\\Imperial\\larvalModel\\Data\\garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
rainfall$time<-(1:nrow(rainfall))

#load in mosquito data
garki<-read.table("Q:\\Imperial\\larvalModel\\Data\\spraycollect.csv",sep=",",head=T)
garki$ag_sum<-rowSums(garki[,c(9:16)])#create sum of A.gambiae spray samples
garki$Date<-as.Date(garki$Date)
#garki$Date<-dmy(as.character(garki$Date))
garki72<-subset(garki,Date >= as.Date("1973-05-27") & Date <= as.Date("1974-06-03"))

garkiSubs<-c(101)
garki72_101<-subset(garki72,E_Station%in%garkiSubs) ##subset by village - needs checking with michael specific villagse
garki72_101<-merge(garki72_101,rainfall,by.x="Date",by.y="date",all=T) #merge rainfall and mosquito data
#garki72_101$ag_sum[is.na(garki72_101$ag_sum)]<-0


garkiObs<-garki72_101[,c(39,37)]
colnames(garkiObs)<-c("time","M")
garkiObs<-subset(garkiObs,M>=0)
garkiObs<-round(aggregate(M~time, data=garkiObs, FUN=mean),0)
garkiObs<-rbind(data.frame(time = 0, M = 0), garkiObs)

delta<-0.1#discrete time period
rF<-as.data.frame(rainfall$rainfall)
rFx<-rF[rep(seq_len(nrow(rF)), each=1/delta),] #split rainfall data into discreet time periods

#########################################################################################################################################
#                                                                                                                                       #
#                                                     start local cluster                                                               #
#                                                                                                                                       #
#########################################################################################################################################

odin_package("Q:\\Imperial\\larvalModel\\packages\\odinPackage")
load_all()
document()

devtools::install(args = "--no-multiarch")

cl <- makeSOCKcluster(c(7),outfile="")
clusterEvalQ(cl, library("odinPackage"))
clusterEvalQ(cl, library("lhs"))
clusterEvalQ(cl, library("VGAM"))
clusterEvalQ(cl, library("deSolve"))
clusterExport(cl, c("rFx","delta","garkiObs"), envir=environment())


##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                            run pMCM local                                                                      #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

#set.seed(44)
system.time(runX200z3 <- mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.25,Y=13,n=25,sf=4,p0=0.5)
                                     ,nburn=2
                                     ,monitoring=2
                                     , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.0,0.01,0.01))
                                     , randInit = F
                                     ,particles=98
                                     , niter = 100))



write.table(runX200,"C:\\res.csv")



#test just particle filter
testParams = c(uoE=0.031,uoL=0.035,uP=0.26,Y=17.13,n=25,sf=4.05,p0=0.58)
res4<-NULL
for (i in 1:100){
  ss<-pFilt(5,simx0,0,modStep3,dataLikFunc,garkiObs,pr=testParams)+ lprior(testParams)
  res4<-rbind(ss,res4)
  print(ss)
}

#  testParams = c(uoE=0.031,uoL=0.035,uP=0.26,Y=17.13,n=25,sf=4.05,p0=0.58)

sense<-function(x){
  testParams = c(uoE=0.031,uoL=0.035,uP=x,Y=17.13,n=25,sf=4.05,p0=0.58)
  print(x)
  pFilt(120,simx0,0,modStep3,dataLikFunc,garkiObs,pr=testParams)#+ lprior(testParams)
}

resuP<-lapply(seq(0.01,0.6,by=0.001),sense)

#########################################################################################################################################
#                                                                                                                                       #
#                                                     start dide cluster                                                                #
#                                                                                                                                       #
#########################################################################################################################################


context::context_log_start()
root <- "contexts"
obj$cluster_load(TRUE)

sources<-c("Q:\\Imperial\\larvalModel\\packages\\odinPackage\\data_functions.R")

ctx <- context::context_save(root,
                             package_source=provisionr::package_sources(
                             local="Q:\\Imperial\\larvalModel\\packages\\odinPackage"),
                             sources=sources,
                             packages = c("lhs","VGAM","deSolve","lubridate","odinPackage","dde","buildr","coda"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",home="//fi--san02/homes/alm210",cores=16,parallel = T))


#########################################################################################################################################
#                                                                                                                                       #
#                                                     run pMCMC on cluster                                                              #
#                                                                                                                                       #
#########################################################################################################################################

runX1_1 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.25,Y=13,n=100,sf=4,p0=0.5)
                                 ,nburn=1000
                                 ,monitoring=0
                                 ,particles=128
                                 , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.0,0.01,0.01))
                                 , randInit = F
                                 , niter = 1000000),name="pMCMC n100")


#########################################################################################################################################
#                                                                                                                                       #
#                                                     visualise results and extract C.I.                                                #
#                                                                                                                                       #
#########################################################################################################################################
library(miscTools)

medFunc<-function(x){
  res<-x$result()
  colMedians(res$results)
}

cIfunc<-function(x){
  quantRes<-as.data.frame(c(1:2))
  res<-x$result()
  for(i in 1:(ncol(res$results)-1)){
    resD<-density(res$results[,i])
    samp <- sample(resD$x, 1e6, replace = TRUE, prob = resD$y)
    quants<-as.data.frame(quantile(samp, c(0.05, 0.95)))
    colnames(quants)<-c(names(res$initParams[i]))
    quantRes<-cbind(quantRes,quants)
  }
  quantRes[,-1]
}


plotMod<-function(mcmcRes){
parms<-medFunc(mcmcRes)
cI<-cIfunc(mcmcRes)

resSimx<-as.vector(pFilt(200,simx0,0,modStep3,dataLikFunc,garkiObs,
                  pr=c(uoE=parms[1],uoL=parms[2],uP=parms[3],Y=parms[4],n=parms[5],sf=parms[6],p0=0.6),resM=T))

resSimxL<- as.vector(pFilt(200,simx0,0,modStep3,dataLikFunc,garkiObs,
                  pr=c(uoE=cI[1,1],uoL=cI[1,2],uP=cI[1,3],Y=cI[1,4],n=cI[2,5],sf=cI[1,6],p0=0.6),resM=T))

resSimxH<- as.vector(pFilt(200,simx0,0,modStep3,dataLikFunc,garkiObs,
                 pr=c(uoE=cI[2,1],uoL=cI[2,2],uP=cI[2,3],Y=cI[2,4],n=cI[1,5],sf=cI[2,6],p0=0.6),resM=T))

rMean<-as.data.frame(cbind(resSimx,resSimxL,resSimxH,garkiObs$time))
colnames(rMean)<-c("M","uM","lM","time")

ggplot(data = rMean,aes(rMean$M))+
  geom_point(x=garkiObs$time,y=garkiObs$M,col="red")+
  expand_limits(y=c(0,150))+
  geom_ribbon(aes(x=rMean$time, ymax=rMean$lM, ymin=rMean$uM), fill="grey", alpha=.5) +
  geom_line(aes(x=time,y = rMean$lM), colour = 'grey') +
  geom_line(aes(x=time,y = rMean$uM), colour = 'grey')+
  geom_line(aes(x=time,rMean$M))+
  theme_bw()

}


plotMod(runX1_7)


resSimMean<-rowMeans(resSimx[,-1])
plot(resSimMean,col="white",ylim=c(0,150))
points(garkiObs$time*10,garkiObs$M,col="red")
lines(resSimMean,col="blue",ylim=c(0,150))
