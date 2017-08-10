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
library(coda)


#########################################################################################################################################
#                                                                                                                                       #
#                                                     load in data                                                                      #
#                                                                                                                                       #
#########################################################################################################################################
delta<-0.25#discrete time period

#simDat<-read.csv("C:\\simDat.csv",sep=" ")
##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("Q:\\Imperial\\larvalModel\\Data\\meteoFUP1.csv",head=F)
colnames(rainfall)<-c("rainfall","date")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1971-01-01") & date <= as.Date("1973-04-27"))
rainfall$time<-(1:nrow(rainfall))
rFx<-rainfall[rep(seq_len(nrow(rainfall)), each=1/delta),]
rFx<-rFx[,1]

rainfall2<-read.csv("Q:\\Imperial\\larvalModel\\Data\\meteoFUP2.csv",head=F)
colnames(rainfall2)<-c("rainfall","date")
rainfall2$date<-dmy(rainfall2$date)
#limit data to one year
rainfall2<-subset(rainfall2, date >= as.Date("1971-01-01") & date <= as.Date("1973-04-27"))
rainfall2$time<-(1:nrow(rainfall2))
rFx2<-rainfall[rep(seq_len(nrow(rainfall)), each=1/delta),]
rFx2<-rFx2[,1]

#load in mosquito data
garki<-read.table("Q:\\Imperial\\larvalModel\\Data\\spraycollect.csv",sep=",",head=T)
garki$ag_sum<-rowSums(garki[,c(9:16)])#create sum of A.gambiae spray samples
garki$Date<-as.Date(garki$Date)
#garki$Date<-dmy(as.character(garki$Date))
garki72<-subset(garki,Date >= as.Date("1972-04-27") & Date <= as.Date("1973-01-01"))
#101,104,113,117,201,203
for (i in c( 101  ,   104   ,  219,  220 )){
  print(i)
  if(i<200) rainTemp<-rainfall else rainTemp<-rainfall2
garki72_101<-subset(garki72,E_Station%in%i) ##subset by village
garki72_101<-merge(garki72_101,rainfall,by.x="Date",by.y="date",all=T) #merge rainfall and mosquito data

garkiObs<-garki72_101[,c(39,37)]
colnames(garkiObs)<-c("time","M")
garkiObs<-subset(garkiObs,M>=0)
garkiObs<-round(aggregate(M~time, data=garkiObs, FUN=mean),0)
#garkiObs<-rbind(data.frame(time = 0, M = 0), garkiObs)
plot(garkiObs$time,garkiObs$M,main=i)
colnames(garkiObs)<-c("time",i)
assign(paste0("garkiObs", i),garkiObs)
}
garkiObsX<-Reduce(function(...) merge(..., all=TRUE), list(garkiObs101, garkiObs104,garkiObs219,garkiObs220))



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
clusterExport(cl, c("rFx","rFx2","delta","garkiObs","nBgP"), envir=environment())


##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                            run pMCM local                                                                      #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

#set.seed(44)
system.time(runX200z5 <- mcmcSampler(initParams = c(uoE=0.03526661,uoL=0.03776644,uP=0.25516500,Y=16.84519894,n=20,p0=0.51322692,o=0.01,Fp=0.97126420,
                                                    o2=1.29593424,sf1=7.31392687,sf2=2.35775867,sf3=10.01558351,sf4=11.56955667)
                                     ,nburn=1000
                                     ,monitoring=2
                                     , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0.1,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1))
                                     , randInit = F
                                     , particles=96
                                     , niter = 10000))


covar=matrix(c(.001,0,0,0,0,0,0,0,0,0,0,0,0,
               0,.001,0,0,0,0,0,0,0,0,0,0,0,
               0,0,.001,0,0,0,0,0,0,0,0,0,0,
               0,0,0,.001,0,0,0,0,0,0,0,0,0,
               0,0,0,0,.001,0,0,0,0,0,0,0,0,
               0,0,0,0,0,.001,0,0,0,0,0,0,0,
               0,0,0,0,0,0,.001,0,0,0,0,0,0,
               0,0,0,0,0,0,0,.001,0,0,0,0,0,
               0,0,0,0,0,0,0,0,.001,0,0,0,0,
               0,0,0,0,0,0,0,0,0,.001,0,0,0,
               0,0,0,0,0,0,0,0,0,0,.001,0,0,
               0,0,0,0,0,0,0,0,0,0,0,.001,0,
               0,0,0,0,0,0,0,0,0,0,0,0,.001),13,13)




#set.seed(44)
system.time(runX200z5 <- mcmcSampler(initParams = c(uoE=0.02783780,uoL=0.04538763,uP=0.25001641,Y=16.86354430,n=60,p0=0.51010071,o=0.09217603,Fp=0.97413706,
                                                    o2=0.91439265,sf1=0.76201783,sf2=0.15444917,sf3=0.85442365,sf4=1.10734896)
                                     ,nburn=1000
                                     ,monitoring=2
                                     , proposer = multiv.proposer(covar)
                                     ,sdProps=c(0.1,0.1,0.1,1,0,0.1,0.1,0.1,0.1,1,1,1,1)
                                     , randInit = F
                                     ,adaptiveMCMC = T
                                     , startAdapt = 150
                                     , particles=126
                                     ,acceptanceRate = 0.3
                                     , niter = 10000))




write.table(runX200z4$results,"Q:\\Imperial\\res2.csv")

meds<-function(x){cbind(median(x[,1]),median(x[,2]),median(x[,3]),median(x[,4]),median(x[,5]),median(x[,6]),median(x[,7]),median(x[,8]),median(x[,9]))}

######################################################################################################################


#test just particle filter
testParams = c(uoE=0.031,uoL=0.035,uP=0.26,Y=17.13,n=25,p0=0.58,Fp=0.3,o=0.03,sf=4.05)
res4<-NULL
system.time(for (i in 1:50){
  ss<-pFilt(35,iState,modStep3,dataLikFunc,garkiObs101,pr=testParams) + lprior(testParams)
  res4<-rbind(ss,res4)
  print(ss)
})

#  testParams = c(uoE=0.031,uoL=0.035,uP=0.26,Y=17.13,n=25,sf=4.05,p0=0.58)

sense<-function(x){
  testParams = c(uoE=0.031,uoL=0.035,uP=x,Y=17.13,n=25,sf=4.05,p0=0.58)
  print(x)
  pFilt(120,simx0,0,modStep3,dataLikFunc,garkiObs,pr=testParams)#+ lprior(testParams)
}

resuP<-lapply(seq(0.001,0.6,by=0.001),sense)

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
                             packages = c("lhs","deSolve","lubridate","odinPackage","dde","buildr","coda","parallel","snow"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",home="//fi--san02/homes/alm210",cores=16,parallel = T))


#########################################################################################################################################
#                                                                                                                                       #
#                                                     run pMCMC on cluster                                                              #
#                                                                                                                                       #
#########################################################################################################################################

runX1_1x1L <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03386983,uoL=0.03652054,uP=0.25801683,Y=14.61421801,n=10,p0=0.51241423,o=0.45443296,Fp=0.99098142,
                                                  o2=1.36693761,sf1=6.09204989,sf2=2.49848543,sf3=9.00862475,sf4=11.17703509)
                                   ,nburn=5000
                                   ,monitoring=2
                                   , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0.1,0.1,0.1,0.01,0.01,0.5,0.5,0.5,0.5))
                                   ,particles = 112
                                 , randInit = F
                                 , niter = 1000000),name="pMMH nx1M")

runX1_1xx1 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03386983,uoL=0.03652054,uP=0.25801683,Y=14.61421801,n=10,p0=0.51241423,o=0.45443296,Fp=0.99098142,
                                                     o2=1.36693761,sf1=6.09204989,sf2=2.49848543,sf3=9.00862475,sf4=11.17703509)
                                    ,nburn=1000
                                    ,monitoring=2
                                    , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,0.1,0.01,0.01,0.5,0.5,0.5,0.5))
                                    ,particles = 96
                                    , randInit = F
                                    , niter = 1000000),name="pMMH n10_1M")



runX1_1y1 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03311636,uoL=0.035,uP=0.25,Y=15.30829479,n=10,p0=0.49791206,o=5,Fp=0.5,
                                                    o2=1,sf1=10,sf2=10,sf3=10,sf4=10)
                                     ,nburn=5000
                                     ,monitoring=2
                                     , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,1,0.01,0.01,0.5,0.5,0.5,0.5))
                                     ,particles = 112
                                     , randInit = F
                                     , niter = 200000),name="pMMH nx")

runX1_2 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03311636,uoL=0.035,uP=0.25,Y=15.30829479,n=10,p0=0.49791206,o=5,Fp=0.5,
                                                  o2=1,sf1=10,sf2=10,sf3=10,sf4=10)
                                   ,nburn=5000
                                   ,monitoring=2
                                   , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,1,0.01,0.01,0.5,0.5,0.5,0.5))
                                   ,particles = 112
                                   , randInit = F
                                   , niter = 100000),name="pMMH n10")


runX1_3 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03311636,uoL=0.035,uP=0.25,Y=15.30829479,n=20,p0=0.49791206,o=5,Fp=0.5,
                                                  o2=1,sf1=10,sf2=10,sf3=10,sf4=10)
                                   ,nburn=5000
                                   ,monitoring=2
                                   , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,1,0.01,0.01,0.5,0.5,0.5,0.5))
                                   ,particles = 112
                                   , randInit = F
                                   , niter = 100000),name="pMMH n20")


runX1_4 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03311636,uoL=0.035,uP=0.25,Y=15.30829479,n=40,p0=0.49791206,o=5,Fp=0.5,
                                                  o2=1,sf1=10,sf2=10,sf3=10,sf4=10)
                                   ,nburn=5000
                                   ,monitoring=2
                                   , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,1,0.01,0.01,0.5,0.5,0.5,0.5))
                                   ,particles = 112
                                   , randInit = F
                                   , niter = 100000),name="pMMH n40")



runX1_5 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03311636,uoL=0.035,uP=0.25,Y=15.30829479,n=60,p0=0.49791206,o=5,Fp=0.5,
                                                  o2=1,sf1=10,sf2=10,sf3=10,sf4=10)
                                   ,nburn=5000
                                   ,monitoring=2
                                   , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,1,0.01,0.01,0.5,0.5,0.5,0.5))
                                   ,particles = 112
                                   , randInit = F
                                   , niter = 100000),name="pMMH n60")


runX1_6 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.03311636,uoL=0.035,uP=0.25,Y=15.30829479,n=80,p0=0.49791206,o=5,Fp=0.5,
                                                  o2=1,sf1=10,sf2=10,sf3=10,sf4=10)
                                   ,nburn=5000
                                   ,monitoring=2
                                   , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0,0.1,1,0.01,0.01,0.5,0.5,0.5,0.5))
                                   ,particles = 112
                                   , randInit = F
                                   , niter = 100000),name="pMMH n80")
