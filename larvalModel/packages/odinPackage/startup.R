library(lhs)
library(lubridate)
library(ggplot2)
library(rmutil)
library(odin)
library(parallel)
library(devtools)
library(snow)
library(pse)
library(provisionr)
library(coda)
library(mnormt)
library(rmutil)
#########################################################################################################################################
#                                                                                                                                       #
#                                                     load in data                                                                      #
#                                                                                                                                       #
#########################################################################################################################################
delta<-0.25#discrete time period

#simDat<-read.csv("Q:\\simDat.csv",sep=" ")
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
rFx2<-rainfall2[rep(seq_len(nrow(rainfall2)), each=1/delta),]
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
system.time(runX200z6x3 <- mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                    o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                     ,nburn=1000
                                     ,monitoring=2
                                     ,proposer = sequential.proposer
                                     ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                     ,maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                     ,randInit = F
                                     ,fixedParam=50
                                     ,adaptiveMCMC = T
                                     ,proposerType = 'seq'
                                     ,startAdapt = 10
                                     ,particles=15
                                     ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                     ,niter = 10000
                                     ,cluster =T
                                     ))




write.table(runX200z4$results,"Q:\\Imperial\\res2.csv")


######################################################################################################################


#test just particle filter
testParams  = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                           o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
res4<-NULL
system.time(for (i in 1:100){
  ss<-pFilt(44,iState,modStep3,dataLikFunc,garkiObs101,pr=testParams,fxdParams=50,rFclust=1) + lprior(testParams)
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

odin_package("Q:\\Imperial\\larvalModel\\packages\\odinPackage")
load_all()
document()

devtools::install(args = "--no-multiarch")


context::context_log_start()
root <- "contexts"
obj$cluster_load(TRUE)

sources<-c("Q:\\Imperial\\larvalModel\\packages\\odinPackage\\data_functions.R")

ctx <- context::context_save(root,
                             package_source=provisionr::package_sources(
                             local="Q:\\Imperial\\larvalModel\\packages\\odinPackage"),
                             sources=sources,
                             packages = c("lhs","deSolve","lubridate","odinPackage","dde","buildr","coda","parallel","snow","mnormt","rmutil"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",home="//fi--san02/homes/alm210",cores=16,parallel = T))


#########################################################################################################################################
#                                                                                                                                       #
#                                                     run pMCMC on cluster                                                              #
#                                                                                                                                       #
#########################################################################################################################################

minSpd=c(0.001,0.001,0.003,0.001,0.001,0.001,0.001,0.00001,0.03)

runX60Inf <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                           o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                            ,nburn=10000
                                            ,monitoring=2
                                            , proposer = sequential.proposer
                                            ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                            , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                            , randInit = F
                                            ,fixedParam=60
                                            ,adaptiveMCMC = F
                                            , startAdapt = 150
                                            , particles=25
                                            ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                            , niter = 100000
                                            ,cluster = T),name="pMMH n60 100k Inf")
  

runX60 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                    o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                     ,nburn=10000
                                     ,monitoring=2
                                     , proposer = sequential.proposer
                                     ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                  , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                  , randInit = F
                                     ,fixedParam=60
                                     ,adaptiveMCMC = F
                                     , startAdapt = 150
                                     , particles=25
                                     ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                     , niter = 100000
                                     ,cluster = T),name="pMMH n60 1mil")

runX50 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                 o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                  ,nburn=10000
                                  ,monitoring=2
                                  , proposer = sequential.proposer
                                  ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                  , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                  , randInit = F
                                  ,fixedParam=50
                                  ,adaptiveMCMC = F
                                  , startAdapt = 150
                                  , particles=25
                                  ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                  , niter = 1000000
                                  ,cluster = T),name="pMMH n50 1mil")


runX40 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                 o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                  ,nburn=10000
                                  ,monitoring=2
                                  , proposer = sequential.proposer
                                  ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                  , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                  , randInit = F
                                  ,fixedParam=40
                                  ,adaptiveMCMC = F
                                  , startAdapt = 150
                                  , particles=25
                                  ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                  , niter = 1000000
                                  ,cluster = T),name="pMMH n40 1mil")

runX30 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                 o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                  ,nburn=10000
                                  ,monitoring=2
                                  , proposer = sequential.proposer
                                  ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                  , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                  , randInit = F
                                  ,fixedParam=30
                                  ,adaptiveMCMC = F
                                  , startAdapt = 150
                                  , particles=25
                                  ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                  , niter = 1000000
                                  ,cluster = T),name="pMMH n30 1mil")

runX20 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                 o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                  ,nburn=10000
                                  ,monitoring=2
                                  , proposer = sequential.proposer
                                  ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                  , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                  , randInit = F
                                  ,fixedParam=20
                                  ,adaptiveMCMC = F
                                  , startAdapt = 150
                                  , particles=25
                                  ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                  , niter = 1000000
                                  ,cluster = T),name="pMMH n20 1mil")

runX10 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                                                 o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
                                  ,nburn=10000
                                  ,monitoring=2
                                  , proposer = sequential.proposer
                                  ,sdProps=c(0.1,0.1,0.1,1,1,0.1,0.1,3,3,3,3)
                                  , maxSddProps=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
                                  , randInit = F
                                  ,fixedParam=10
                                  ,adaptiveMCMC = F
                                  , startAdapt = 150
                                  , particles=25
                                  ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.3,0.6,0.45,0.6,0.6)
                                  , niter = 1000000
                                  ,cluster = T),name="pMMH n10 1mil")
