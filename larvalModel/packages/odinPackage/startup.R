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
system.time(runX200z2 <- mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.25,Y=13,n=25,sf=4,p0=0.5)
                                     ,nburn=2
                                     ,monitoring=2
                                     , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.0,0.01,0.01))
                                     , randInit = F
                                     , niter = 100))



write.table(runX200,"C:\\res.csv")



#test just particle filter
testParams = c(uoE=0.03156379,uoL=0.03581366,uP=0.26077825,Y=17.13787936,n=25.0000000020,sf=4.05889653,p0=0.58250050)
set.seed(10)
res3<-NULL
for (i in 1:10){
  ss<-pFilt(200,simx0,0,modStep3,dataLikFunc,garkiObs,pr=testParams)+ lprior(testParams)
  res3<-rbind(ss,res3)
  print(ss)
}


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

runX1 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.25,Y=13,n=10,sf=4,p0=0.5)
                                 ,nburn=2000
                                 ,monitoring=0
                                 ,particles=128
                                 , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.0,0.01,0.01))
                                 , randInit = F
                                 , niter = 100000),name="pMCMC")

