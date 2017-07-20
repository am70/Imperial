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
#simDat<-read.csv("C:\\simDat.csv",sep=" ")
##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("C:\\Imperial\\larvalModel\\Data\\meteoFUP1.csv",head=F)
colnames(rainfall)<-c("rainfall","date")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1971-04-27") & date <= as.Date("1972-01-01"))
rainfall$time<-(1:nrow(rainfall))

#load in mosquito data
garki<-read.table("C:\\Imperial\\larvalModel\\Data\\spraycollect.csv",sep=",",head=T)
garki$ag_sum<-rowSums(garki[,c(9:16)])#create sum of A.gambiae spray samples
garki$Date<-as.Date(garki$Date)
#garki$Date<-dmy(as.character(garki$Date))
garki72<-subset(garki,Date >= as.Date("1971-04-27") & Date <= as.Date("1972-01-01"))
#101,104,113,117,201,203
for (i in c( 101  ,   104   ,  111  ,   113)){
  print(i)
garki72_101<-subset(garki72,E_Station%in%i) ##subset by village
garki72_101<-merge(garki72_101,rainfall,by.x="Date",by.y="date",all=T) #merge rainfall and mosquito data

garkiObs<-garki72_101[,c(39,37)]
colnames(garkiObs)<-c("time","M")
garkiObs<-subset(garkiObs,M>=0)
garkiObs<-round(aggregate(M~time, data=garkiObs, FUN=mean),0)
garkiObs<-rbind(data.frame(time = 0, M = 0), garkiObs)
plot(garkiObs$time,garkiObs$M,main=i)
colnames(garkiObs)<-c("time",i)
assign(paste0("garkiObs", i),garkiObs)
}
garkiObsX<-Reduce(function(...) merge(..., all=TRUE), list(garkiObs101, garkiObs104,garkiObs111,garkiObs113))

delta<-0.1#discrete time period

rFx<-rainfall[rep(seq_len(nrow(rainfall)), each=1/delta),]
rFx<-rFx[,1]

#########################################################################################################################################
#                                                                                                                                       #
#                                                     start local cluster                                                               #
#                                                                                                                                       #
#########################################################################################################################################

odin_package("C:\\Imperial\\larvalModel\\packages\\odinPackage")
load_all()
document()

devtools::install(args = "--no-multiarch")

cl <- makeSOCKcluster(c(7),outfile="")
clusterEvalQ(cl, library("odinPackage"))
clusterEvalQ(cl, library("lhs"))
clusterEvalQ(cl, library("VGAM"))
clusterEvalQ(cl, library("deSolve"))
clusterExport(cl, c("rFx","delta","garkiObs","nBgP"), envir=environment())


##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                            run pMCM local                                                                      #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

#set.seed(44)
system.time(runX200z4 <- mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.25,Y=13,n=10,p0=0.5,o=0.5,sf1=4,sf2=4.1,sf3=4.2,sf4=4.3)
                                     ,nburn=1000
                                     ,monitoring=2
                                     , proposer = sequential.proposer(sdProps=c(0.001,0.001,0.01,0.1,0.1,0.01,0.01,0.1,0.1,0.1,0.1))
                                     , randInit = F
                                     ,particles=35
                                     , niter = 25000))



write.table(runX200z4$results,"C:\\res.csv")

meds<-function(x){cbind(median(x[,1]),median(x[,2]),median(x[,3]),median(x[,4]),median(x[,5]),median(x[,6]),median(x[,7]),median(x[,8]),median(x[,9]))}

######################################################################################################################


#test just particle filter
testParams = c(uoE=0.031,uoL=0.035,uP=0.26,Y=17.13,n=25,p0=0.58,o=0.03,sf=4.05)
res4<-NULL
system.time(for (i in 1:50){
  ss<-pFilt(5,iState,modStep3,dataLikFunc,garkiObs101,pr=testParams)+ lprior(testParams)
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

sources<-c("C:\\Imperial\\larvalModel\\packages\\odinPackage\\data_functions.R")

ctx <- context::context_save(root,
                             package_source=provisionr::package_sources(
                             local="C:\\Imperial\\larvalModel\\packages\\odinPackage"),
                             sources=sources,
                             packages = c("lhs","VGAM","deSolve","lubridate","odinPackage","dde","buildr","coda","parallel","snow"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",home="//fi--san02/homes/alm210",cores=16,parallel = T))


#########################################################################################################################################
#                                                                                                                                       #
#                                                     run pMCMC on cluster                                                              #
#                                                                                                                                       #
#########################################################################################################################################

runX1_1 <- obj$enqueue(mcmcSampler(initParams = c(uoE=0.035,uoL=0.035,uP=0.25,Y=13,n=100,sf=4,p0=0.5)
                                 ,nburn=5
                                 ,monitoring=0
                                 ,particles=64
                                 , proposer = sequential.proposer(sdProps=c(0.01,0.01,0.01,0.1,0.0,0.01,0.01))
                                 , randInit = F
                                 , niter = 50),name="pMCMC n100")
