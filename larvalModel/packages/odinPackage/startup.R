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

##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall<-read.csv("Q:\\ImperialMalaria\\larvalModel\\Data\\garkiRainfall.csv",head=F)
colnames(rainfall)<-c("date","rainfall")
rainfall$date<-dmy(rainfall$date)
#limit data to one year
rainfall<-subset(rainfall, date >= as.Date("1973-05-27") & date <= as.Date("1974-06-03"))
rainfall$time<-(1:nrow(rainfall))

#load in mosquito data
garki<-read.table("Q:\\ImperialMalaria\\larvalModel\\Data\\spraycollect.csv",sep=",",head=T)
garki$ag_sum<-rowSums(garki[,c(8:16)])#create sum of A.gambiae spray samples
garki$Date<-as.Date(garki$Date)
#garki$Date<-dmy(as.character(garki$Date))
garki72<-subset(garki,Date >= as.Date("1973-05-27") & Date <= as.Date("1974-06-03"))
garki72_101<-subset(garki72,E_Station==101) ##subset by village - needs checking with michael specific villagse
garki72_101<-merge(garki72_101,rainfall,by.x="Date",by.y="date",all=T) #merge rainfall and mosquito data
#garki72_101$ag_sum[is.na(garki72_101$ag_sum)]<-0


garkiObs<-garki72_101[,c(39,37)]
colnames(garkiObs)<-c("time","M")
garkiObs<-subset(garkiObs,M>=0)



delta<-0.1#discrete time period
rF<-as.data.frame(rainfall$rainfall)
rFx<-rF[rep(seq_len(nrow(rF)), each=1/delta),] #split rainfall data into discreet time periods

#########################################################################################################################################
#                                                                                                                                       #
#                                                     start local cluster                                                               #
#                                                                                                                                       #
#########################################################################################################################################

odin_package("Q:\\ImperialMalaria\\larvalModel\\packages\\odinPackage")
load_all()
document()

devtools::install(args = "--no-multiarch")

cl <- makeSOCKcluster(c(7),outfile="")
clusterEvalQ(cl, library("odinPackage"))
clusterEvalQ(cl, library("lhs"))
clusterEvalQ(cl, library("VGAM"))
clusterEvalQ(cl, library("deSolve"))
clusterExport(cl, c("rFx","delta","garkiObs"), envir=environment())


#########################################################################################################################################
#                                                                                                                                       #
#                                                     start dide cluster                                                                #
#                                                                                                                                       #
#########################################################################################################################################


context::context_log_start()
root <- "contexts"
obj$cluster_load(TRUE)

sources<-c("Q:\\ImperialMalaria\\larvalModel\\packages\\odinPackage\\data_functions.R")

ctx <- context::context_save(root,
                             package_source=provisionr::package_sources(github=c("richfitz/odin","richfitz/dde"),
                             local="Q:\\ImperialMalaria\\larvalModel\\packages\\odinPackage", expire = 0.01),
                             sources=sources,
                             packages = c("lhs","VGAM","deSolve","odin","lubridate","odinPackage","dde"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",cores = 16,home="//fi--san02/homes/alm210"))

vals=seq(0.001,0.04,length.out=10)

p1<-obj$enqueue(parRunS(1:16),name="particle filter 1.1")
p2<-obj$enqueue(parRunS(17:32),name="particle filter 2.1")
p3<-obj$enqueue(parRunS(33:48),name="particle filter 3.1")
p4<-obj$enqueue(parRunS(49:64),name="particle filter 4.1")
p5<-obj$enqueue(parRunS(65:80),name="particle filter 5.1")
p6<-obj$enqueue(parRunS(81:96),name="particle filter 6.1")

pRes2<-c(unlist(p1$result()),unlist(p2$result()),unlist(p3$result()),unlist(p4$result()),unlist(p5$result()),unlist(p6$result()))
