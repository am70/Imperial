library(lubridate)
library(odin)
library(parallel)
library(snow)
library(mnormt)

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
rainfall<-subset(rainfall, date >= as.Date("1971-04-27") & date <= as.Date("1973-01-01"))
rainfall$time<-(1:nrow(rainfall))
rFx<-rainfall[rep(seq_len(nrow(rainfall)), each=1/delta),]
rFx<-rFx[,1]

rainfall2<-read.csv("Q:\\Imperial\\larvalModel\\Data\\meteoFUP2.csv",head=F)
colnames(rainfall2)<-c("rainfall","date")
rainfall2$date<-dmy(rainfall2$date)
#limit data to one year
rainfall2<-subset(rainfall2, date >= as.Date("1971-04-27") & date <= as.Date("1973-01-01"))
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

