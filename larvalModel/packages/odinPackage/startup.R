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
rainfall<-read.csv("Q:\\Imperial\\larvalModel\\Data\\garkiRainfall.csv",head=F)
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

odin_package("Q:\\Imperial\\larvalModel\\packages\\odinPackage")
load_all()
document()

devtools::install(args = "--no-multiarch")

cl <- makeSOCKcluster(c(7),outfile="")
clusterEvalQ(cl, library("odinPackage"))
clusterEvalQ(cl, library("lhs"))
clusterEvalQ(cl, library("VGAM"))
clusterEvalQ(cl, library("deSolve"))
clusterExport(cl, c("rFx","delta","garkiObs","data.point"), envir=environment())


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
                             package_source=provisionr::package_sources(github=c("richfitz/odin","richfitz/dde"),
                             local="Q:\\Imperial\\larvalModel\\packages\\odinPackage"),
                             sources=sources,
                             packages = c("lhs","VGAM","deSolve","odin","lubridate","odinPackage","dde"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",cores = 16,home="//fi--san02/homes/alm210"))



#################################particle filter tests###################################

p1<-obj$enqueue(parRunS(1:16,500000,seq(0.002, 0.002, length.out = 500000),"uoE"),name="particle filter 100k1")
p2<-obj$enqueue(parRunS(17:32,500000,seq(0.002, 0.002, length.out = 500000),"uoE"),name="particle filter 100k2")
p3<-obj$enqueue(parRunS(33:48,500000,seq(0.002, 0.002, length.out = 500000),"uoE"),name="particle filter 100k3")
p4<-obj$enqueue(parRunS(49:64,500000,seq(0.002, 0.002, length.out = 500000),"uoE"),name="particle filter 100k4")
p5<-obj$enqueue(parRunS(65:80,500000,seq(0.002, 0.002, length.out = 500000),"uoE"),name="particle filter 100k5")
p6<-obj$enqueue(parRunS(81:96,500000,seq(0.002, 0.002, length.out = 500000),"uoE"),name="particle filter 100k6")
pRes100k<-c(unlist(p1$result()),unlist(p2$result()),unlist(p3$result()),unlist(p4$result()),unlist(p5$result()),unlist(p6$result()))


puoE2<-obj$enqueue(parRunS(1:500,1000,seq(0.002, 0.002, length.out = 100),"uoE"),name="particle filter uoE")
puoL<-obj$enqueue(parRunS(1:16,2000,seq(0.01, 0.01, length.out = 16),"uoL"),name="particle filter uoL")
puP<-obj$enqueue(parRunS(1:1000,5000,seq(0.001, 0.5, length.out = 1000),"uP"),name="particle filter uP")
pn<-obj$enqueue(parRunS(1:500,5000,seq(0, 20, length.out = 500),"n"),name="particle filter n")
psf<-obj$enqueue(parRunS(1:500,5000,seq(0, 20, length.out = 500),"sf"),name="particle filter sf")
pY<-obj$enqueue(parRunS(1:500,5000,seq(0, 20, length.out = 500),"y"),name="particle filter Y")


##################particle filter MCMC tester#######################################

particleFilterMCMC(fitmodel=larvalModP, 
                   theta=fit.params,
                   init.state = init.state,
                   data = garkiObs,
                   n.particles = 200000)


###########################plots######################################

puoEplot<-ggplot()+
  xlab("uoE")+
  ylab("Marginal log Like")+
  geom_point(aes(y=unlist(puoE$result()), x=seq(0.001, 0.1, length.out = 500)), color='red',alpha=0.9)+
  theme_bw()


puoLplot<-ggplot()+
  xlab("uoL")+
  ylab("Marginal log Like")+
  geom_point(aes(y=unlist(puoL$result()), x=seq(0.001, 0.1, length.out = 500)), color='red',alpha=0.9)+
  theme_bw()

puPplot<-ggplot()+
  xlab("uP")+
  ylab("Marginal log Like")+
  geom_point(aes(y=unlist(puP$result()), x=seq(0.001, 0.3, length.out = 1000)), color='red',alpha=0.9)+
  theme_bw()


psfplot<-ggplot()+
  xlab("sf")+
  ylab("Marginal log Like")+
  geom_point(aes(y=unlist(psf$result()), x=seq(0, 20, length.out = 500)), color='red',alpha=0.9)+
  theme_bw()

pnplot<-ggplot()+
  xlab("n")+
  ylab("Marginal log Like")+
  geom_point(aes(y=unlist(pn), x=seq(0, 20, length.out = 500)), color='red',alpha=0.9)+
  theme_bw()

pYplot<-ggplot()+
  xlab("Y")+
  ylab("Marginal log Like")+
  geom_point(aes(y=unlist(pY$result()), x=seq(0, 20, length.out = 500)), color='red',alpha=0.9)+
  theme_bw()


resAllplot<-ggplot()+
  xlab("particles")+
  ylab("S.D.")+
  geom_point(aes(y=resAllMSD$sSD, x=resAllMSD$i), color='red',alpha=0.9)+
  geom_line(aes(y=resAllMSD$sSD, x=resAllMSD$i), color='red',alpha=0.9)+
  theme_bw()

#(sf*(1/(trx*(1-exp(-step/trx)))))*(sum(rF[(step-trx):step]))
#(sf*((1/trx)*(sum(rF[(step-trx):step]))))


resAll<-NULL
resAll<-rbind(p1kR,p2kR,p3kR,p4kR,p5kR,p6kR,p7kR,p8kR,p9kR,
              pRes10k,pRes20k,pRes30k,pRes40k,pRes50k,pRes100k,pRes200k,pRes1Mil)


resAllMSD<-NULL
for (i in unique(resAll$particles)){
  
  s<-subset(resAll,particles==i)
  sMean<-mean(s$ll)
  sSD<-sd(s$ll)
  d<-cbind(sMean,sSD,i)
  resAllMSD<-rbind(resAllMSD,d)
}




