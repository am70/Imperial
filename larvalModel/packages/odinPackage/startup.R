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

delta <- 0.25#discrete time period
rainfall<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rfall.csv")
rainfall$Date<-as.Date(rainfall$Date)
rainfall <-subset(rainfall,Date >= "1970-01-01" & Date <= "1975-04-27")
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf01.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx1 <- rFx1[, 1]

rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf02.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx2 <- rFx1[, 1]
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf03.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx3 <- rFx1[, 1]
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf04.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx4 <- rFx1[, 1]
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf05.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx5 <- rFx1[, 1]
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf06.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx6 <- rFx1[, 1]
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf07.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx7 <- rFx1[, 1]
rainfall1<-read.csv("Q:\\Imperial\\lModCpp\\Data\\rf08.txt",head=F)
rainfall1$time <- (1:nrow(rainfall1))
rFx1 <- rainfall1[rep(seq_len(nrow(rainfall1)), each = 1 / delta), ]
rFx1$V1<-rFx1$V1*delta
rFx8 <- rFx1[, 1]

#load in mosquito data
garki <-read.table("Q:\\Imperial\\lModCpp\\Data\\anophTable.csv",sep = ",",head = T)
villageEstation <-read.table("Q:\\Imperial\\lModCpp\\Data\\villageEstation.csv",sep = ",",head = T)
colnames(villageEstation)<-c("E_station","village")
villClust <-read.table("Q:\\Imperial\\lModCpp\\Data\\villClust.csv",sep = ",",head = T)
garki<-merge(garki,villageEstation,by.x="E_Station",by.y="E_station")
garki$Date<-as.Date(garki$Date)
garki<-merge(garki,villClust,by.x="village",by.y="Village")
garki<-subset(garki,Date >= "1970-05-28" & Date <= "1972-01-01")
for (i in c(154,   202,  218,  304,553,802)) {
  garkiTemp<-subset(garki,village==i)
  rainfallTemp<-subset(rainfall,fup==head(garkiTemp$Fup,1))
  rainfallTemp$time<-c(1:nrow(rainfallTemp))
  garkiTemp<-aggregate(sumGamFem ~ Date, garkiTemp, sum)
  garkiTemp<-merge(garkiTemp,rainfallTemp,by.x="Date",by.y="Date")
  print(garkiTemp)
  garkiObs<-cbind(garkiTemp$time,garkiTemp$sumGamFem)
  colnames(garkiObs) <- c("time", i)
  plot(garkiObs)
  assign(paste0("garkiObs", i), garkiObs)
}

garkiObsX <-
  Reduce(
    function(...)
      merge(..., all = TRUE),
    list(garkiObs154, garkiObs202, garkiObs218, garkiObs304,garkiObs553,garkiObs802)
  )

#########################################################################################################################################
#                                                                                                                                       #
#                                                     start local cluster                                                               #
#                                                                                                                                       #
#########################################################################################################################################


odin_package("Q:\\Imperial\\larvalModel\\packages\\odinPackage")
load_all()
document()

devtools::install(args = "--no-multiarch")

cl <- makeSOCKcluster(c(7), outfile = "")
clusterEvalQ(cl, library("odinPackage"))
clusterEvalQ(cl, library("lhs"))
clusterEvalQ(cl, library("VGAM"))
clusterEvalQ(cl, library("deSolve"))
clusterExport(cl, c("rFx", "rFx2", "rFx3","delta", "garkiObs", "nBgP"), envir =
                environment())


##################################################################################################################################################
#                                                                                                                                                #
#                                                                                                                                                #
#                                                            run pMMH local                                                                      #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################




## Log-Prior function - specific to this model
# @parms current parameters
# @return sum of log priors
lprior <- function(parms) {
  
  priorSum<-(dnorm(parms[1], mean = 0.035, sd = 0.015, log = T)
             +dunif(parms[1], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[2], mean = 0.035, sd = 0.015, log = T)
             +dunif(parms[2], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[3], mean = 0.25, sd = 0.11, log = T)
             +dunif(parms[3], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[4], mean = 13.06, sd = 4.53, log = T)
            # +dnorm(parms[5], mean = 75, sd = 10, log = T)
             +dunif(parms[6], min = 1e-06, max = 1, log = T)
             +dunif(parms[5], min = 0.01, max = 1e+10, log = T)
             +dunif(parms[7], min = 1, max = 1e+20, log = T)
             +dunif(parms[8],min=1,max=93.6,log=T)
  )
  return(priorSum)
}



#set.seed(44)
system.time(runX200z6x3uoE0.1 <-  mcmcSampler(initParams = c(uoE=0.03,uoL=0.03156214,uP=0.2499843,Y=11.5012,z1=75,z2=75,z3=75,z4=75,w=0.01
                                                             ,sf1=57110,sf2=12660,sf3=8097,sf4=64342,n=50)
                                              ,nburn=5000
                                              ,monitoring=2
                                              , proposer = sequential.proposer
                                              ,sdProps=c(0.01,0.01,0.1,1,1,1,1,1,0.1,3,3,3,3,3)
                                              ,maxSddProps=c(0.1,0.1,0.1,Inf,Inf,Inf,Inf,Inf,0.1,1500000,1500000,1500000,1500000,Inf)
                                              , randInit = F
                                              ,fixedParam=4
                                              ,adaptiveMCMC = T
                                              ,proposerType = 'seq'
                                              , startAdapt = 150
                                              ,tell = 1
                                              , particles=2
                                              ,acceptanceRate =c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.6,0.45,0.6,0.6,0.4)
                                              , niter = 20
                                              ,cluster = F
                                              ,oDat=garkiObsX
                                              ,priorFunc=lprior))

write.table(runX200z4$results,"Q:\\Imperial\\res2.csv")


######################################################################################################################


#test just particle filter
testParams  =  c(uoE=0.03,uoL=0.03156214,uP=0.2499843,Y=11.5012,z1=5000,w=0.01,
                 sf=64342,n=50)
res4<-NULL
system.time(for (i in 1:100){
  ss<-pFilt(50,iState,modStep3,dataLikFunc,garkiObs101,pr = testParams,rFclust = 1,fxdParams = 7,resM = F, cluster = F) #+ lprior(testParams)
  print(ss)
})

#  testParams = c(uoE=0.031,uoL=0.035,uP=0.26,Y=17.13,n=25,sf=4.05,p0=0.58)

sense<-function(x){
  testParams  = c(uoE=0.035,uoL=0.035,uP=0.24359410,Y=17.29277120,Lx=20,w=0.01,
                  o=1.04163691,sf1=100,sf2=100,sf3=100,sf4=100)
  print(x)
  pFilt(120,simx0,0,modStep3,dataLikFunc,garkiObs,pr=testParams)#+ lprior(testParams)
}

resuP<-lapply(seq(0.001,0.6,by=0.001),sense)


ss2<-subset(res4,res4>-100)



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
                             packages = c("lhs","deSolve","lubridate","odinPackage","dde","buildr","coda","parallel","snow","mnormt","rmutil","Rcpp"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="mrc",home="//fi--san02/homes/alm210",cores=64,parallel = T))


#########################################################################################################################################
#                                                                                                                                       #
#                                                     run pMMH on cluster                                                               #
#                                                                                                                                       #
#########################################################################################################################################


## Log-Prior function - specific to this model
# @parms current parameters
# @return sum of log priors
lprior <- function(parms) {
  
  priorSum<-(dnorm(parms[1], mean = 0.035, sd = 0.015, log = T)
             +dunif(parms[1], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[2], mean = 0.035, sd = 0.015, log = T)
             +dunif(parms[2], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[3], mean = 0.25, sd = 0.11, log = T)
             +dunif(parms[3], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[4], mean = 13.06, sd = 10, log = T)
          #   +dnorm(parms[5], mean = 75, sd = 10, log = T)
             +dunif(parms[6], min = 1e-06, max = 1, log = T)
             +dunif(parms[5], min = 0.01, max = 1e+10, log = T)
             +dunif(parms[7], min = 1, max = 1e+20, log = T)
             +dunif(parms[8],min=1,max=93.6,log=T)
  )
  return(priorSum)
}


runXTrn<- obj$enqueue( mcmcSampler(initParams = c(uoE=0.03,uoL=0.03156214,uP=0.2499843,Y=11.5012,z1=75,z2=75,z3=75,z4=75,w=0.01
                                                    ,sf1=57110,sf2=12660,sf3=8097,sf4=64342,n=50)
                                     ,nburn=5000
                                     ,monitoring=2
                                     , proposer = sequential.proposer
                                     ,sdProps=c(0.01,0.01,0.1,1,1,1,1,1,0.1,3,3,3,3,3)
                                     ,maxSddProps=c(0.1,0.1,0.1,Inf,Inf,Inf,Inf,Inf,0.1,1500000,1500000,1500000,1500000,Inf)
                                     , randInit = F
                                     ,fixedParam=4
                                     ,adaptiveMCMC = T
                                     ,proposerType = 'seq'
                                     , startAdapt = 150
                                     ,tell = 20
                                     , particles=75
                                     ,acceptanceRate =c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.6,0.45,0.6,0.6,0.4)
                                     , niter = 250000
                                     ,cluster = T
                                     ,oDat=garkiObsX
                                     ,priorFunc=lprior),name="pMMH Z_Infy9 uoE0.015 250k zMax tr4 64")




runXTrn50<- obj$enqueue( mcmcSampler(initParams = c(uoE=0.03,uoL=0.03156214,uP=0.2499843,Y=11.5012,z1=75,z2=75,z3=75,z4=75,w=0.01
                                                    ,sf1=57110,sf2=12660,sf3=8097,sf4=64342,n=50)
                                     ,nburn=5000
                                     ,monitoring=2
                                     , proposer = sequential.proposer
                                     ,sdProps=c(0.01,0.01,0.1,1,1,1,1,1,0.1,3,3,3,3,3)
                                     ,maxSddProps=c(0.1,0.1,0.1,Inf,Inf,Inf,Inf,Inf,0.1,1500000,1500000,1500000,1500000,Inf)
                                     , randInit = F
                                     ,fixedParam=14
                                     ,adaptiveMCMC = T
                                     ,proposerType = 'seq'
                                     , startAdapt = 150
                                     ,tell = 20
                                     , particles=50
                                     ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.45,0.6,0.6,0.5)
                                     , niter = 50000
                                     ,cluster = T
                                     ,oDat=garkiObsX
                                     ,priorFunc=lprior),name="pMMH Z_Inf uoE0.015 50k tr14")





## Log-Prior function - specific to this model
# @parms current parameters
# @return sum of log priors
lprior <- function(parms) {
  
  priorSum<-(dnorm(parms[1], mean = 0.035, sd = 0.015, log = T)
             +dunif(parms[1], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[2], mean = 0.035, sd = 0.015, log = T)
             +dunif(parms[2], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[3], mean = 0.25, sd = 0.11, log = T)
             +dunif(parms[3], min = 0.001, max = 0.99, log = T)
             +dnorm(parms[4], mean = 13.06, sd = 4.53, log = T)
             +dunif(parms[5], min = 1, max = 1e+5, log = T)
                +dnorm(parms[5], mean = 75, sd = 10, log = T)
             +dunif(parms[6], min = 1e-06, max = 1, log = T)
             +dunif(parms[5], min = 0.01, max = 100000, log = T)
             +dunif(parms[7], min = 1, max = 1e+20, log = T)
             +dunif(parms[8],min=1,max=93.6,log=T)
  )
  return(priorSum)
}


runXTrn50<- obj$enqueue( mcmcSampler(initParams = c(uoE=0.03,uoL=0.03156214,uP=0.2499843,Y=11.5012,z1=75,z2=75,z3=75,z4=75,w=0.01
                                                    ,sf1=57110,sf2=12660,sf3=8097,sf4=64342,n=50)
                                     ,nburn=5000
                                     ,monitoring=2
                                     , proposer = sequential.proposer
                                     ,sdProps=c(0.01,0.01,0.1,1,1,1,1,1,0.1,3,3,3,3,3)
                                     ,maxSddProps=c(0.1,0.1,0.1,Inf,Inf,Inf,Inf,Inf,0.1,1500000,1500000,1500000,1500000,Inf)
                                     , randInit = F
                                     ,fixedParam=7
                                     ,adaptiveMCMC = T
                                     ,proposerType = 'seq'
                                     , startAdapt = 150
                                     ,tell = 20
                                     , particles=50
                                     ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.45,0.6,0.6,0.5)
                                     , niter = 50000
                                     ,cluster = T
                                     ,oDat=garkiObsX,
                                     priorFunc=lprior),name="pMMH Z_75_sd10 uoE0.015 50k tr7")




runXTrn50<- obj$enqueue( mcmcSampler(initParams = c(uoE=0.03,uoL=0.03156214,uP=0.2499843,Y=11.5012,z1=75,z2=75,z3=75,z4=75,w=0.01
                                                    ,sf1=57110,sf2=12660,sf3=8097,sf4=64342,n=50)
                                     ,nburn=5000
                                     ,monitoring=2
                                     , proposer = sequential.proposer
                                     ,sdProps=c(0.01,0.01,0.1,1,1,1,1,1,0.1,3,3,3,3,3)
                                     ,maxSddProps=c(0.1,0.1,0.1,Inf,Inf,Inf,Inf,Inf,0.1,1500000,1500000,1500000,1500000,Inf)
                                     , randInit = F
                                     ,fixedParam=14
                                     ,adaptiveMCMC = T
                                     ,proposerType = 'seq'
                                     , startAdapt = 150
                                     ,tell = 20
                                     , particles=50
                                     ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.45,0.6,0.6,0.5)
                                     , niter = 50000
                                     ,cluster = T
                                     ,oDat=garkiObsX,
                                     priorFunc=lprior),name="pMMH Z_75_sd10 uoE0.015 50k tr14")


