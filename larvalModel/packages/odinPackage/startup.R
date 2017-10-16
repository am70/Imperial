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

#simDat<-read.csv("Q:\\simDat.csv",sep=" ")
##load in Garki rainfall data NOT YET VILLAGE SPECIFIC
rainfall <-
  read.csv("Q:\\Imperial\\larvalModel\\Data\\meteoFUP1.csv", head = F)
colnames(rainfall) <- c("rainfall", "date")
rainfall$date <- dmy(rainfall$date)
#limit data to one year
rainfall <-
  subset(rainfall,
         date >= as.Date("1971-01-01") & date <= as.Date("1973-04-27"))
rainfall$time <- (1:nrow(rainfall))
rFx <- rainfall[rep(seq_len(nrow(rainfall)), each = 1 / delta), ]
rFx <- rFx[, 1]

rainfall2 <-
  read.csv("Q:\\Imperial\\larvalModel\\Data\\meteoFUP2.csv", head = F)
colnames(rainfall2) <- c("rainfall", "date")
rainfall2$date <- dmy(rainfall2$date)
#limit data to one year
rainfall2 <-
  subset(rainfall2,
         date >= as.Date("1971-01-01") & date <= as.Date("1973-04-27"))
rainfall2$time <- (1:nrow(rainfall2))
rFx2 <- rainfall2[rep(seq_len(nrow(rainfall2)), each = 1 / delta), ]
rFx2 <- rFx2[, 1]

#load in mosquito data
garki <-
  read.table(
    "Q:\\Imperial\\larvalModel\\Data\\spraycollect.csv",
    sep = ",",
    head = T
  )
garki$ag_sum <-
  rowSums(garki[, c(9:16)])#create sum of A.gambiae spray samples
garki$Date <- as.Date(garki$Date)
garki72 <-
  subset(garki,#1972-05-27
         Date >= as.Date("1972-05-30") & Date <= as.Date("1973-01-01"))
for (i in c(101  ,   104   ,  219,  220)) {
  print(i)
  if (i < 200)
    rainTemp <- rainfall
  else
    rainTemp <- rainfall2
  garki72_101 <- subset(garki72, E_Station %in% i) ##subset by village
  garki72_101 <-
    merge(
      garki72_101,
      rainfall,
      by.x = "Date",
      by.y = "date",
      all = T
    ) 
  #merge rainfall and mosquito data
  garkiObs <- garki72_101[, c(39, 37)]
  colnames(garkiObs) <- c("time", "M")
  garkiObs <- subset(garkiObs, M >= 0)
  garkiObs <- round(aggregate(M ~ time, data = garkiObs, FUN = mean), 0)
  #garkiObs<-rbind(data.frame(time = 0, M = 0), garkiObs)
  plot(garkiObs$time, garkiObs$M, main = i)
  colnames(garkiObs) <- c("time", i)
  assign(paste0("garkiObs", i), garkiObs)
}
garkiObsX <-
  Reduce(
    function(...)
      merge(..., all = TRUE),
    list(garkiObs101, garkiObs104, garkiObs219, garkiObs220)
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
clusterExport(cl, c("rFx", "rFx2", "delta", "garkiObs", "nBgP"), envir =
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
             +dunif(parms[5], min = 1, max = 1e+5, log = T)
            # +dnorm(parms[5], mean = 75, sd = 10, log = T)
             +dunif(parms[6], min = 1e-06, max = 1, log = T)
             +dunif(parms[5], min = 0.01, max = 100000, log = T)
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
                                              ,fixedParam=7
                                              ,adaptiveMCMC = T
                                              ,proposerType = 'seq'
                                              , startAdapt = 150
                                              ,tell = 1
                                              , particles=50
                                              ,acceptanceRate =c(0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.45,0.6,0.6,0.5)
                                              , niter = 30000
                                              ,cluster = F
                                              ,oDat=garkiObsX
                                              ,priorFunc=lprior))

write.table(runX200z4$results,"Q:\\Imperial\\res2.csv")


######################################################################################################################


#test just particle filter
testParams  = c(3.000000e-02, 3.156214e-02 ,2.499843e-01 ,1.150120e+01 ,3.513033e+01, 1.000000e-02 ,7.827497e+00 ,3.771100e+05, 5.000000e+01)
res4<-NULL
system.time(for (i in 1:100){
  ss<-pFilt(50,iState,modStep3,dataLikFunc,garkiObs101,pr = testParams,rFclust = 1,fxdParams = 7,resM = F, cluster = F)#+ lprior(testParams)
  res4<-rbind(ss,res4)
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
                             packages = c("lhs","deSolve","lubridate","odinPackage","dde","buildr","coda","parallel","snow","mnormt","rmutil"))

obj <- didehpc::queue_didehpc(ctx,didehpc::didehpc_config(cluster="didemrchnb",home="//fi--san02/homes/alm210",cores=32,parallel = T))


#########################################################################################################################################
#                                                                                                                                       #
#                                                     run pMMH on cluster                                                               #
#                                                                                                                                       #
#########################################################################################################################################

minSpd=c(0.001,0.001,0.003,0.001,0.001,0.001,0.001,0.00001,0.03)





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
          #   +dnorm(parms[5], mean = 75, sd = 10, log = T)
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
                                     ,oDat=garkiObsX
                                     ,priorFunc=lprior),name="pMMH Z_Inf uoE0.015 50k tr7")




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


