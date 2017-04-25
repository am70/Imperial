
likeFunc<-function(pr){ #add parameters
  parameters<-read.table(text = pr, sep = ",", colClasses = "numeric")
  reps<-parameters$V7
  logx<-NULL
  modReff<-NULL
  while(reps>0){
  modS <- modelRun(parms=mosParams(uoE=parameters$V1,uoL=parameters$V2,uP=parameters$V3,
                                    Y=parameters$V4,sf=parameters$V5,n=parameters$V6),obsDat = garkiObs)
  matchedTimes <- modS$step %in% garkiObs$time

  log_like <- sum(-dzipois(garkiObs$M,modS$M[matchedTimes], log = TRUE))
  reps<-reps-1
  logx<-cbind(log_like,logx)
  modReff<-cbind(modS$Reff,modReff)
  }
  
  logx <- logx[!is.na(logx)]
  modReff <- modReff[!is.na(modReff)]

  if (!is.na(log_like)){
    res<-c(parameters$V1,parameters$V2,parameters$V3,parameters$V4,parameters$V5,parameters$V6,sum(logx)/length(logx),sum(modReff)/length(modReff))
  }
  return(res)
}

