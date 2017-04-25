modelRun <-
  function(parms, obsDat) {
    modR <- larvalMod(user=parms)
    simDat <- as.data.frame(modR$run(0:2000))
    simDat<-simDat[seq(1, NROW(simDat), by = 1/delta),]
    simDat$step<-simDat$step*delta
    matchedTimes <- simDat$step %in% garkiObs$time
    simDat$M[simDat$M<=0|is.na(simDat$M)] = 1e-3
    return(simDat)
  }
