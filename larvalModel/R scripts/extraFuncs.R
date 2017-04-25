
runFilt<-function(runs){
  
  x<-particleFilter(larvalModP,
                    theta(),
                    init.state,
                    data = garkiObs,
                    n.particles = 1)
  return(x)
}



parRunS<-function(x){parallel::parLapply(NULL,x,runFilt)}