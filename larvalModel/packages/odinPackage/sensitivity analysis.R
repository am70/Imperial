

#########################################################################################################################################
#                                                                                                                                       #
#                                                     sensitivity analysis                                                              #
#                                                                                                                                       #
#########################################################################################################################################

library(boot)
library(sensitivity)
library(gridExtra)
#########################################################################################################################################
#                                                                                                                                       #
#                                                              LHS                                                                      #
#                                                                                                                                       #
#########################################################################################################################################

HC<-NULL
##create latin hypercube
HCtemp<- randomLHS(10000, 6)
HC$uoE <- qnorm(HCtemp[,1],mean=0.04193237,sd=0.0485)
HC$uoL <-  qnorm(HCtemp[,2],mean=0.04355339,sd=0.0485)
HC$uP <-  qnorm(HCtemp[,3],mean=0.137071,sd=0.0357)
HC$Ys <-qnorm(HCtemp[,4],mean=13.06,sd=1.997)
HC$sf <-  qunif(HCtemp[,5],min=1,max=100)
HC$n<-30#qnorm(HCtemp[,6],mean=12.75,sd=2)
HC2<-HC
HC<-as.data.frame(HC)
#(24.86-0.64)/3.92

HC<- paste(HC$uoE,HC$uoL,HC$uP, HC$Y,HC$sf, HC$n, 50,sep=",")#concatonate into single vector 
system.time(lhcResults<-data.frame(t(sapply(tempRes<-parLapply(cl,HC,likeFunc), `[`))))#run in parallel and convert results to dataframe

colnames(lhcResults)<-c("uoE","uoL","uP","Y","sf","n","logLike","Reff")

mins<-lhcResults[which(lhcResults$logLike == min(lhcResults$logLike)), ]#find minimum values

modS <- modelRun(parms=mosParams(uoE=mins$uoE,uoL=mins$uoL,uP=mins$uP,
                                  Y=mins$Y,sf=mins$sf,n=mins$n),obsDat=garkiObs)

plot(garkiObs$M~garkiObs$time,col="red")
lines(modS$M)

write.table(lhcResults,"C:\\Users\\Aaron\\Documents\\ClumpyNewResultsConstrainedBandY.csv")

write.table(lhcResults,"C:\\Users\\ALM210\\Dropbox (SPH Imperial College)\\deterministicModResults.csv")



#########################################################################################################################################
#                                                                                                                                       #
#                                                    prcc                                                                               #
#                                                                                                                                       #
#########################################################################################################################################

g<-cbind(as.data.frame(HC2),lhcResults[,c(7:8)])
#colnames(g)<-c("uoE","uoL","uP","Y","sf","n","logLike","Reff")
f<-g#subset(g,Reff>0) B

pccRes <- pcc(f[,c(1:6)],f$logLike, nboot = 50,rank=T)
#
plot(pccRes)
abline(h=0,lty=2)

pccrDat<-pccRes$PRCC
pccrDat$names<-rownames(pccrDat)
pccrDat$rows<-c(1:nrow(pccrDat))

colnames(pccrDat)<-c("prcc","bias","std.error","min","max","names","rows")
#pccrDat$names<-c("uoE","uoL","uP","Y","sf","B")

ggplot(pccrDat,aes(x=names,y=prcc))+
  xlab("Parameters")+
  expand_limits(y=c(-1,1))+
   ggtitle("Partial (rank) correlation coefficient - logLike - non clumpy laying")+
  geom_errorbar(data=pccrDat,aes(ymin=min, ymax=max), width=.1)+
  geom_point(data=pccrDat,aes(x=names, y=prcc), color='blue',alpha=0.5)+
   geom_hline(yintercept=0,lty=2)+
  theme_bw()
  



#########################################################################################################################################
#                                                                                                                                       #
#                                                     scatterplots of LHS                                                               #
#                                                                                                                                       #
#########################################################################################################################################

df<-f#cbind(lhcResults,as.data.frame(HC2))
df<-subset(f,B>=0)
#colnames(df)<-c("uoE","uoL","uP","Y","sf","B","logLike","Reff")
p<-ggplot(data=df)+ 
    ylab("logLike")+
# coord_cartesian(ylim=c(0,90))+
    theme_bw()
  
p1<-p+geom_point(data=df,aes(x=uoE, y=log(logLike)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=uoE, y=log(logLike)),se = FALSE,col="black",lty=2)
p2<-p+geom_point(data=df,aes(x=uoL, y=log(logLike)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=uoL, y=log(logLike)),se = FALSE,col="black",lty=2)
p3<-p+geom_point(data=df,aes(x=uP, y=log(logLike)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=uP, y=log(logLike)),se = FALSE,col="black",lty=2)
p4<-p+geom_point(data=df,aes(x=n, y=log(logLike)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=n, y=log(logLike)),se = FALSE,col="black",lty=2)
p5<-p+geom_point(data=df,aes(x=Ys, y=log(logLike)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=Ys, y=log(logLike)),se = FALSE,col="black",lty=2)
p6<-p+geom_point(data=df,aes(x=sf, y=log(logLike)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=sf, y=log(logLike)),se = FALSE,col="black",lty=2)

grid.arrange(p1,p2,p3,p4,p5,p6)



df<-f#cbind(lhcResults,as.data.frame(HC2))
colnames(df)<-c("uoE","uoL","uP","Y","sf","n","logLike","Reff")
p<-ggplot(data=df)+ 
  ylab("Reff")+
   coord_cartesian(ylim=c(0,120))+
  theme_bw()

p1<-p+geom_point(data=df,aes(x=uoE, y=(Reff)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=uoE, y=(Reff)),se = FALSE,col="black",lty=2)
p2<-p+geom_point(data=df,aes(x=uoL, y=(Reff)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=uoL, y=(Reff)),se = FALSE,col="black",lty=2)
p3<-p+geom_point(data=df,aes(x=uP, y=(Reff)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=uP, y=(Reff)),se = FALSE,col="black",lty=2)
p4<-p+geom_point(data=df,aes(x=n, y=(Reff)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=n, y=(Reff)),se = FALSE,col="black",lty=2)
p5<-p+geom_point(data=df,aes(x=Y, y=(Reff)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=Y, y=(Reff)),se = FALSE,col="black",lty=2)
p6<-p+geom_point(data=df,aes(x=sf, y=(Reff)), color='blue',alpha=0.5)+
  geom_smooth(method = "lm",aes(x=sf, y=(Reff)),se = FALSE,col="black",lty=2)

grid.arrange(p1,p2,p3,p4,p5,p6)

