#param
dE <- user()
dL <- user()
dP <- user()
uoE <- user()
uoL <- user()
uP <- user()
uM <- user()
Y <- user()
S <- user()
tr <- user()
sf <- user()
dt<-user()
n<-user()
rF[] <- user()
dim(rF) <- user()
trx<-(tr/dt)
Emax<-user()

#initial values
initial(Be)<-0
initial(Bl)<-0
initial(Bp)<-0

initial(E) <- 177
initial(L) <- 8
initial(P) <- 1
initial(M) <- 7
initial(Bm)<-0
initial(nt)<-0
initial(Reff)<-0

K <-if (step<=trx) (1+(sf*((1/trx)))) else (1+(sf*((1/trx)*(sum(rF[(step-trx):step])))))

uE<-uoE*dt*(1+((E+L)/(K)))
uL<-uoL*dt*(1+(Y*(E+L)/(K)))

update(Be)<-rbinom(E,(dE+uE)*dt)
update(Bl)<-rbinom(L, (dL+uL)*dt)
update(Bp)<-rbinom(P,(dP+uP)*dt)
update(Bm)<-rbinom(M,uM*dt)
update(nt)<-rbinom(M,(dt/S))

update(Reff)<-0.5*(Emax/(exp(uM*S)-1))*(1/(1+uE/dE))*(1/(1+uL/dL))*(1/(1+(uP*dt)/dP))

update(E)<-if(E-Be>0)E-Be+rpois(nt*n) else rpois(nt*n)
update(L)<-if(L-Bl>0)L-Bl+rbinom(Be,(dE/(uE+dE))) else rbinom(Be,(dE/(uE+dE)))
update(P)<-if(P-Bp>0)P-Bp+rbinom(Bl,(dL/(uL+dL))) else rbinom(Bl,(dL/(uL+dL)))
update(M)<-if(M-Bm>0)M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))-Bm else M+(0.5*(rbinom(Bp,(dP/(uP+dP)))))

