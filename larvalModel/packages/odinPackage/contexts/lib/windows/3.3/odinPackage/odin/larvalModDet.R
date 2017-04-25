##Parameters
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
B<-user()
rF[] <- user()
dim(rF) <- user()
trx<-(tr/dt)
Emax<-user()

## Derivatives
deriv(E) <- (B*M)-(E*dE)-((uoE*((1+(E+L)/(K))))*E)
deriv(L) <- (E*dE)-(L*dL)-((uoL*(1+(Y*(E+L)/(K))))*L)
deriv(P) <- (L*dL)-(P*dP)-(uP*P)
deriv(M) <- 0.5*(P/dP)-(uM*M)
deriv(day)<-+1

uE<- (uoE*((1+(E+L)/(K))))
uL <- (uoL*(1+(Y*(E+L)/(K))))

K <-if (day<=trx) (1+(sf*((1/trx)))) else (1+(sf*((1/trx)*(sum(rF[(day-trx):day])))))


deriv(Reff)<-0.5*(Emax/(exp(uM*S)-1))*(1/(1+uE*dE))*(1/(1+uL*dL))*(1/(1+uP*dP))-Reff

## Initial conditions
initial(E) <- 177
initial(L) <- 8
initial(P) <- 1
initial(M) <- 7
initial(Reff)<-0
initial(day)<-1