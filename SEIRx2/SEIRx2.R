library(deSolve)
SEIR<-function(t, state, parameters) {
   with(as.list(c(state, parameters)),{
   # rate of change
   N = S+E+I+R
   dS <- -beta*I*S/N
   dE <- beta*I*S/N - E/d1
   dI <- E/d1 - I/d2
   dR <- I/d2
   # return the rate of change
   list(c(dS, dE, dI, dR))
   }) # end with(as.list ...
}  
if (F){
params = c(beta=1.2/5,d1=5, d2=5)
state  = c(S=17000,E=6,I=6,R=0)
times <- seq(0, 730, by = 1)
out <- ode(y = state, times = times, func = SEIR, parms = params, method="adams")
plot(out[,1],out[,2],type="l",col="blue",ylim=c(0,1000),xlab="days",ylab="I")
lines(out[,1],out[,4],type="l",col="red")
lines(out[,1],out[,5],type="l",col="darkgreen")
View(out)
}

SEIR2<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    if (t>tinc){
      b1 = b1inc
      b2 = b2inc
      w  = winc
    }
    # rate of change
    N1 = S1+E1+I1+R1
    N2 = S2+E2+I2+R2
    dS1 = -w*b1*S1*I1/N1 - (1-w)*b1*S1*I2/N2
    dE1 =  w*b1*S1*I1/N1 + (1-w)*b1*S1*I2/N2 - E1/d
    dI1 = E1/d - I1/d
    dC1 = fcrit1*I1/d-C1/dcrit
    dR1 = (1-fcrit1)*I1/d+C1/dcrit
    dS2 = -(b2-b1*(1-w)*N1/N2)*S2*I2/N2 - b1*(1-w)*N1/N2*S2*I1/N1
    dE2 =  (b2-b1*(1-w)*N1/N2)*S2*I2/N2 + b1*(1-w)*N1/N2*S2*I1/N1 - E2/d
    dI2 = E2/d - I2/d
    dC2 = fcrit2*I2/d-C2/dcrit
    dR2 = (1-fcrit2)*I2/d+C2/dcrit
    # return the rate of change
    list(c(dS1, dE1, dI1, dC1, dR1, dS2, dE2, dI2, dC2, dR2))
  }) # end with(as.list ...
}  


# 60% * 340K  = 204K 
# 204 K * 0.5 % = 1020

# 180 K * 0.5 % =  900
# 66  K * 4   % = 2640 

# huidige maatregelen: flatten the peak: OK, maar
# na 1 jaar R0 weer op oude niveau: alsnog epidemie

dur=5
R0 = 1.5
# b1 = 1.50 => 1.455 (assortative) + 0.045 (disass)
# b2 = 0.60 => 0.555 + 0.045
R0inc = 2
# b1 = 2.00 => 1.8 (assortative) + 0.2 (disass)
# b2 = 1.20 => 1.0 + 0.2


N1=10000
N2=7000

params2 = c(d=dur,b1=R0/dur,b2=0.4*R0/dur, w=0.97, tinc=300, b1inc=R0inc/dur, b2inc=0.6*R0inc/dur ,winc=0.90, fcrit1=0.12/100,fcrit2=6/100,dcrit=7)
state2  = c(S1=10000,E1=3.5,I1=3.5,C1=0,R1=0,S2=7000,E2=2.5,I2=2.5,C2=0,R2=0)
times <- seq(0, 700, by = 1)
out <- ode(y = state2, times = times, func = SEIR2, parms = params2, method="adams")
png("fig1.png")
plot(out[,1],out[,4],type="l",col="red",ylim=c(0,300),xlab="dagen",ylab="aantal infecties (duizenden)",lwd=2)
lines(out[,1],out[,9],type="l",lty=2,col="magenta",lwd=2)
legend(450, 300, legend=c("leeftijd < 50", "leeftijd >= 50"),
       col=c("red", "magenta"), lty=1:2, cex=0.8)
dev.off()


png("bedden.png")
plot(out[,1],out[,5]*1000,type="l",col="red",ylim=c(0,2500),xlab="dagen",ylab="benodigd aantal bedden op de IC")
lines(out[,1],out[,10]*1000,type="l",col="magenta",lty=2)
lines(out[,1],(out[,5]+out[,10])*1000,type="l",col="blue",lwd=2)
legend(450, 2400, legend=c("totaal","leeftijd < 50", "leeftijd >= 50"),
       col=c("blue","red", "magenta"), lty=c(1,1,2), cex=0.8)
dev.off()

png("immune.png")
plot(out[,1],out[,6],type="l",col="red",ylim=c(0,6000),xlab="dagen",ylab="hersteld en immuun (duizenden)")
lines(out[,1],out[,11],type="l",col="magenta",lty=2)
legend(450, 6000, legend=c("leeftijd < 50", "leeftijd >= 50"),
       col=c("red", "magenta"), lty=c(1,1,2), cex=0.8)
dev.off()

plot(out[,1],out[,6],type="l",col="red",ylim=c(0,6000),xlab="days",ylab="R")
lines(out[,1],out[,11],type="l",col="magenta")
View(out)


dur=5
params2 = c(d=dur,b1=1.5/dur,b2=0.6/dur,bx=0.15/dur, tinc=300, b1inc=2.0/dur, b2inc=1.2/dur ,bxinc=0.5/dur)
state2  = c(S1=10000,E1=3.5,I1=3.5,R1=0,S2=7000,E2=2.5,I2=2.5,R2=0)
times <- seq(0, 800, by = 1)
out <- ode(y = state2, times = times, func = SEIR2, parms = params2, method="adams")

plot(out[,1],out[,4],type="l",col="red",ylim=c(0,1000),xlab="days",ylab="I1 vs I2")
lines(out[,1],out[,8],type="l",col="magenta")

plot(out[,1],out[,5],type="l",col="red",ylim=c(0,10000),xlab="days",ylab="I1 vs I2")
lines(out[,1],out[,9],type="l",col="magenta")

plot(out[,1],out[,2],type="l",col="blue",ylim=c(0,10000),xlab="days",ylab="I")
lines(out[,1],out[,4],type="l",col="red")
lines(out[,1],out[,5],type="l",col="darkgreen")
lines(out[,1],out[,8],type="l",col="orange")
View(out)

