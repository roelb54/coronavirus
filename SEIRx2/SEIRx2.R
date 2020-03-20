# code by Roel Bakker
# roel.bakker@gmail.com


library(deSolve)

SEIR2 = function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    if (t>tinc){
      b1 = b1inc
      b2 = b2inc
      w  = winc
    }
    N1 = S1+E1+I1+R1
    N2 = S2+E2+I2+R2
    labda1 = b1*(w*I1/N1 + (1-w)*I2/N2)
    dS1 = -labda1*S1
    dE1 =  labda1*S1 - E1/dinc
    dI1 =  E1/dinc - I1/dinf
    dH1 =    fh1 *I1/dinf    - fc1 *H1/(0.75*dhosp) - (1-fc1)*H1/dhosp
    dC1 =                      fc1 *H1/(0.75*dhosp)        - C1/dicu
    dR1 = (1-fh1)*I1/dinf + (1-fc1)*H1/dhosp        + (1-mu1)*C1/dicu
    dD1 =  mu1*C1/dicu
    labda2 = (b2-b1*(1-w)*N1/N2)*I2/N2 + (1-w)*b1*N1/N2*I1/N1
    dS2 = -labda2*S2
    dE2 =  labda2*S2 - E2/dinc
    dI2 =  E2/dinc   - I2/dinf
    dH2 =  fh2*I2/dinf - fc2 *H2/(0.75*dhosp) - (1-fc2)*H2/dhosp
    dC2 =                fc2 *H2/(0.75*dhosp) - C2/dicu
    dR2 = (1-fh2)*I2/dinf + (1-fc2)*H2/dhosp        + (1-mu2)*C2/dicu
    dD2 =  mu2*C2/dicu
    # return the rate of change
    list(c(dS1, dE1, dI1, dH1, dC1, dR1, dD1, dS2, dE2, dI2, dH2, dC2, dR2, dD2))
  }) 
}  
dur = 4.6
R0young  = 0.7 * 0.7 * 2.84 # 0.70 * 0.70
R0old    = 0.6 * 0.7 * 1.75 # 0.60 * 0.70 
R0young2 = 1.0 * 0.7 * 2.84 # 0.70 * 0.70
R0old2   = 0.7 * 0.7 * 1.75
# INTERVENTIE
# b1 = 1.39 => 0.97 * 1.39 = 1.35 (assortative) + 0.03 * 1.39 = 0.042 (disass)
# b2 = 0.735 => 0.693 (assor) + 0.042 (disass)
# NA 300 DAGEN
# b1 = 1.99 => 1.72 (assortative) + 0.27 (disass)
# b2 = 0.86 => 0.59 + 0.27
params2 = c(dinc=dur,dinf=dur,dhosp=8,dicu=10,fh1=2.03/100,fh2=16.88/100,fc1=5.674/100,fc2=0.355,mu1=0.5,mu2=0.5,
            b1=R0young/dur,b2=R0old/dur, w=0.97, tinc=300, b1inc=R0young2/dur, b2inc=R0old2/dur ,winc=0.863)
state2  = c(S1=10000-7.5,E1=3.75,I1=3.75,H1=0,C1=0,R1=0,D1=0,
            S2= 7000-2.5,E2=1.25,I2=1.25,H2=0.3,C2=0.2,R2=0,D2=0)
times <- seq(0, 700, by = 1)
out = as.data.frame(ode(y = state2, times = times, func = SEIR2, parms = params2, method="adams"))

# View(out)

png("besmettingen.png")
plot(out$time,out$I1,type="l",col="red",ylim=c(0,250),xlab="dagen",ylab="aantal infecties (duizenden)",lwd=2)
lines(out$time,out$I2,type="l",lty=2,col="magenta",lwd=2)
legend(500, 250, legend=c("leeftijd < 50", "leeftijd >= 50"),
       col=c("red", "magenta"), lty=1:2, cex=0.8)
dev.off()

png("ziekenhuisbedden.png")
plot(out$time,out$H1*1000,type="l",col="red",ylim=c(0,16000),xlab="dagen",ylab="benodigd aantal ziekenhuisbedden")
lines(out$time,out$H2*1000,type="l",col="magenta",lty=2)
lines(out$time,(out$H1+out$H2)*1000,type="l",col="blue",lwd=2)
legend(450, 15000, legend=c("totaal","leeftijd < 50", "leeftijd >= 50"),
       col=c("blue","red", "magenta"), lty=c(1,1,2), cex=0.8)
dev.off()

png("IC-bedden.png")
plot(out$time,out$C1*1000,type="l",col="red",ylim=c(0,4000),xlab="dagen",ylab="benodigd aantal IC bedden")
lines(out$time,out$C2*1000,type="l",col="magenta",lty=2)
lines(out$time,(out$C1+out$C2)*1000,type="l",col="blue",lwd=2)
legend(500, 4000, legend=c("totaal","leeftijd < 50", "leeftijd >= 50"),
       col=c("blue","red", "magenta"), lty=c(1,1,2), cex=0.8)
dev.off()

png("immune.png")
plot(out$time,out$R1,type="l",col="red",ylim=c(0,5500),xlab="dagen",ylab="hersteld en immuun (duizenden)")
lines(out$time,out$R2,type="l",col="magenta",lty=2)
legend(500, 5500, legend=c("leeftijd < 50", "leeftijd >= 50"),
       col=c("red", "magenta"), lty=c(1,1,2), cex=0.8)
dev.off()

