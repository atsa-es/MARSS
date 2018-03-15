###################################################
### code chunk number 3: Cs00_required_packages
###################################################
library(MARSS)
library(xtable)
library(survival)
library(Formula)
library(ggplot2)
library(Hmisc)
library(datasets)
library(broom)


###################################################
### code chunk number 4: Cs01_plotdata
###################################################
#load the datasets package
library(datasets)
data(Nile)   #load the data
plot(Nile,ylab="Flow volume",xlab="Year")


###################################################
### code chunk number 5: Cs02_mod-nile-0
###################################################
mod.nile.0 = list( 
Z=matrix(1), A=matrix(0), R=matrix("r"),
B=matrix(1), U=matrix(0), Q=matrix(0),
x0=matrix("a") )


###################################################
### code chunk number 6: Cs03_fit-data-0
###################################################
#The data is in a ts format, and we need a matrix
dat = t(as.matrix(Nile))
#Now we fit the model
kem.0 = MARSS(dat, model=mod.nile.0)


###################################################
### code chunk number 7: Cs04_mod-nile-1
###################################################
mod.nile.1 = list(
Z=matrix(1), A=matrix(0), R=matrix("r"),
B=matrix(1), U=matrix("u"), Q=matrix(0),
x0=matrix("a") )


###################################################
### code chunk number 8: Cs05_fit-data-1
###################################################
kem.1 = MARSS(dat, model=mod.nile.1, silent=TRUE)
glance(kem.1)


###################################################
### code chunk number 9: Cs06_mod-nile-2
###################################################
mod.nile.2 = list(
Z=matrix(1), A=matrix(0), R=matrix("r"),
B=matrix(1), U=matrix(0), Q=matrix("q"),
x0=matrix("pi") )


###################################################
### code chunk number 11: Cs07_fit-data-2
###################################################
kem.2em = MARSS(dat, model=mod.nile.2, silent=TRUE)
kem.2 = MARSS(dat, model=mod.nile.2,
  inits=kem.2em$par, method="BFGS", silent=TRUE)
tidy(kem.2, conf.int=FALSE)


###################################################
### code chunk number 13: Cs08_plotfit
###################################################
library(Hmisc)
par(mfrow=c(3,1), mar=c(4,4,0.5,0.5), oma=c(1,1,1,1))
x=seq(tsp(Nile)[1],tsp(Nile)[2],tsp(Nile)[3])
#model 0
plot(Nile,ylab="Flow volume",xlab="",xaxp=c(1870,1970,10),bty="L")
minor.tick(nx=10,ny=0,tick.ratio=.3)
kem=kem.0 #model 0 results
lines(x,kem$states[1,],col="red",lwd=2)
legend("topright", paste("model 0, AICc=",format(kem.0$AICc,digits=1)), bty="n")

#model 1
plot(Nile,ylab="Flow volume",xlab="",xaxp=c(1870,1970,10),bty="n")
minor.tick(nx=10,ny=0,tick.ratio=.3)
kem=kem.1 #model 1 results
lines(x,kem$states[1,],col="red",lwd=2)
legend("topright", paste("model 1, AICc=",format(kem.1$AICc,digits=1)),bty="n")

#model 2
plot(Nile,ylab="Flow volume",xlab="",xaxp=c(1870,1970,10),bty="L")
minor.tick(nx=10,ny=0,tick.ratio=.3)
kem=kem.2 #model 0 results
lines(x,kem$states[1,],col="red",lwd=2)
lines(1871:1970,kem$states[1,]-2*kem$states.se[1,],col="red",lty=2)
lines(1871:1970,kem$states[1,]+2*kem$states.se[1,],col="red",lty=2)
legend("topright", paste("model 2, AICc=",format(kem$AICc,digits=1)),bty="n")


###################################################
### code chunk number 14: Cs09_compute-resids
###################################################
resids.0=residuals(kem.0)$std.residuals
resids.1=residuals(kem.1)$std.residuals
resids.2=residuals(kem.2)$std.residuals


###################################################
### code chunk number 15: Cs10_plotoutliertests
###################################################
library(Hmisc)
par(mfrow=c(3,1),mar=c(3,4,1.5,2))
x=seq(tsp(Nile)[1],tsp(Nile)[2],tsp(Nile)[3])
plot(x,resids.0[1,],ylab="std. residuals",xlab="",type="l",
   ylim=c(-4,4),xaxp=c(1870,1970,10),bty="L")
minor.tick(nx=10,ny=0,tick.ratio=.3)
abline(h=c(1.97,-1.97,0),lty=2)
title("model 0--flat level")

plot(x,resids.1[1,],ylab="std. residuals",xlab="",type="l",
   ylim=c(-4,4),xaxp=c(1870,1970,10),bty="L")
minor.tick(nx=10,ny=0,tick.ratio=.3)
abline(h=c(1.97,-1.97,0),lty=2)
title("model 1--linearly declining level")

plot(x,resids.2[1,],ylab="std. residuals",xlab="",type="l",
   ylim=c(-4,4),xaxp=c(1870,1970,10),bty="L")
minor.tick(nx=10,ny=0,tick.ratio=.3)
abline(h=c(1.97,-1.97,0),lty=2)
title("model 2--stochastic level")


###################################################
### code chunk number 16: Cs11_plotresids
###################################################
par(mfrow=c(2,1),mar=c(4,3,2,1))
x=seq(tsp(Nile)[1],tsp(Nile)[2],tsp(Nile)[3])
plot(x,resids.2[1,],ylab="",xlab="",type="l",ylim=c(-4,4),xaxp=c(1870,1970,10))
minor.tick(nx=10,ny=0,tick.ratio=.3)
abline(h=c(1.97,-1.97),lty=2)
title("test for outliers")

plot(x,resids.2[2,],ylab="",xlab="",type="l",ylim=c(-4,4),xaxp=c(1870,1970,10))
minor.tick(nx=10,ny=0,tick.ratio=.3)
abline(h=c(1.97,-1.97),lty=2)
title("test for level changes")
mtext("standardized residuals", side=2, outer=TRUE, line=-1) 


