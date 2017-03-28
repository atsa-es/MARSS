### R code from vignette source 'CIs.Rnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)
options(prompt=" ", continue=" ", width=60)


###################################################
### code chunk number 2: readdat
###################################################
# Data compiled by 
# Dr. Arthur M. Shapiro
# http://butterfly.ucdavis.edu/
library(stringr)
bdat=read.csv("butterfly_data.csv")
sppcol=which(str_detect(colnames(bdat),"Skipper"))
sppcol=which(str_detect(colnames(bdat),"Common.Skipper"))
nspp = length(sppcol)
sppnames=colnames(bdat)[sppcol]
#badyr=bdat$year==1983 | bdat$year==1985
badyr=1983
for(i in sppcol) bdat[,i]=bdat[,i]-mean(bdat[,i])
#express as relative to mean date of first spring flight
ave.max.win = (bdat$ave.max.mar.temp+bdat$ave.max.feb.temp+bdat$ave.max.jan.temp)/3
dat=data.frame(
  y=unlist(bdat[,sppcol]),
  x=rep(ave.max.win,nspp),
  year=rep(bdat$year, each=nspp),
  name=rep(sppnames,each=dim(bdat)[1]))
dat=dat[!(dat$year %in% badyr),]


###################################################
### code chunk number 3: Cs_000_fig0
###################################################
par(mfrow=c(1,2),mar=c(5,5,2,2))
n = dim(dat)[1]
fit = lm(y~x, data=dat)
sigma = sqrt(sum(fit$residual^2)/fit$df.residual) #unbiased not MLE
alpha=coef(fit)[1]
beta=coef(fit)[2]

plot(dat$x, dat$y, xlim=c(min(ave.max.win),max(ave.max.win)), ylim=c(min(bdat[,sppcol]),max(bdat[,sppcol])),
     ylab="Day of first spring flight\nminus the mean", xlab="Average max Jan-Mar temperature (F)",
     bty="L")
abline(fit, col="red", lwd=2)
m=10000
x=runif(m,min(ave.max.win),max(ave.max.win))
y=alpha+beta*x+rnorm(m,0,sigma)
points(x,y,pch=".",col="red")
points(dat$x,dat$y,pch=16,col="blue")
points(rep(ave.max.win[bdat$year %in% badyr],length(sppcol)),bdat[bdat$year %in% badyr,sppcol],pch=1,col="blue")

#panel b
hist(dat$x,xlab="Ave max winter temperature", main="")



###################################################
### code chunk number 4: Cs_001_fit-j
###################################################
set.seed(123)
nsamp=10 #sample size
x.j = runif(nsamp, min(dat$x), max(dat$x))
y.j =  alpha + beta*x.j + rnorm(nsamp,0,sigma)
dat.j = data.frame(x=x.j, y=y.j)
#we will use the fit to the j-data throughout the chapter
fit.j = lm(y~x, data=dat.j)
df.j = fit.j$df.residual
alpha.j = coef(fit.j)[1]
beta.j = coef(fit.j)[2]
sigma.j = sqrt(sum(fit.j$residual^2)/df.j) #unbiased not MLE


###################################################
### code chunk number 5: Cs_001_fit-i
###################################################
nsim = 5000
i.results=matrix(NA,nsim,3)
for(i in 1:nsim){
  x = runif(nsamp, min(dat$x), max(dat$x))
  y = alpha + beta*x + rnorm(nsamp,0,sigma)
  dat.i=data.frame(x=x, y=y)
  fit.i=lm(y~x, data=dat.i)
  i.results[i,]=c(coef(fit.i), 
                  sqrt(sum(fit.i$residual^2)/fit.i$df.residual))
}


###################################################
### code chunk number 6: Cs_001_fig1
###################################################
par(mfrow=c(2,2), mar=c(4,4,2,2))
#truth
n = dim(dat)[1]
fit = lm(y~x, data=dat)
sigma = sqrt(sum(fit$residual^2)/fit$df.residual)
alpha=coef(fit)[1]
beta=coef(fit)[2]

plot(dat$x, dat$y, ylab="first flight relative to mean", xlab="", col="red")
abline(fit, col="red", lwd=2)

#regression line from the j-th sample of 10
points(dat.j$x, dat.j$y, pch=17)
legend("topright","(a)",bty="n")
abline(fit.j)
title("(red) true regression line\n(black) regression line with j sample of 10",cex.main=.75)

#distribution of all those regression lines from a sample of 10
#j sample of 10
nsim = 5000
predx=seq(min(dat$x),max(dat$x),(max(dat$x)-min(dat$x))/7)[2:7]; nx = length(predx)
j.results=matrix(NA,nsim,3)
j.ci99 = j.ci95 = j.ci75 = matrix(NA,nsim,3*nx)
for(i in 1:nsim){
  x = runif(nsamp, min(dat$x), max(dat$x))
  y = alpha + beta*x + rnorm(nsamp,0,sigma)
  tmp.dat=data.frame(x=x, y=y)
  tmp.fit=lm(y~x, data=tmp.dat)
  j.results[i,]=c(coef(tmp.fit), sum(tmp.fit$residual^2)/tmp.fit$df.residual)
  j.ci99[i,] = as.vector(predict(tmp.fit, newdata=data.frame(x=predx), interval="confidence", level=0.99))
  j.ci95[i,] = as.vector(predict(tmp.fit, newdata=data.frame(x=predx), interval="confidence", level=0.95))
  j.ci75[i,] = as.vector(predict(tmp.fit, newdata=data.frame(x=predx), interval="confidence", level=0.75))
}
#plot out 1000 of the regression lines
plot(dat$x, dat$y, type="n", ylab="first flight relative to mean", xlab="")  
title("1000 bootstrapped regression lines\nfrom samples of 10",cex.main=.75)
for(i in 1:1000){ lines(c(min(dat$x),max(dat$x)), j.results[i,1]+j.results[i,2]*c(min(dat$x),max(dat$x))) }
legend("topright","(b)",bty="n")
abline(fit, col="red", lwd=2)

#confidence intervals
plot(dat$x, dat$y, type="n", ylab="first flight relative to mean", xlab="ave. max. temperature")
abline(fit, col="red", lwd=2)
points(dat.j$x, dat.j$y, pch=17)
abline(fit.j)
tmp.cis = predict(fit.j, newdata=data.frame(x=predx), interval="confidence")
for(i in 1:nx){
  x=predx[i]
  lines(c(x,x), tmp.cis[i,2:3])
}
legend("topright","(c)",bty="n")
title("Confidence intervals for the\nj sample of 10",cex.main=.75)

#Coverage
true.y = predict(fit, newdata=data.frame(x=predx))
ci.obs95 = ci.obs75 = ci.obs99 = rep(NA,nx)
for(i in 1:nx){
  ci.obs99[i] = sum(j.ci99[,i+nx]<true.y[i] & j.ci99[,i+2*nx]>true.y[i])/nsim
  ci.obs95[i] = sum(j.ci95[,i+nx]<true.y[i] & j.ci95[,i+2*nx]>true.y[i])/nsim
  ci.obs75[i] = sum(j.ci75[,i+nx]<true.y[i] & j.ci75[,i+2*nx]>true.y[i])/nsim
}

plot(dat$x, dat$y, type="n", ylim=c(-.02,.02), 
     ylab="fraction of coverage\nminus correct fraction", xlab="ave. max. temperature")
abline(h=0,col="red")
#points(predx, ci.obs75)
lines(predx, ci.obs99-0.99, lty=1)
lines(predx, ci.obs95-0.95, lty=2)
lines(predx, ci.obs75-0.75, lty=3)

legend("topright","(d)",bty="n")
legend("bottomright",c("99%","95%","75%"),lty=1:3,bty="n")
title(paste("Coverage of the",nsim, "confidence intervals\nversus objective (red line)"),cex.main=.75)


###################################################
### code chunk number 7: Cs_002_fig_cartoon
###################################################
nsim=1000
ij.results=matrix(NA,nsim,2)
for(i in 1:nsim){
  y.ij=alpha+beta*x.j+rnorm(nsamp,0,sigma)
  fit.ij=lm(y.ij ~ x.j)
  ij.results[i,]=coef(fit.ij)
}

par(mfcol=c(3,2),oma=c(1.5,1.5,3,0))
## Left side of figure
fig.cartoon.base=function(){
  plot(dat$x,dat$y, ylim=c(-30,30), type="n", ylab="", xlab="", xlim=c(56,63))
  for(i in 1:nsim){
    lines(c(56,63), ij.results[i,1]+ij.results[i,2]*c(56,63),col="grey")
  }
  abline(v=x.j,col="green")
  abline(fit, col="red", lwd=2)
}

#panel a
fig.cartoon.base()
title("Regression lines from samples of 10\nwith x values at green lines",cex.main=.75)
legend("topright","(a)",bty="n")

#panel b
fig.cartoon.base()
x1=58
lines( c(x1,x1), quantile(ij.results[,1]+ij.results[,2]*x1,probs=c(0.025,.975)), lwd=4, col="blue" )
title("blue line contains 95% of the regression lines\nfrom samples at the green x",cex.main=.75)
legend("topright","(b)",bty="n")

#panel c
plot(dat$x,dat$y, ylim=c(-30,30), type="n", ylab="", xlab="", xlim=c(56,63))
abline(v=x.j, col="green")
abline(fit, col="red", lwd=2)
points(x.j, y.j, pch=3)
abline(fit.j, col="black", lwd=2)
x1=58
lines( c(x1,x1), (coef(fit.j)[1]-alpha+(coef(fit.j)[2]-beta)*x1)+quantile(ij.results[,1]+ij.results[,2]*x1,probs=c(0.025,.975)), lwd=4, col="blue" )
title("the blue line centered on any regression line from\na sample of 10 at green x covers red line 95% of the time",cex.main=.75)
legend("topright","(c)",bty="n")

# #panel d
# plot(dat$x,dat$y, ylim=c(-30,30), type="n", ylab="", xlab="", xlim=c(56,63))
# set.seed(222)
# x.k = runif(nsamp, min(dat$x), max(dat$x))
# y.k = alpha + beta*x.k + rnorm(nsamp,0,sqrt(sigma2))
# dat.k=data.frame(x=x.k, y=y.k)
# fit.k = lm(y~x, data=dat.k)
# ik.results=matrix(NA,100,2)
# for(i in 1:100){
#   y.ik=alpha+beta*x.k+rnorm(10,0,sqrt(sigma2))
#   fit.ik=lm(y.ik ~ x.k)
#   ik.results[i,]=coef(fit.ik)
# }
# abline(fit.k, col="black", lwd=2)
# points(x.k, y.k)
# abline(fit, col="red", lwd=2)
# abline(v=x.k, col="green")
# x1=58
# lines( c(x1,x1), (coef(fit.k)[1]-alpha+(coef(fit.k)[2]-beta)*x1)+quantile(ik.results[,1]+ik.results[,2]*x1,probs=c(0.025,.975)), lwd=4, col="blue" )
# title("a 95% CI constructed in this way from samples\nwith other x values will also cover\nthe red line 95% of the time",cex.main=.75)

ik.results=matrix(NA,nsim,2)
for(i in 1:nsim){
  x.k = sample(dat$x,nsamp,replace=TRUE)
  y.ik=alpha+beta*x.k+rnorm(nsamp,0,sigma)
  fit.ik=lm(y.ik ~ x.k)
  ik.results[i,]=coef(fit.ik)
}

## Right side of figure
fig.cartoon.base=function(){
  plot(dat$x,dat$y, ylim=c(-30,30), type="n", ylab="", xlab="", xlim=c(56,63))
  for(i in 1:nsim){
    lines(c(56,63),ij.results[i,1]+ij.results[i,2]*c(56,63),col="grey")
  }
  abline(fit, col="red", lwd=2)
}

#panel d
fig.cartoon.base()
title("Regression lines from samples of 10\nwith random x values",cex.main=.75)
legend("topright","(d)",bty="n")

#panel b
fig.cartoon.base()
x1=58
lines( c(x1,x1), quantile(ik.results[,1]+ik.results[,2]*x1,probs=c(0.025,.975)), lwd=4, col="blue" )
title("blue line contains 95% of the regression lines\nfrom samples with x drawn from historical range",cex.main=.75)
legend("topright","(e)",bty="n")

#panel c
plot(dat$x,dat$y, ylim=c(-30,30), type="n", ylab="", xlab="", xlim=c(56,63))
abline(fit, col="red", lwd=2)
points(x.j, y.j, pch=3)
abline(fit.j, col="black", lwd=2)
x1=58
lines( c(x1,x1), (coef(fit.j)[1]-alpha+(coef(fit.j)[2]-beta)*x1)+quantile(ik.results[,1]+ik.results[,2]*x1,probs=c(0.025,.975)), lwd=4, col="blue" )
title("the blue line centered on any regression line from\na sample of 10 will cover red line 95% of the time",cex.main=.75)
legend("topright","(f)",bty="n")

mtext(side=3, outer=TRUE, "Two approaches to constructing CIs",line=1)
mtext(side=2, outer=TRUE, "first flight day minus average",line=-1.5)
mtext(side=1, outer=TRUE, "ave. max. temperature",line=-1.5)


###################################################
### code chunk number 8: Cs_001_dist_of_param
###################################################
par(mfrow=c(1,2), mar=c(5,5,4,2))
nsim=5000
nsamp=nrow(dat.j)
true.dist=matrix(NA,nsim,2)
for(i in 1:nsim){
  y.ij=alpha+beta*x.j+rnorm(nsamp,0,sigma)
  fit.ij=lm(y.ij ~ x.j)
  true.dist[i,]=coef(fit.ij)
}

data.dist=matrix(NA,nsim,2)
for(i in 1:nsim){
  y.ij=alpha.j+beta.j*x.j+rnorm(nsamp,0,sigma.j)
  fit.ij=lm(y.ij ~ x.j)
  data.dist[i,]=coef(fit.ij)
}

nbreaks=50
# now we do a kernel density estimate
true.kde = kde2d(true.dist[,1], true.dist[,2], n = nbreaks)
data.kde = kde2d(data.dist[,1], data.dist[,2], n = nbreaks)

cont.lev = c(.002,.001,.0001)
contour(true.kde,levels=cont.lev,xlim=c(0,700),ylim=c(-12,0))
title(
  xlab=expression(paste("intercept (",alpha,")",sep="")),
  ylab=expression(paste("slope (",beta,")",sep=""))
)
title("True distribution")

contour(data.kde, col="red",levels=cont.lev,xlim=c(0,700),ylim=c(-12,0))
title(
  xlab=expression(paste("intercept (",alpha,")",sep="")),
  ylab=expression(paste("slope (",beta,")",sep=""))
)
title("Distribution estimated\nfrom one data set")



###################################################
### code chunk number 9: Cs_0031
###################################################
nsim=5000
ij.results=matrix(NA,nsim,2)
for(i in 1:nsim){
  y.ij=alpha+beta*x.j+rnorm(nsamp,0,sigma)
  fit.ij=lm(y.ij ~ x.j)
  ij.results[i,]=coef(fit.ij)
}
x1=58
CI=quantile(ij.results[,1]+ij.results[,2]*x1,probs=c(0.025,.975))


###################################################
### code chunk number 10: Cs_004_analytical.mle.var
###################################################
Xj=cbind(alpha=1,beta=x.j)
analytical.Sigma=sigma^2*solve(t(Xj)%*%Xj)
analytical.Sigma


###################################################
### code chunk number 11: Cs_007
###################################################
Xj = cbind(1, dat.j$x)
XjXj.inv = solve(t(Xj)%*%Xj)
Sigma = sigma^2*XjXj.inv
theta = rbind(alpha, beta)
x=58; X = cbind(1, x)

#the analytical CI
X%*%theta + qnorm(c(0.025,.975)) * sqrt(X%*%Sigma%*%t(X))


###################################################
### code chunk number 12: Cs_0071
###################################################
quantile(ij.results[,1]+ij.results[,2]*x, probs=c(0.025, 0.975))


###################################################
### code chunk number 13: Cs_0072
###################################################
  Xj = cbind(1, x.j)
  XjXj.inv = solve(t(Xj)%*%Xj)
  Sigma.j = sigma.j^2*XjXj.inv
  theta.j = matrix(c(alpha.j,beta.j), ncol=1)
 #Compute CI at x=58
x=58; X = cbind(1, x)
 EyX = X%*%theta.j
  corrected.ci.j = EyX + qt(c(0.025,.975), df=fit.df) * sqrt(X%*%Sigma.j%*%t(X))
corrected.ci.j


###################################################
### code chunk number 14: Cs_013_predict.lm.con
###################################################
corrected.ci.j
predict(lm(y~x, data=dat.j), new=data.frame(x=x),interval="confidence")


###################################################
### code chunk number 15: Cs_008
###################################################
Xj = cbind(1, x.j)
XjXj.inv = solve(t(Xj)%*%Xj)
true.Sigma = sigma^2*XjXj.inv
#we are going to compute CI at x=58
x=58; X = cbind(1, x)
#holders
nsim=5000
i.CIs.bad=i.CIs.true=i.CIs.corr=matrix(NA,nsim,2)
for(i in 1:nsim){
  tilde.y=alpha+beta*x.j+rnorm(nsamp,0,sigma)
  fit.i=lm(tilde.y ~ x.j)
  hat.theta=matrix(coef(fit.i),ncol=1)
  hat.sigma = sqrt(sum(fit.i$residual^2)/fit.i$df.residual)
  hat.Sigma = hat.sigma^2*XjXj.inv
  meanCI = X%*%hat.theta
  normaldist = qnorm(c(0.025,0.975))
  tdist = qt(c(0.025,0.975),df=fit.i$df.residual)
  #CI using asymptotic equation and estimated Sigma
  i.CIs.bad[i,] = meanCI + normaldist * sqrt(X%*%hat.Sigma%*%t(X))
  #CI using asymptotic equation and true Sigma
  i.CIs.true[i,] = meanCI + normaldist * sqrt(X%*%true.Sigma%*%t(X))
  #Corrected CI using t-distribution
  i.CIs.corr[i,] = meanCI + tdist * sqrt(X%*%hat.Sigma%*%t(X))
}


###################################################
### code chunk number 16: Cs_009
###################################################
x=58
true.val = alpha+beta*x
100*sum(i.CIs.true[,1]<true.val & i.CIs.true[,2]>true.val)/nsim


###################################################
### code chunk number 17: Cs_010
###################################################
100*sum(i.CIs.bad[,1]<true.val & i.CIs.bad[,2]>true.val)/nsim


###################################################
### code chunk number 18: Cs_011
###################################################
100*sum(i.CIs.corr[,1]<true.val & i.CIs.corr[,2]>true.val)/nsim


###################################################
### code chunk number 19: Cs_020_linear.reg.LL
###################################################
# Define the log likelihood function for a linear regression
# parm is the alpha, beta, sigma vector
NLL <- function(parm,  dat=NULL){
  #parm is alpha, beta, sigma
  resids = dat$y - dat$x * parm[2] - parm[1]
  dresids = suppressWarnings(dnorm(resids, 0, parm[3], log = TRUE))
  -sum(dresids)
}


###################################################
### code chunk number 20: Cs_023_linear.ref.optim
###################################################
start.pars = c(alpha.j, beta.j, sigma.j)
ofit.j=optim(start.pars, NLL, dat=dat.j, hessian=TRUE)
parSigma = solve(ofit.j$hessian)[1:2,1:2]
parMean = ofit.j$par[1:2]


###################################################
### code chunk number 21: Cs_023_hessian_cis
###################################################
x1=58
#add 0 to deal with sigma in parMean
X = cbind(1, x1)
EyX = X%*%parMean
hessian.cis = c(EyX, EyX + qnorm(c(0.025,.975)) * sqrt(X%*%parSigma%*%t(X)))


###################################################
### code chunk number 22: Cs_024_hessian_cis
###################################################
rbind(hessian=hessian.cis, correct=correct.ci.j) 


###################################################
### code chunk number 23: Cs_014_boot_CIs_fun
###################################################
#takes a set of boot parameters and makes CIs from them at x
boot.CI=function(boot.params, x, alp=0.05){
  CIs=apply(
    boot.params[,c("alpha","beta"),drop=FALSE]%*%rbind(1,x),
    2,quantile, c(0.5, alp/2, 1-alp/2)
  )
  colnames(CIs)=x
  t(CIs) #to look like predict output
}


###################################################
### code chunk number 24: Cs_015_parametricboot
###################################################
parametric.boot=function(dat, nboot=1000){
  #first fit model to data
  fit=lm(y~x, data=dat)
  #x's at which to generate new data
  x=dat$x
  #matrix to store the estimates
  boot.params=matrix(NA,nboot,3)
  sigma = sqrt(sum(fit$residual^2)/fit$df.residual)
  alpha=coef(fit)[1]
  beta=coef(fit)[2]
  for(i in 1:nboot){
    y=alpha + beta*x + rnorm(nrow(dat),0,sigma)
    tmp.fit=lm(y~x)
    boot.params[i,]=c(coef(tmp.fit), 
                      sqrt(sum(tmp.fit$residual^2)/tmp.fit$df.residual))
  }
  colnames(boot.params)=c("alpha","beta","sigma")
  boot.params
}

parametric.boot.cis=function(dat, x, nboot=1000){
  boot.params=parametric.boot(dat, nboot=nboot)
  boot.CI(boot.params, x)
}


###################################################
### code chunk number 25: Cs_016_parametricboot.cis
###################################################
rbind(
  parametric=parametric.boot.cis(dat.j, x=58)[1,],
  correct.ci.j
  )


###################################################
### code chunk number 26: Cs_021_fig_biv_comparison
###################################################
par(mfrow=c(1,1), mar=c(5,5,2,2))
nbreaks=50
#generate bivariate normal alpha and beta from the estimated distribution
hessian.par = mvrnorm(1000, mu = parMean, Sigma = parSigma)
# now we do a kernel density estimate
hessian.kde = kde2d(hessian.par[,1], hessian.par[,2], n = nbreaks)

# now we do a kernel density estimate using the bootstrapped parameters
boot.params=parametric.boot(dat.j)
boot.kde = kde2d(boot.params[,1], boot.params[,2], n = nbreaks)
cont.lev = c(.002,.001,.0001)
contour(boot.kde,levels=cont.lev)
contour(hessian.kde, add=TRUE, col="red",levels=cont.lev)
title(
  xlab=expression(paste("intercept (",alpha,")",sep="")),
  ylab=expression(paste("slope (",beta,")",sep=""))
)
legend("topright",c("actual distribution","distribution from Hessian"),lty=1,col=c("black","red"),bty="n")


###################################################
### code chunk number 27: Cs_017_residuals_boot_function
###################################################
residuals.boot=function(dat, nboot=1000){
  fit=lm(y~x, data=dat)
  resids=residuals(fit)  
  alpha=coef(fit)[1]
  beta=coef(fit)[2]
  boot.params=matrix(NA,nboot,3)
  n = nrow(dat) #number of data points
  for(i in 1:nboot){
    tmp = sample(n, replace=TRUE)
    tmp.y=alpha + beta*dat$x + resids[tmp]
    tmp.fit=lm(tmp.y~dat$x)
    boot.params[i,]=c(coef(tmp.fit), 
                      sqrt(sum(tmp.fit$residual^2)/tmp.fit$df.residual))
  }
  colnames(boot.params)=c("alpha","beta","sigma")
  boot.params
}

residuals.boot.cis=function(dat, x, nboot=1000){
  boot.params=residuals.boot(dat, nboot=nboot)
  boot.CI(boot.params, x)
}


###################################################
### code chunk number 28: Cs_018_residuals_boot_cis
###################################################
rbind(
  residuals=residuals.boot.cis(dat.j, x=58)[1,],
  correct.ci.j
)


###################################################
### code chunk number 29: Cs_024_non_param_boot_function
###################################################
resampling.boot=function(dat, nboot=1000){
  boot.params=matrix(NA,nboot,3)
  n = nrow(dat) #number of data points
  for(i in 1:nboot){
    #sample with replacement
    tmp = sample(n, replace=TRUE)
    #tmp.fit is the fit to this bootstrapped data
    tmp.fit=lm(y~x, data=dat, subset=tmp)
    boot.params[i,]=c(
      coef(tmp.fit),
      sqrt(sum(tmp.fit$residuals^2)/tmp.fit$df.residual))
  }
  colnames(boot.params)=c("alpha","beta","sigma")
  boot.params
}

resampling.boot.cis=function(dat, x, nboot=5000){
  boot.params=resampling.boot(dat, nboot=nboot)
  boot.CI(boot.params, x)
}


###################################################
### code chunk number 30: Cs_025_resampling_cis
###################################################
rbind(
  resampling=resampling.boot.cis(dat.j, x=58)[1,],
  correct=correct.ci.j
)


###################################################
### code chunk number 31: Cs_026_resampling_cis_big
###################################################
#create a larger sample
nsamp=50
x = runif(nsamp, min(dat$x), max(dat$x))
y = alpha + beta*x + rnorm(nsamp,0,sigma)
dat.big=data.frame(x=x, y=y)

rbind(
  residuals=resampling.boot.cis(dat.big, x=x1)[1,],
  correct=predict(lm(y~x, data=dat.big), new=data.frame(x=58),interval="conf")
)


###################################################
### code chunk number 32: Cs_025_ci_adj
###################################################
#the correct adjustment:
qt(c(0.025),df=fit.j$df.residual)/qnorm(c(0.025))

#adjustment computed for the hessian bootstrap CI
x1=58
ci.adj(dat.j,x1)


###################################################
### code chunk number 33: Cs_5
###################################################
fit = lm(mpg~wt, data=mtcars)
s2 = sum(fit$residual^2)/fit$df.residual
alpha=coef(fit)[1]
beta=coef(fit)[2]
plot(mtcars$wt,mtcars$mpg,ylim=c(0,60),xlab="car weight",ylab="mpg")
x = runif(1000,1,6)
y = alpha+beta*x+rnorm(1000,0,sqrt(s2))
points(x,y)
abline(fit, lwd=2, col="red")


###################################################
### code chunk number 34: Cs_6
###################################################
lines(c(0,10),alpha+beta*c(0,10)+1.96*sqrt(s2),col="red",lty=2,lwd=1)
lines(c(0,10),alpha+beta*c(0,10)-1.96*sqrt(s2),col="red",lty=2,lwd=1)


###################################################
### code chunk number 35: Cs_7 (eval = FALSE)
###################################################
## npred=1000
## plot(mtcars$wt,mtcars$mpg,ylim=c(0,60),type="n",xlab="car weight",ylab="mpg")
## x = runif(npred,1,6)
## y = alpha+beta*x+rnorm(npred,0,sqrt(s2))
## points(x,y, col="grey")
## abline(fit, lwd=2, col="red")
## abline(fit.j, lwd=2)
## points(mtcars.j$wt, mtcars.j$mpg, pch=3)
## preds = predict(lm(mpg~wt,data=mtcars.j), newdata= data.frame(wt=pred.wt), interval="prediction")
## lines(pred.wt,preds[,2],lty=2)
## lines(pred.wt,preds[,3],lty=2)


###################################################
### code chunk number 36: Cs_8 (eval = FALSE)
###################################################
## #get the 95% pred intervals for each x
## preds = predict(fit.j, newdata= data.frame(wt=x), interval="prediction")
## #see how many y fall outside that
## 1-sum(y>preds[,3] | y<preds[,2])/npred
## #see how many fall inside at different x values
## 1-tapply(y>preds[,3] | y<preds[,2],cut(x,breaks=1:6),mean)


###################################################
### code chunk number 37: Cs_9 (eval = FALSE)
###################################################
## #j sample of 10 cars
## nsim = 5000
## pi.coverage=rep(NA, nsim)
## for(i in 1:nsim){
##   tmp.fit=lm(mpg~wt, data=mtcars, subset=sample(dim(mtcars)[1], 10))
##   preds = predict(tmp.fit, newdata= data.frame(wt=x), interval="prediction")
##   pi.coverage[i] = 1-sum(y>preds[,3] | y<preds[,2])/npred
## }
## sum(pi.coverage)/nsim


###################################################
### code chunk number 38: Cs_10 (eval = FALSE)
###################################################
## y = alpha+beta*x+sample(residuals(fit), npred, replace=TRUE)
## nsim = 5000
## pi.coverage=rep(NA, nsim)
## for(i in 1:nsim){
##   tmp.fit=lm(mpg~wt, data=mtcars, subset=sample(dim(mtcars)[1], 10))
##   preds = predict(tmp.fit, newdata= data.frame(wt=x), interval="prediction")
##   pi.coverage[i] = 1-sum(y>preds[,3] | y<preds[,2])/npred
## }
## sum(pi.coverage)/nsim


###################################################
### code chunk number 39: fig-PIs (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## set.seed(32421)
## #TRUTH (the relationship in the 32-car dataset
## fit = lm(mpg~wt, data=mtcars)
## s2 = sum(fit$residual^2)/fit$df.residual
## alpha=coef(fit)[1]
## beta=coef(fit)[2]
## x = runif(1000,1,6)
## y = alpha+beta*x+rnorm(1000,0,sqrt(s2))
## ylims=c(0,60); xlims=c(1,6)
## 
## #panel 1
## #using MLEs
## plot(mtcars$wt,mtcars$mpg,ylim=ylims,xlim=xlims,xlab="car weight",ylab="mpg",type="n")
## points(x,y, col="grey")
## abline(fit.j, lwd=2, col="red")
## points(mtcars.j$wt, mtcars.j$mpg, pch=3,col="red")
## lines(c(0,10),alpha.j+beta.j*c(0,10)+1.96*sqrt(s2.j),col="red",lty=2,lwd=1)
## lines(c(0,10),alpha.j+beta.j*c(0,10)-1.96*sqrt(s2.j),col="red",lty=2,lwd=1)
## title("using MLEs (a)",cex.main=1)
## #coverage using MLEs
## nsim = 1000
## pi.coverage=rep(NA, nsim)
## for(i in 1:nsim){
##   npred=length(x)
##   tmp.fit=lm(mpg~wt, data=mtcars, subset=sample(dim(mtcars)[1], 10))
##   s2.tmp = sum(tmp.fit$residual^2)/fit.j$df.residual
##   alpha.tmp=coef(tmp.fit)[1]
##   beta.tmp=coef(tmp.fit)[2]
##   preds = cbind(
##     alpha.tmp+x*beta.tmp+1.96*sqrt(s2.tmp),
##     alpha.tmp+x*beta.tmp-1.96*sqrt(s2.tmp)
##   )
##   pi.coverage[i] = 1-sum(y>preds[,1] | y<preds[,2])/npred
## }
## legend("top",legend=paste("coverage",round(100*sum(pi.coverage)/nsim,digits=1)),bty="n")
## 
## #panel 2
## #analytical
## plot(mtcars$wt,mtcars$mpg,ylim=ylims,xlim=xlims,xlab="car weight",ylab="mpg",type="n")
## points(x,y, col="grey")
## abline(fit.j, lwd=2,col="red")
## points(mtcars.j$wt, mtcars.j$mpg, pch=3,col="red")
## preds = predict(lm(mpg~wt,data=mtcars.j), newdata= data.frame(wt=pred.wt), interval="prediction")
## lines(pred.wt,preds[,2],lty=2,col="red")
## lines(pred.wt,preds[,3],lty=2,col="red")
## title("analytical (b)",cex.main=1)
## #coverage using analytical
## nsim = 1000
## pi.coverage=rep(NA, nsim)
## for(i in 1:nsim){
##   npred=length(x)
##   tmp.fit=lm(mpg~wt, data=mtcars, subset=sample(dim(mtcars)[1], 10))
##   preds = predict(tmp.fit, newdata= data.frame(wt=x), interval="prediction")
##   pi.coverage[i] = 1-sum(y>preds[,3] | y<preds[,2])/npred
## }
## legend("top",legend=paste("coverage",round(100*sum(pi.coverage)/nsim,digits=1)),bty="n")
## 
## #panel 3
## #parametric bootstrap
## #function to do get PIs from parametric bootstrap
## par.boot.PIs=function(predvar, pars, PIx, nboot=1000){
##   par.results=matrix(NA,nboot,3)
##   for(i in 1:nboot){
##     tmp=pars[1] + pars[2]*predvar + rnorm(length(predvar),0,sqrt(pars[3]))
##     tmp.fit=lm(tmp~predvar)
##     par.results[i,]=c(coef(tmp.fit),sum(tmp.fit$residual^2)/tmp.fit$df.residual)
##   }
##   res=par.results
##   PIs=apply(matrix(res[,1],length(PIx),dim(res)[1],byrow=TRUE)+
##               matrix(PIx,ncol=1)%*%matrix(res[,2],nrow=1)+
##               matrix(rnorm(dim(res)[1],0,sqrt(res[,3])),length(PIx),dim(res)[1]),
##             1,quantile,probs=c(0.025,0.975))
##   
##   PIs
## }
## plot(mtcars$wt,mtcars$mpg,ylim=ylims,xlim=xlims,xlab="car weight",ylab="mpg",type="n")
## points(x,y, col="grey")
## abline(fit.j, lwd=2,col="red")
## points(mtcars.j$wt, mtcars.j$mpg, pch=3,col="red")
## PIs=par.boot.PIs(mtcars.j$wt, c(alpha.j, beta.j, s2.j), pred.wt)
## lines(pred.wt, PIs[1,],lty=2,col="red")
## lines(pred.wt, PIs[2,],lty=2,col="red")
## title("parametric bootstrap (c)",cex.main=1)
## #coverage using par boot.  Just use x at 1:6
## nsim = 100
## ys=c()
## for(i in 1:6) ys=cbind(ys,alpha+beta*i+rnorm(1000,0,sqrt(sigma2)))
## pi.coverage=rep(NA, nsim)
## for(i in 1:nsim){
##   tmp.samp=sample(dim(mtcars)[1], 10)
##   tmp.fit=lm(mpg~wt, data=mtcars, subset=tmp.samp)
##   s2.tmp = sum(tmp.fit$residual^2)/fit.j$df.residual
##   alpha.tmp=coef(tmp.fit)[1]
##   beta.tmp=coef(tmp.fit)[2]
##   PIs=par.boot.PIs(mtcars$wt[tmp.samp], c(alpha.tmp, beta.tmp, s2.tmp), 1:6, nboot=1000)
##   pi.coverage[i] = 1-sum(ys>matrix(PIs[2,],1000,6,byrow=TRUE) | ys<matrix(PIs[1,],1000,6,byrow=TRUE))/6000
## }
## legend("top",legend=paste("coverage",round(100*sum(pi.coverage)/nsim,digits=1)),bty="n")
## 
## 
## #panel 4
## #bootstrap using hessian
## hess.boot.PIs=function(predvar, respvar, PIx, nboot=1000){
##   LL <- function(parm, y=NULL, x=NULL) {
##     resids = y - x * parm[2] - parm[1]
##     resids = suppressWarnings(dnorm(resids, 0, sqrt(parm[3]), log = TRUE))
##     -sum(resids)
##   }
##   tmp.fit=lm(respvar~predvar)
##   pars=c(coef(tmp.fit), sum(tmp.fit$residual^2)/tmp.fit$df.residual)
##   ofit.j=optim(pars, LL, y=respvar, x=predvar, hessian=TRUE)
##   parSigma = solve(ofit.j$hessian)
##   parMean = ofit.j$par
##   hess.results = mvrnorm(nboot, mu = parMean, Sigma = parSigma)
##   res=hess.results[hess.results[,3]>0,]
##   PIs=apply(matrix(res[,1],length(PIx),dim(res)[1],byrow=TRUE)+
##               matrix(PIx,ncol=1)%*%matrix(res[,2],nrow=1)+
##               matrix(rnorm(dim(res)[1],0,sqrt(res[,3])),length(PIx),dim(res)[1]),
##             1,quantile,probs=c(0.025,0.975))
##   PIs
## }
## plot(mtcars$wt,mtcars$mpg,ylim=ylims,xlim=xlims,xlab="car weight",ylab="mpg",type="n")
## points(x,y, col="grey")
## abline(fit.j, lwd=2,col="red")
## points(mtcars.j$wt, mtcars.j$mpg, pch=3,col="red")
## PIs=hess.boot.PIs(mtcars.j$wt, mtcars.j$mpg, pred.wt, nboot=1000)
## lines(pred.wt, PIs[1,],lty=2,col="red")
## lines(pred.wt, PIs[2,],lty=2,col="red")
## title("bootstrap using hessian (d)",cex.main=1)
## #coverage using hess boot.  Just use x at 1:6
## nsim = 100
## ys=c()
## for(i in 1:6) ys=cbind(ys,alpha+beta*i+rnorm(1000,0,sqrt(sigma2)))
## pi.coverage=rep(NA, nsim)
## for(i in 1:nsim){
##   tmp.samp=sample(dim(mtcars)[1], 10)
##   tmp.fit=lm(mpg~wt, data=mtcars, subset=tmp.samp)
##   s2.tmp = sum(tmp.fit$residual^2)/fit.j$df.residual
##   alpha.tmp=coef(tmp.fit)[1]
##   beta.tmp=coef(tmp.fit)[2]
##   PIs=hess.boot.PIs(mtcars$wt[tmp.samp], mtcars$mpg[tmp.samp], 1:6, nboot=1000)
##   pi.coverage[i] = 1-sum(ys>matrix(PIs[2,],1000,6,byrow=TRUE) | ys<matrix(PIs[1,],1000,6,byrow=TRUE))/6000
## }
## legend("top",legend=paste("coverage",round(100*sum(pi.coverage)/nsim,digits=1)),bty="n")


###################################################
### code chunk number 40: Cs_9
###################################################
#show example


###################################################
### code chunk number 41: Cs_005_ret.hessian.s2 (eval = FALSE)
###################################################
## #Here's a function to return the Hessian when sigma is known
## ret.hessian.known.s2=function(dat, sigma){
##   NLL <- function(parm, dat=NULL) {
##     resids = dat$y - parm[1] - dat$x * parm[2]
##     resids = suppressWarnings(dnorm(resids, 0, sigma, log = TRUE))
##     -sum(resids)
##   }
##   pars=coef(lm(y ~ x, data=dat))
##   names(pars)=c("alpha","beta")
##   ofit=optim(pars, NLL, dat=dat, hessian=TRUE)
##   parSigma = solve(ofit$hessian)
##   parMean = ofit$par
##   return(list(parMean=parMean, parSigma=parSigma))
## }


###################################################
### code chunk number 42: Cs_019_fig3 (eval = FALSE)
###################################################
## #create many parameter estimates
## nsim=5000
## #ij.results is alpha, beta, sigma from parametric boot of dat.j
## ij.results=matrix(NA,nsim,3)
## for(i in 1:nsim){
##   y.i = alpha.j + beta.j*x.j + rnorm(nsamp,0,sigma.j)
##   fit.i=lm(y.i~x.j)
##   ij.results[i,]=c(coef(fit.i), 
##                   sqrt(sum(fit.i$residual^2)/fit.i$df.residual))
## }
## 
## library(RColorBrewer)
## rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
## r = rf(32)
## 
## nbreaks=50
## library(MASS)
## df = data.frame(x=ij.results[,1],y=ij.results[,2])
## h1 = hist(df$x, breaks=nbreaks, plot=FALSE)
## h2 = hist(df$y, breaks=nbreaks, plot=FALSE)
## top = max(h1$counts, h2$counts)
## k = kde2d(df$x, df$y, n=nbreaks)
## 
## # margins
## oldpar <- par()
## par(mar=c(5,5,1,1))
## layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
## image(k, col=r) #plot the image
## title(
##   xlab=expression(paste("intercept (",alpha,")",sep="")),
##   ylab=expression(paste("slope (",beta,")",sep=""))
## )
## par(mar=c(0,2,1,0))
## barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='red')
## par(mar=c(2,0,0.5,1))
## barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='red', horiz=T)


###################################################
### code chunk number 43: Cs_100
###################################################
CI.analytical.true.s2=function(dat, x, sigma, alp=0.05){
  fit=lm(y~x, data=dat)
  Xi = cbind(1, dat$x)
  XiXi.inv = solve(t(Xi)%*%Xi)
  theta = matrix(coef(fit), ncol=1)
  if(is.null(sigma)){
    sigma = sqrt(sum(fit$residual^2)/fit.df)
  }
  Sigma = sigma^2*XiXi.inv
  X = cbind(1, x)
  EyX = X%*%theta
  #the analytical CI
  CI = EyX + qnorm(c(alp/2,1-alp/2)) * sqrt(X%*%Sigma%*%t(X))
  c(fit=EyX, lwr=CI[1], upr=CI[2])
}

This function computes the analytical CIs with the observed Fisher information matrix using the estimated $\sigma$. This is biased.


###################################################
### code chunk number 44: Cs_100
###################################################
CI.analytical.est.s2=function(dat, x, alp=0.05){
  fit=lm(y~x, data=dat)
  Xi = cbind(1, dat$x)
  XiXi.inv = solve(t(Xi)%*%Xi)
  theta = matrix(coef(fit), ncol=1)
  sigma = sqrt(sum(fit$residual^2)/fit.df)
  Sigma = sigma^2*XiXi.inv
  X = cbind(1, x)
  EyX = X%*%theta
  #the analytical CI
  CI = EyX + qnorm(c(alp/2,1-alp/2)) * sqrt(X%*%Sigma%*%t(X))
  c(fit=EyX, lwr=CI[1], upr=CI[2])
}


###################################################
### code chunk number 45: Cs_100
###################################################
CI.analytical.corrected=function(dat, x, alp=0.05){
  fit=lm(y~x, data=dat)
  fit.df = fit$df.residual
  Xi = cbind(1, dat$x)
  XiXi.inv = solve(t(Xi)%*%Xi)
  hat.sigma = sqrt(sum(fit$residual^2)/fit.df)
  hat.Sigma = hat.sigma^2*XiXi.inv
  hat.theta = matrix(coef(fit), ncol=1)
  X = cbind(1, x)
  EyX = X%*%hat.theta
  CI = EyX + qt(c(alp/2,1-alp/2), df=fit.df) * sqrt(X%*%hat.Sigma%*%t(X))
  c(fit=EyX, lwr=CI[1], upr=CI[2])
}


###################################################
### code chunk number 46: Cs_101
###################################################
#Get the MLEs and MLE Sigma
hessian.parm=function(dat){
  library(MASS)
  NLL <- function(parm, y=NULL, x=NULL) {
    resids = y - x * parm[2] - parm[1]
    dresids = suppressWarnings(dnorm(resids, 0, parm[3], log = TRUE))
    -sum(dresids)
  }
  fit=lm(y~x, data=dat)
  sigma = sqrt(sum(fit$residual^2)/fit$df.residual)
  alpha=coef(fit)[1]
  beta=coef(fit)[2]
  pars=c(alpha, beta, sigma)
  names(pars)=c("alpha","beta","sigma")
  
  fit.tmp=optim(pars, NLL, y=dat$y, x=dat$x, hessian=TRUE)
  parSigma = solve(fit.tmp$hessian)
  parMean = matrix(fit.tmp$par, ncol=1)
  names(parMean)=c("alpha","beta","sigma")
  
  list(parMean=parMean, parSigma=parSigma)
}


###################################################
### code chunk number 47: Cs_102
###################################################
hessian.boot=function(dat, nboot=1000){
  hes=hessian.parm(dat)  
  #generate alpah and beta from Sigma
  boot.params = mvrnorm(nboot, mu = hes$parMean, Sigma = hes$parSigma)
  colnames(boot.params)=c("alpha","beta","sigma")
  boot.params
}


###################################################
### code chunk number 48: Cs_103
###################################################
parametric.boot


###################################################
### code chunk number 49: Cs_103
###################################################
residuals.boot


###################################################
### code chunk number 50: Cs_103
###################################################
resampling.boot


###################################################
### code chunk number 51: CIs.Rnw:1272-1273
###################################################
boot.CI


###################################################
### code chunk number 52: Cs_103
###################################################
hessian.boot.cis=function(dat, x, nboot=1000){
  boot.params=hessian.boot(dat, nboot=nboot)
  boot.CI(boot.params, x)
}

parametric.boot.cis

residuals.boot.cis

resampling.boot.cis


###################################################
### code chunk number 53: Cs_027_ci_correction
###################################################
ci.adj=function(dat, x, nboot1=5000, nboot2=1000, type="hessian"){
  #x is the x where the ci is computed; one value
  #type is analytical, hessian, parametric or residuals
  nsamp=nrow(dat)
  fit=lm(y~x, data=dat)
  sigma = sqrt(sum(fit$residuals^2)/fit$df.residual)
  alpha=coef(fit)[1]
  beta=coef(fit)[2]
  #this is the correct width of the 95\% CI
  ci=parametric.boot.cis(dat,x,nboot=nboot1)[2:3]
  ci.width=ci[2]-ci[1]
  #ci.width.boot will hold a set of bootstrapped CIs
  ci.width.boot=rep(NA,nboot2)
  for(i in 1:nboot2){
    if(type %in% c("parametric", "hessian")){
      y=alpha + beta*dat$x + rnorm(nsamp,0,sigma)
      tmp.dat=data.frame(x=dat$x, y=y)
      if(type=="analytical") tmp.ci=CI.analytical.est.s2(tmp.dat,x)[2:3]
      if(type=="parametric") tmp.ci=parametric.boot.cis(tmp.dat,x,nboot1)[2:3]
      if(type=="hessian") tmp.ci=hessian.boot.cis(tmp.dat,x,nboot1)[2:3]
    }
    if(type == "residuals"){
      y=alpha + beta*dat$x + sample(fit$residuals, nsamp)
      tmp.ci=residuals.boot.cis(tmp.dat,x,nboot1)[2:3]
    }
    if(type == "resampling"){
      tmp.dat = dat[sample(nsamp, replace=TRUE),]
      tmp.ci=resampling.boot.cis(tmp.dat,x,nboot1)[2:3]
    }
    ci.width.boot[i]=tmp.ci[2]-tmp.ci[1]
  }
  #Trim to deal with random ouliers
  mean(ci.width/ci.width.boot, trim=.1)
}


