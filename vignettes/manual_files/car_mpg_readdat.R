# Data compiled by 
# Rhett Allain
# Associate Professor of Physics
# Department of Chemistry and Physics
# email: rallain@selu.edu
# http://www.wired.com/2012/08/fuel-economy-vs-mass/

library(stringr)
dat=read.csv("vignettes/manual_files/mpg.csv")
dat=dat[!str_detect(dat$model, "Aston Marton"),]
dat=dat[!str_detect(dat$model, "Lamborghini"),]
dat=dat[dat$wt!=0,]
dat$wt=dat$wt/1000
fit=lm(city~wt,data=dat)
n=dim(dat)[1]
alpha=coef(fit)[1]
beta=coef(fit)[2]
sigma=sqrt(sum(fit$residual^2)/fit$df.residual)
