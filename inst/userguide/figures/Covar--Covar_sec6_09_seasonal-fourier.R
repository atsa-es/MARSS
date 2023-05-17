###################################################
### code chunk number 19: Covar_sec6_09_seasonal-fourier
###################################################
cos.t <- cos(2 * pi * seq(TT) / period)
sin.t <- sin(2 * pi * seq(TT) / period)
c.Four <- rbind(cos.t, sin.t)


