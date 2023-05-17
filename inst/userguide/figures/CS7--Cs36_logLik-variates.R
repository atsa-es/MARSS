###################################################
### code chunk number 40: Cs36_logLik-variates
###################################################
tmp <- kem.plank.5
tmp$marss$data[5, ] <- NA
LL.variates <- MARSSkf(tmp)$logLik


