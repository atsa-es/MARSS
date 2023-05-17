###################################################
### code chunk number 33: Cs30_marsskf
###################################################
kf <- MARSSkf(kem)
names(kf)
# if you only need the logLik, 
MARSSkf(kem, only.logLik = TRUE)
# or
logLik(kem)


