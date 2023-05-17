###################################################
### code chunk number 41: Cs38_sim-like
###################################################
kem.sim.2 <- kem.sim.1
kem.sim.2$model$data <- sim.data[, , 2]
MARSSkf(kem.sim.2)$logLik


