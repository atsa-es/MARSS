###################################################
### code chunk number 22: Cs305_obs-filtering
###################################################
ytt1_fit <- fitted(fit_marss, type = "ytt1")$.fitted
ytt1_hatyt <- MARSShatyt(fit_marss, only.kem = FALSE)$ytt1
cbind(m = kf_kfas$m[1:n], fitted = ytt1_fit[1:n], MARSShatyt = ytt1_hatyt[1:n])


