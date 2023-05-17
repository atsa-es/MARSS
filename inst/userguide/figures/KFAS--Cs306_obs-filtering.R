###################################################
### code chunk number 23: Cs306_obs-filtering
###################################################
var.Eytt1_fit <-
  fitted(fit_marss, type = "ytt1", interval = "confidence")$.se^2
var.Eytt1_hatyt <-
  MARSShatyt(fit_marss, only.kem = FALSE)$var.Eytt1
cbind(
  P_mu = kf_kfas$P_mu[1:n], fitted = var.Eytt1_fit[1:n],
  MARSShatyt = var.Eytt1_hatyt[1:n]
)


