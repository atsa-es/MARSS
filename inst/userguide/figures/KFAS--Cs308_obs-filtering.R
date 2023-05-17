###################################################
### code chunk number 25: Cs308_obs-filtering
###################################################
var.EytT_fit <-
  fitted(fit_marss, type = "ytT", interval = "confidence")$.se^2
var.EytT_hatyt <-
  MARSShatyt(fit_marss, only.kem = FALSE)$var.EytT
cbind(
  V_mu = kf_kfas$V_mu[1:n], fitted = var.EytT_fit[1:n],
  MARSShatyt = var.EytT_hatyt[1:n]
)


