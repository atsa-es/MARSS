###################################################
### code chunk number 24: Cs307_obs-filtering
###################################################
ytT_fit <- fitted(fit_marss, type = "ytT")$.fitted
ytT_hatyt <- MARSShatyt(fit_marss)$ytT
cbind(
  a = kf_kfas$muhat[1:n], fitted = ytT_fit[1:n],
  MARSShatyt = ytT_hatyt[1:n], Nile = Nile[1:n]
)


