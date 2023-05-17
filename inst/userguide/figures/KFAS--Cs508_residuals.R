###################################################
### code chunk number 45: Cs508_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- residuals(kfs, type = "state")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "state")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


