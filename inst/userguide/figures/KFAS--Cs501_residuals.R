###################################################
### code chunk number 38: Cs501_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- residuals(kfs, type = "recursive")
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


