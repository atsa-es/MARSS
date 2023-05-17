###################################################
### code chunk number 44: Cs507_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- residuals(kfs, type = "response")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


