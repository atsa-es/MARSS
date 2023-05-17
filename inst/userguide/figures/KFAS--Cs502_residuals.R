###################################################
### code chunk number 39: Cs502_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- rstandard(kfs,
  type = "recursive",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


