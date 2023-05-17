###################################################
### code chunk number 46: Cs509_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "state")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 47: Cs509_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_marss <- subset(resid_marss, name == "state")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


