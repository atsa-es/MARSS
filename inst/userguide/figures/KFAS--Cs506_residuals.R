###################################################
### code chunk number 43: Cs506_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- rstandard(kfs,
  type = "pearson",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "Cholesky"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


