###################################################
### code chunk number 48: Cs510_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
test <- cbind(
  b = fit_kfas$model$Q[1, 1, 1] - kfs$V_eta[1, 1, ],
  a = MARSSresiduals(fit_marss, type = "tT")$var.residuals[2, 2, ]
)
test <- as.data.frame(test)
test$diff <- test$b - test$a
head(test)
tail(test)


