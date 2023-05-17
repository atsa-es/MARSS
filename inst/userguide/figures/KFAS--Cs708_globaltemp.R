###################################################
### code chunk number 66: Cs708_globaltemp (eval = FALSE)
###################################################
## # test; these should be identical
## kfas_test <- kfas_temp_stoch
## mod.list <- list(
##   Z = matrix(1, 2, 1),
##   R = kfas_test$model$H[, , 1],
##   U = matrix(0),
##   A = matrix(0, 2, 1),
##   Q = matrix(kfas_test$model$Q[, , 1]),
##   tinitx = 1
## )
## mod.list$x0 <- matrix(0)
## mod.list$V0 <- kfas_test$model$P1
## marss_test <- MARSS(t(GlobalTemp), model = mod.list)
## kfas_test_mod <- MARSSkfas(marss_test,
##   return.kfas.model = TRUE,
##   return.lag.one = FALSE
## )$kfas.model
## 
## kfs <- KFS(kfas_test_mod, smoothing = "disturbance")
## test <- cbind(b = kfas_test_mod$Q[1, 1, 1] - kfs$V_eta[1, 1, ], a = MARSSresiduals(marss_test, type = "tT")$var.residuals[3, 3, ])
## test <- as.data.frame(test)
## test$diff <- test$b - test$a
## head(test)
## tail(test)
## 
## MARSSkfas(marss_test)$VtT[, , 1]
## MARSSkfss(marss_test)$VtT[, , 1]
## 
## 
## var.EytT_fit <-
##   fitted(marss_test, type = "ytT", interval = "confidence")$.se^2
## cbind(V_mu = KFS(kfas_test$model)$V_mu[1, 1, ], fitted = var.EytT_fit)
## 
## kfs <- KFS(kfas_test_mod, smoothing = "disturbance")
## resid_kfas <- rstandard(kfs,
##   type = "state",
##   standardization_type = "marginal"
## )
## resid_marss <- residuals(marss_test,
##   type = "tT",
##   standardization = "marginal"
## )
## resid_marss <- subset(resid_marss, name == "state")
## df <- data.frame(
##   MARSS = resid_marss$.std.resids,
##   KFAS = as.vector(resid_kfas),
##   diff = resid_marss$.std.resids - as.vector(resid_kfas)
## )
## head(df)
## 
## resid_kfas <- residuals(kfs, type = "state")
## resid_marss <- residuals(marss_test, type = "tT")
## resid_marss <- subset(resid_marss, name == "state")
## df <- cbind(
##   MARSS = resid_marss$.resids,
##   KFAS = as.vector(resid_kfas),
##   diff = resid_marss$.resids - as.vector(resid_kfas)
## )
## head(df)
## 
## 
## kfs <- KFS(kfas_temp_stoch$model)
## resid_kfas <- residuals(kfs, type = "pearson")
## resid_marss <- MARSSresiduals(marss_test, type = "tT")
## df <- cbind(
##   as.data.frame(t(resid_marss$model.residuals)),
##   as.data.frame(resid_kfas),
##   as.data.frame(t(resid_marss$model.residuals) - resid_kfas)
## )


