###################################################
### code chunk number 4: Cs102_fitting-models
###################################################
model_Nile_stoch <- model_Nile
model_Nile_stoch$a1[1, 1] <- 0
model_Nile_stoch$P1[1, 1] <- 1000
model_Nile_stoch$P1inf[1, 1] <- 0
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_stoch <- fitSSM(model_Nile_stoch, kinits, method = "BFGS")
kfs_kfas_stoch <- KFS(fit_kfas_stoch$model)


