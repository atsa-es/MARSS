###################################################
### code chunk number 10: Cs201_state-filtering
###################################################
fit_kfas <- fit_kfas_stoch
fit_marss <- fit_em_stoch
fit_marss$par$Q[1, 1] <- fit_kfas$model$Q
fit_marss$par$R[1, 1] <- fit_kfas$model$H


