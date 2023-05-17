###################################################
### code chunk number 11: Cs_09_univDLM-fore-Var
###################################################
# variance of regr parameters; 1x2xT array
Phi <- kf.out$Vtt1
# obs variance; 1x1 matrix
R.est <- coef(dlm1, type = "matrix")$R
# ts of Var(forecasts)
fore.var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1) # transpose of Z
  fore.var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R.est
}


