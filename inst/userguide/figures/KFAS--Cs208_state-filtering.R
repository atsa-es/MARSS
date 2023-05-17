###################################################
### code chunk number 17: Cs208_state-filtering
###################################################
cbind(
  v = kf_kfas$v[1:n], Innov = kf_marss$Innov[1:n],
  F = kf_kfas$F[1:n], Sigma = kf_marss$Sigma[1:n]
)


