###################################################
### code chunk number 18: Cs209_state-filtering
###################################################
cbind(
  P = kf_kfas$P[1:n], Vtt1 = kf_marss$Vtt1[1:n],
  Ptt = kf_kfas$Ptt[1:n], Vtt = kf_marss$Vtt[1:n]
)


