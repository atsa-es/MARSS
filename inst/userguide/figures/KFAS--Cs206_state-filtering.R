###################################################
### code chunk number 15: Cs206_state-filtering
###################################################
cbind(
  a = kf_kfas$a[1:n], xtt1 = kf_marss$xtt1[1:n],
  att = kf_kfas$att[1:n], xtt = kf_marss$xtt[1:n]
)


