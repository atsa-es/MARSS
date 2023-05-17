###################################################
### code chunk number 11: Cs202_state-filtering
###################################################
kf_kfas <- KFS(fit_kfas$model,
  filtering = "state",
  smoothing = "state", simplify = FALSE
)


