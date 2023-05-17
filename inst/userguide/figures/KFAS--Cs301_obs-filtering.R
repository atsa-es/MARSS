###################################################
### code chunk number 19: Cs301_obs-filtering
###################################################
kf_kfas <- KFS(fit_kfas$model,
  filtering = "signal",
  smoothing = "signal", simplify = FALSE
)


