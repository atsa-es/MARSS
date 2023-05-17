###################################################
### code chunk number 49: Cs501_plotting
###################################################
ts.plot(cbind(Nile, pred_kfas[, c("fit", "lwr", "upr")], conf_kfas[, c("lwr", "upr")]),
  col = c(1:2, 3, 3, 4, 4),
  ylab = "Predicted annual flow", main = "River Nile"
)


