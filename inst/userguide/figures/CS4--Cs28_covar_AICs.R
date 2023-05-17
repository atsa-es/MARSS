###################################################
### code chunk number 42: Cs28_covar_AICs
###################################################
print(cbind(
  model = c("no covars", "Temp", "TP", "Temp & TP"),
  AICc = round(c(
    the.fit$AICc, kemz.temp$AICc, kemz.TP$AICc,
    kemz.both$AICc
  ))
), quote = FALSE)


