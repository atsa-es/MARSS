###################################################
### code chunk number 16: Cs07_add.delta.aicc
###################################################
out.tab.1 <- cbind(out.tab.1,
  delta.AICc = out.tab.1$AICc - out.tab.1$AICc[1]
)


