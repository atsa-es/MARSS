###################################################
### code chunk number 17: Cs08_add.delta.aicc
###################################################
out.tab.1 <- cbind(out.tab.1,
  rel.like = exp(-1 * out.tab.1$delta.AICc / 2)
)


