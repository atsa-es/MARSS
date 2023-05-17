###################################################
### code chunk number 18: Cs09_aic.weight
###################################################
out.tab.1 <- cbind(out.tab.1,
  AIC.weight = out.tab.1$rel.like / sum(out.tab.1$rel.like)
)


