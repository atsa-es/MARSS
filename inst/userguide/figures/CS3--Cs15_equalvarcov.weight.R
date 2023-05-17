###################################################
### code chunk number 25: Cs15_equalvarcov.weight
###################################################
c(
  sum(out.tab.2$AIC.weight[out.tab.2$Q == "equalvarcov"]),
  sum(out.tab.2$AIC.weight[out.tab.2$Q == "unconstrained"]),
  sum(out.tab.2$AIC.weight[out.tab.2$Q == "diagonal and equal"])
)


