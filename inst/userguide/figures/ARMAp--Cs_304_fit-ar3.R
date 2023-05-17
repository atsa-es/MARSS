###################################################
### code chunk number 26: Cs_304_fit-ar3
###################################################
print(cbind(
  true = true.3[c("b1", "b2", "b3", "q")],
  estimates.no.miss = coef(ar3, type = "vector")
))


