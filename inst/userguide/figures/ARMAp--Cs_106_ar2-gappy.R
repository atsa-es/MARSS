###################################################
### code chunk number 9: Cs_106_ar2-gappy
###################################################
print(cbind(
  true = true.2[2:4],
  estimates.no.miss = coef(ar2, type = "vector"),
  estimates.w.miss = coef(ar2.gappy, type = "vector")
))


