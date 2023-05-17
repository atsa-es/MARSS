###################################################
### code chunk number 33: Cs_407_fit-arss2-model
###################################################
print(cbind(
  true = true.2ss,
  model.no.error = c(NA, coef(ar2ss2, type = "vector")),
  model.w.error = coef(ar2ss, type = "vector")
))


