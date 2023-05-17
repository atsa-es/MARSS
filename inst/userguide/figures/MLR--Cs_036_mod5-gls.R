###################################################
### code chunk number 39: Cs_036_mod5-gls
###################################################
sleep.mod5.gls <- gls(Reaction ~ 0 + Subject + Days:Subject,
  data = sleepstudy,
  correlation = corAR1(form = ~ 1 | Subject), method = "ML"
)


