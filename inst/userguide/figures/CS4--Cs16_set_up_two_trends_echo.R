###################################################
### code chunk number 21: Cs16_set_up_two_trends_echo
###################################################
model.list <- list(m = 2, R = "diagonal and unequal")
kemz.2 <- MARSS(dat.spp.1980,
  model = model.list,
  z.score = TRUE, form = "dfa", control = cntl.list
)


