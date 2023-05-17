###################################################
### code chunk number 22: set_up_two_trends_big_maxit
###################################################
if (!saved.res) {
  model.list <- list(m = 2, R = "diagonal and unequal")
  kemz.2 <- MARSS(dat.spp.1980,
    model = model.list,
    z.score = TRUE, form = "dfa", control = big.maxit.cntl.list
  )
}


