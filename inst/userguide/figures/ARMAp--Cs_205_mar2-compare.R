###################################################
### code chunk number 20: Cs_205_mar2-compare
###################################################
model.list.2$x0 <- matrix(sim.mar2[1, 2:1], 2, 1)
mar2a <- MARSS(sim.mar2[1, 2:TT], model = model.list.2)
model.list.2$x0 <- matrix(sim.mar2[2, 2:1], 2, 1)
mar2b <- MARSS(sim.mar2[2, 2:TT], model = model.list.2)


