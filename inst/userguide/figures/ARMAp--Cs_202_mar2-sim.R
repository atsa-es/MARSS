###################################################
### code chunk number 15: Cs_202_mar2-sim
###################################################
temp2 <- arima.sim(n = TT, list(ar = true.2[c("b1", "b2")]), sd = sqrt(true.2["q"]))
sim.mar2 <- rbind(temp1, temp2)


