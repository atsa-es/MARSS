###################################################
### code chunk number 8: Cs_105_ar2-gappy
###################################################
gappy.data <- sim.ar2[3:TT]
gappy.data[floor(runif(TT / 2, 2, TT))] <- NA
ar2.gappy <- MARSS(gappy.data, model = model.list.2, fun.kf="MARSSkfss")


