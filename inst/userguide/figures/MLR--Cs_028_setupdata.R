###################################################
### code chunk number 31: Cs_028_setupdata
###################################################
# number of subjects
nsub <- length(unique(sleepstudy$Subject))
ndays <- length(sleepstudy$Days) / nsub
dat <- matrix(sleepstudy$Reaction, nsub, ndays, byrow = TRUE)
rownames(dat) <- paste("sub", unique(sleepstudy$Subject), sep = ".")
exp.var <- matrix(sleepstudy$Days, 1, ndays, byrow = TRUE)


