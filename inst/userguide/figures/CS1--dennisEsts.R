###################################################
### code chunk number 13: dennisEsts
###################################################
den91 <- lm(delta.pop ~ -1 + tau)
# note: the "-1" specifies no intercept
den91.u <- den91$coefficients
den91.Q <- var(resid(den91))
# type summary(den91) to see other info about our regression fit


