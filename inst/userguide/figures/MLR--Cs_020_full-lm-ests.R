###################################################
### code chunk number 23: Cs_020_full-lm-ests
###################################################
mod3.lm <- lm(Employed ~ 1 + GNP.deflator + GNP + Unemployed
  + Armed.Forces + Population + Year, data = longley)
c(coef(mod3.lm), logLik = logLik(mod3.lm))


