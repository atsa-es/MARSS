###################################################
### code chunk number 12: Cs_010_example1-lm
###################################################
mod1.lm <- lm(Employed ~ GNP + Population, data = longley)
coef(mod1.lm)


