###################################################
### code chunk number 23: Covar_sec6_13_show-aics
###################################################
data.frame(
  Model = c("Fixed", "Cubic", "Fourier"),
  AICc = round(c(
    seas.mod.1$AICc,
    seas.mod.2$AICc,
    seas.mod.3$AICc
  ), 1),
  stringsAsFactors = FALSE
)


