###################################################
### code chunk number 31: Cs22_rotations_cis
###################################################
# Add CIs to marssMLE object
the.fit <- MARSSparamCIs(the.fit)
# Use coef() to get the upper and lower CIs
Z.low <- coef(the.fit, type = "Z", what = "par.lowCI")
Z.up <- coef(the.fit, type = "Z", what = "par.upCI")
Z.rot.up <- Z.up %*% H.inv
Z.rot.low <- Z.low %*% H.inv
df <- data.frame(
  est = as.vector(Z.rot), 
  conf.up = as.vector(Z.rot.up), 
  conf.low = as.vector(Z.rot.low)
  )


