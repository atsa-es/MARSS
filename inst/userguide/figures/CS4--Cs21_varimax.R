###################################################
### code chunk number 29: Cs21_varimax
###################################################
# get the inverse of the rotation matrix
Z.est <- coef(the.fit, type = "matrix")$Z
H.inv <- 1
if (ncol(Z.est) > 1) 
  H.inv <- varimax(coef(the.fit, type = "matrix")$Z)$rotmat


