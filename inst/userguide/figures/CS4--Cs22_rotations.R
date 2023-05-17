###################################################
### code chunk number 30: Cs22_rotations
###################################################
# rotate factor loadings
Z.rot <- Z.est %*% H.inv
# rotate trends
trends.rot <- solve(H.inv) %*% the.fit$states


