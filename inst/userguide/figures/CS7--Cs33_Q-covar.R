###################################################
### code chunk number 36: Cs33_Q-covar
###################################################
Q <- matrix(list(0), 5, 5)
Q[1:4, 1:4] <- paste(rep(1:4, times = 4), rep(1:4, each = 4), sep = "")
Q[5, 5] <- "fish"
Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
print(Q)


