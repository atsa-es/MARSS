###################################################
### code chunk number 24: Cs_302_set-up-ar3-model
###################################################
Z <- matrix(c(1, 0, 0), 1, 3)
B <- matrix(list("b1", 1, 0, "b2", 0, 1, "b3", 0, 0), 3, 3)
U <- matrix(0, 3, 1)
Q <- matrix(list(0), 3, 3)
Q[1, 1] <- "q"
A <- matrix(0, 1, 1)
R <- matrix(0, 1, 1)
pi <- matrix(sim.ar3[3:1], 3, 1)
V <- matrix(0, 3, 3)
model.list.3 <- list(
  Z = Z, B = B, U = U, Q = Q, A = A,
  R = R, x0 = pi, V0 = V, tinitx = 1
)


