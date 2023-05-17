###################################################
### code chunk number 29: Cs_403_fig-arss-model
###################################################
Z <- matrix(c(1, 0), 1, 2)
B <- matrix(list("b1", 1, "b2", 0), 2, 2)
U <- matrix(0, 2, 1)
Q <- matrix(list("q", 0, 0, 0), 2, 2)
A <- matrix(0, 1, 1)
R <- matrix("r")
V <- matrix(0, 2, 2)
pi <- matrix(mean(noisy.data), 2, 1)
model.list.2ss <- list(
  Z = Z, B = B, U = U, Q = Q, A = A,
  R = R, x0 = pi, V0 = V, tinitx = 0
)


