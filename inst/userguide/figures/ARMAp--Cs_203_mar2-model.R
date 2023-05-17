###################################################
### code chunk number 16: Cs_203_mar2-model
###################################################
Z <- matrix(c(1, 0, 0, 1, 0, 0, 0, 0), 2, 4)
B1 <- matrix(list(0), 2, 2)
diag(B1) <- "b1"
B2 <- matrix(list(0), 2, 2)
diag(B2) <- "b2"
B <- matrix(list(0), 4, 4)
B[1:2, 1:2] <- B1
B[1:2, 3:4] <- B2
B[3:4, 1:2] <- diag(1, 2)
U <- matrix(0, 4, 1)
Q <- matrix(list(0), 4, 4)
Q[1, 1] <- "q"
Q[2, 2] <- "q"
A <- matrix(0, 2, 1)
R <- matrix(0, 2, 2)
pi <- matrix(c(sim.mar2[, 2], sim.mar2[, 1]), 4, 1)
V <- matrix(0, 4, 4)
model.list.2m <- list(
  Z = Z, B = B, U = U, Q = Q, A = A,
  R = R, x0 = pi, V0 = V, tinitx = 1
)


