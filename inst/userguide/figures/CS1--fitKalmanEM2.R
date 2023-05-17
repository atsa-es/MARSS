###################################################
### code chunk number 9: fitKalmanEM2
###################################################
kem.params <- matrix(NA, nrow = 11, ncol = 3, dimnames = list(c(paste("sim", 1:9), "mean sim", "true"), c("kem.U", "kem.Q", "kem.R")))
kem.states <- matrix(NA, nrow = 9, ncol = nYr)
for (i in 1:9) {
  kem <- MARSS(y.tss[i, ], silent = TRUE)
  kem.params[i, ] <- coef(kem, type = "vector")[c(2, 3, 1)]
  kem.states[i, ] <- kem$states
}
kem.params[10, ] <- apply(kem.params[1:9, ], 2, mean)
kem.params[11, ] <- c(sim.u, sim.Q, sim.R)


