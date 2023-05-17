###################################################
### code chunk number 31: Cs21_hood.uqr.models
###################################################
Q3 <- matrix(list("offdiag"), 5, 5)
diag(Q3) <- "q"
Q3[, 5] <- 0
Q3[5, ] <- 0
Q3[5, 5] <- "q.hc"
Q.models <- list("equalvarcov", "unconstrained", Q3)
names(Q.models) <- c("equalvarcov", "unconstrained", "hood.independent")


