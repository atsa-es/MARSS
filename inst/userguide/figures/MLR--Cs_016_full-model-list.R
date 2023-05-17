###################################################
### code chunk number 19: Cs_016_full-model-list
###################################################
eVar.names <- colnames(longley)[-7]
eVar <- t(longley[, eVar.names])
longley.model <- list()
longley.model$U <- longley.model$Q <- "zero"
longley.model$C <- "zero"
longley.model$B <- longley.model$Z <- "identity"
longley.model$A <- matrix("intercept")
longley.model$R <- matrix("r")
longley.model$D <- matrix(eVar.names, nrow = 1)
longley.model$d <- eVar
longley.model$x0 <- "zero"
longley.model$tinitx <- 0


