###################################################
### code chunk number 22: testing
###################################################
U.model <- array(list(), dim = c(1, 1, length(rockfish[, "Year"])))
U.model[1, 1, rockfish[, "Year"] < 1980] <- "pre-1980"
U.model[1, 1, rockfish[, "Year"] >= 1980] <- "post-1980"
model4 <- list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem4 <- MARSS(fishdat, model = model4)

Z.model <- factor(c("trawl", "trawl", "rec", "rec", "rec", "rec", "rec", "scuba", "wdfw.trawl"))
Q.model <- "diagonal and unequal"
U.model <- "equal"
model5 <- list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem5 <- MARSS(fishdat, model = model5)


