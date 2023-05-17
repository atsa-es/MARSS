###################################################
### code chunk number 29: Cs2_Code5_7
###################################################
# Two subpopulations with different population parameters
Z.model <- factor(c(1, 1, 2, 2, 2))
U.model <- "unequal"
Q.model <- "diagonal and unequal"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z.model, U = U.model, Q = Q.model, R = R.model))

# Hood Canal covaries with the other regions
Z.model <- factor(c(1, 1, 1, 1, 2))
U.model <- "unequal"
Q.model <- "equalvarcov"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z.model, U = U.model, Q = Q.model, R = R.model))

# Three subpopulations with shared parameter values
Z.model <- factor(c(1, 1, 2, 2, 3))
U.model <- "unequal"
Q.model <- "diagonal and unequal"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z.model, U = U.model, Q = Q.model, R = R.model))


