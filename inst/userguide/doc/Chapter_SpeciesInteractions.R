###################################################
### code chunk number 2: Cs01_load.wolf.data
###################################################
yr1960to2011 <- isleRoyal[, "Year"] >= 1960 & isleRoyal[, "Year"] <= 2011
royale.dat <- log(t(isleRoyal[yr1960to2011, c("Wolf", "Moose")]))


###################################################
### code chunk number 3: Cs02_plotwolfmoosedata
###################################################
x <- isleRoyal[, "Year"]
y <- log(isleRoyal[, c("Wolf", "Moose")])
graphics::matplot(x, y,
  ylab = "Log count", xlab = "Year", type = "l",
  lwd = 3, bty = "L", col = "black"
)
legend("topright", c("Wolf", "Moose"), lty = c(1, 2), bty = "n")


###################################################
### code chunk number 5: Cs03_z.score.wolf.data
###################################################
# if missing values are in the data, they should be NAs
z.royale.dat <- zscore(royale.dat)


###################################################
### code chunk number 6: Cs04_fit.model
###################################################
royale.model.1 <- list(
  Z = "identity", B = "unconstrained",
  Q = "diagonal and unequal", R = "diagonal and unequal",
  U = "zero", tinitx = 1
)
cntl.list <- list(allow.degen = FALSE, maxit = 200)
kem.1 <- MARSS(z.royale.dat, model = royale.model.1, control = cntl.list)


###################################################
### code chunk number 7: Cs05_fit.model.R0
###################################################
royale.model.2 <- list(
  Z = "identity", B = "unconstrained",
  Q = "diagonal and unequal", R = "zero", U = "zero"
)
kem.2 <- MARSS(z.royale.dat, model = royale.model.2)


###################################################
### code chunk number 14: Cs05_fit.model.tinitx1
###################################################
royale.model.4 <- list(
  B = "unconstrained", U = "zero", Q = "diagonal and unequal",
  Z = "identity", R = "zero", tinitx = 1
)
kem.4 <- MARSS(z.royale.dat, model = royale.model.4)


###################################################
### code chunk number 8: Cs06_print-wolf.B
###################################################
wolf.B <- coef(kem.2, type = "matrix")$B
rownames(wolf.B) <- colnames(wolf.B) <- rownames(royale.dat)
print(wolf.B, digits = 2)


###################################################
### code chunk number 9: Cs07_prep-cov-wolf-moose
###################################################
clim.variables <- c(
  "jan.feb.ave.temp", "jan.feb.ave.precip",
  "july.sept.ave.temp"
)
yr1959to2010 <- isleRoyal[, "Year"] >= 1959 & isleRoyal[, "Year"] <= 2010
clim.dat <- t(isleRoyal[yr1959to2010, clim.variables])
z.score.clim.dat <- zscore(clim.dat)


###################################################
### code chunk number 10: Cs08_cov.wolf.moose.model
###################################################
royale.model.3 <- list(
  Z = "identity", B = "unconstrained",
  Q = "diagonal and unequal", R = "zero", U = "zero",
  C = matrix(list(
    0, "Moose win temp", 0, "Moose win precip",
    0, "Moose sum temp"
  ), 2, 3),
  c = z.score.clim.dat
)


###################################################
### code chunk number 11: Cs09_fit-cov-wolf-moose-model
###################################################
kem.3 <- MARSS(z.royale.dat, model = royale.model.3)


###################################################
### code chunk number 12: Cs10_figwolfcov
###################################################
cor.fun <- function(x, y) {
  text(0.5, 0.5, format(cor(x, y), digits = 2), cex = 2)
}
pairs(t(z.score.clim.dat), lower.panel = cor.fun)


###################################################
### code chunk number 13: Cs11_bad-data-2 (eval = FALSE)
###################################################
## bad.data <- z.royale.dat + matrix(rnorm(100, 0, sqrt(.2)), 2, 50)
## kem.bad <- MARSS(bad.data, model = model)


###################################################
### code chunk number 15: Cs12_load-plankton-data
###################################################
# only use the plankton, daphnia, & non-daphnia
plank.spp <- c("Large Phyto", "Small Phyto", "Daphnia", "Non-daphnia")
plank.dat <- ivesDataByWeek[, plank.spp]
# The data are not logged
plank.dat <- log(plank.dat)
# Transpose to get time going across the columns
plank.dat <- t(plank.dat)
# make a demeaned version
d.plank.dat <- (plank.dat - apply(plank.dat, 1, mean, na.rm = TRUE))


###################################################
### code chunk number 16: Cs13_plot-plankton-data
###################################################
graphics::matplot((1:(52 * 6))[27:295], t(d.plank.dat), type = "l", lty = c(1, 1, 1, 1), lwd = c(1, 1, 3, 3), xlab = "week of study", ylab = "log biomass", xaxt = "n", xlim = c(11, 52 * 6 - 11), bty = "L")
# axis(1,at=(1:(52*6))[seq(27,295,2)])
axis(1, at = seq(1, 52 * 6, 2))
abline(v = c(52 * (1:6)))
abline(h = 0)


###################################################
### code chunk number 17: Cs14_set-up-plankton-model
###################################################
Q <- matrix(list(0), 4, 4)
diag(Q) <- c("Phyto", "Phyto", "Zoo", "Zoo")
R <- matrix(list(0), 4, 4)
diag(R) <- c("Phyto", "Phyto", "Zoo", "Zoo")
plank.model.0 <- list(
  B = "unconstrained", U = "zero", Q = Q,
  Z = "identity", A = "zero", R = R,
  x0 = "unequal", tinitx = 1
)


###################################################
### code chunk number 18: Cs15_fit-plank-model-0
###################################################
kem.plank.0 <- MARSS(d.plank.dat, model = plank.model.0)


###################################################
### code chunk number 19: Cs16_print-B-0
###################################################
# Cleaning up the B matrix for printing
B.0 <- coef(kem.plank.0, type = "matrix")$B[1:4, 1:4]
rownames(B.0) <- colnames(B.0) <- c("LP", "SP", "D", "ND")
print(B.0, digits = 2)


###################################################
### code chunk number 20: Cs17_print-B-Ives
###################################################
# Cleaning up the B matrix for printing
B.Ives.ML <- matrix(c(.5, NA, NA, NA, -.39, .076, NA, .1, NA, -.02, .77, NA, NA, -.1, NA, .55), 4, 4)
B.Ives.Obs <- matrix(c(.48, NA, NA, NA, -.39, .25, NA, .1, NA, -.17, .74, 0, NA, -.11, 0, .6), 4, 4)
B.Ives <- B.Ives.Obs
rownames(B.Ives) <- colnames(B.Ives) <- c("LP", "SP", "D", "ND")
print(B.Ives, digits = 2, na.print = "--")


###################################################
### code chunk number 21: Cs18_test-rm-NAs (eval = FALSE)
###################################################
## # Example code to see what would happen if we removed the NAs
## test.dat <- d.plank.dat[, !is.na(d.plank.dat[1, ])]
## test <- MARSS(test.dat, model = plank.model.0)


###################################################
### code chunk number 22: Cs19_fit-plank-model-1
###################################################
plank.model.1 <- plank.model.0
plank.model.1$Q <- "unconstrained"
kem.plank.1 <- MARSS(d.plank.dat, model = plank.model.1)


###################################################
### code chunk number 23: Cs20_print-B-1
###################################################
# Cleaning up the B matrix for printing
B <- coef(kem.plank.1, type = "matrix")$B[1:4, 1:4]
rownames(B) <- colnames(B) <- c("LP", "SP", "D", "ND")
B[B == 0] <- NA
B.1 <- B
print(B, digits = 2, na.print = "--")


###################################################
### code chunk number 24: Cs21_B-2
###################################################
B.2 <- matrix(list(0), 4, 4) # set up the list matrix
diag(B.2) <- c("B11", "B22", "B33", "B44") # give names to diagonals
# and names to the estimated non-diagonals
B.2[1, 2] <- "B12"
B.2[2, 3] <- "B23"
B.2[2, 4] <- "B24"
B.2[4, 2] <- "B42"
print(B.2)


###################################################
### code chunk number 25: Cs22_fit-plank-model-2
###################################################
# model 2
plank.model.2 <- plank.model.1
plank.model.2$B <- B.2
kem.plank.2 <- MARSS(d.plank.dat, model = plank.model.2)


###################################################
### code chunk number 26: Cs23_print-B-2
###################################################
# Cleaning up the B matrix for printing
B <- coef(kem.plank.2, type = "matrix")$B[1:4, 1:4]
rownames(B) <- colnames(B) <- c("LP", "SP", "D", "ND")
B[B == 0] <- NA
B.2 <- B
print(B, digits = 2, na.print = "--")


###################################################
### code chunk number 27: Cs24_fit-plank-model-3
###################################################
# model 3
plank.model.3 <- plank.model.2
plank.model.3$R <- diag(c(.04, .04, .16, .16))
kem.plank.3 <- MARSS(d.plank.dat, model = plank.model.3)


###################################################
### code chunk number 28: Cs25_prep-covariates
###################################################
# transpose to make time go across columns
# drop=FALSE so that R doesn't change our matrix to a vector
phos <- t(log(ivesDataByWeek[, "Phosph", drop = FALSE]))
d.phos <- (phos - apply(phos, 1, mean, na.rm = TRUE))


###################################################
### code chunk number 29: Cs26_add-covar-model-3
###################################################
plank.model.4 <- plank.model.3
plank.model.4$C <- matrix(list("C11", "C21", 0, 0), 4, 1)
plank.model.4$c <- d.phos


###################################################
### code chunk number 30: Cs27_plank-model-4
###################################################
kem.plank.4 <- MARSS(d.plank.dat, model = plank.model.4)


###################################################
### code chunk number 31: Cs28_print-C
###################################################
# Cleaning up the C matrix for printing
Cmat <- coef(kem.plank.4, type = "matrix")$C[1:4, 1, drop = FALSE]
rownames(Cmat) <- c("LP", "SP", "D", "ND")
Cmat[Cmat == 0] <- NA
print(Cmat, digits = 2, na.print = "--")


###################################################
### code chunk number 32: Cs29_add-fish-to-data
###################################################
# transpose to make time go across columns
# drop=FALSE so that R doesn't change our matrix to a vector
fish <- t(log(ivesDataByWeek[, "Fish biomass", drop = FALSE]))
d.fish <- (fish - apply(fish, 1, mean, na.rm = TRUE))
# plank.dat.w.fish = rbind(plank.dat,fish)
d.plank.dat.w.fish <- rbind(d.plank.dat, d.fish)


###################################################
### code chunk number 33: Cs30_B-covar
###################################################
B <- matrix(list(0), 5, 5)
diag(B) <- list("B11", "B22", "B33", "B44", "Bfish")
B[1, 2] <- "B12"
B[2, 3] <- "B23"
B[2, 4] <- "B24"
B[4, 2] <- "B42"
B[1:4, 5] <- list(0, 0, "C32", "C42")
print(B)


###################################################
### code chunk number 34: Cs31_C-covar
###################################################
C <- matrix(list("C11", "C21", 0, 0, 0), 5, 1)


###################################################
### code chunk number 35: Cs32_R.covar
###################################################
R <- matrix(list(0), 5, 5)
diag(R) <- list(0.04, 0.04, 0.16, 0.16, 0.36)


###################################################
### code chunk number 36: Cs33_Q-covar
###################################################
Q <- matrix(list(0), 5, 5)
Q[1:4, 1:4] <- paste(rep(1:4, times = 4), rep(1:4, each = 4), sep = "")
Q[5, 5] <- "fish"
Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
print(Q)


###################################################
### code chunk number 37: Cs34_fit-covar-model
###################################################
plank.model.5 <- plank.model.4
plank.model.5$B <- B
plank.model.5$C <- C
plank.model.5$Q <- Q
plank.model.5$R <- R
kem.plank.5 <- MARSS(d.plank.dat.w.fish, model = plank.model.5)


###################################################
### code chunk number 38: Cs35_print-B
###################################################
# Cleaning up the B matrix for printing
B.5 <- coef(kem.plank.5, type = "matrix")$B[1:4, 1:4]
rownames(B.5) <- colnames(B.5) <- c("LP", "SP", "D", "ND")
B.5[B.5 == 0] <- NA
print(B.5, digits = 2, na.print = "--")


###################################################
### code chunk number 40: Cs36_logLik-variates
###################################################
tmp <- kem.plank.5
tmp$marss$data[5, ] <- NA
LL.variates <- MARSSkf(tmp)$logLik


###################################################
### code chunk number 41: Cs37_BQ.5
###################################################
B <- coef(kem.plank.5, type = "matrix")$B[1:4, 1:4]
Q <- coef(kem.plank.5, type = "matrix")$Q[1:4, 1:4]


###################################################
### code chunk number 42: Cs38_max.eigen
###################################################
max(eigen(B)$values)


###################################################
### code chunk number 43: Cs39_max.eig.kron.b
###################################################
max(eigen(kronecker(B, B))$values)


###################################################
### code chunk number 44: Cs40_Vinfty
###################################################
m <- nrow(B)
vecV <- solve(diag(m * m) - kronecker(B, B)) %*% as.vector(Q)
V_inf <- matrix(vecV, nrow = m, ncol = m)


###################################################
### code chunk number 45: Cs41_det.b.squared
###################################################
abs(det(B))^2


###################################################
### code chunk number 46: Cs42_det.b.scaled
###################################################
abs(det(B))^(2 / nrow(B))


###################################################
### code chunk number 47: Cs43_covar.sigma.Vinf
###################################################
-sum(diag(Q)) / sum(diag(V_inf))


###################################################
### code chunk number 48: Cs44_worse.case.reactivity
###################################################
max(eigen(t(B) %*% B)$values) - 1


