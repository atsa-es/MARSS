###################################################
### code chunk number 16: Covar_sec6_06_fit-month-factor-with-MARSS
###################################################
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A, R = R,
  C = C, c = c.in, D = D, d = d
)
seas.mod.1 <- MARSS(dat, model = model.list, control = list(maxit = 1500))

# Get the estimated seasonal effects
# rows are taxa, cols are seasonal effects
seas.1 <- coef(seas.mod.1, type = "matrix")$C
rownames(seas.1) <- phytos
colnames(seas.1) <- month.abb


