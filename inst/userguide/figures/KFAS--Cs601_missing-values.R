###################################################
### code chunk number 52: Cs601_missing-values
###################################################
NileNA <- Nile
NileNA[c(21:40, 61:80)] <- NA
model_NileNA_stoch <-
  SSModel(NileNA ~ SSMtrend(
    degree = 1,
    Q = list(matrix(NA))
  ),
  H = matrix(NA)
  )
model_NileNA_stoch$a1[1, 1] <- 0
model_NileNA_stoch$P1[1, 1] <- model_Nile_stoch$P1[1, 1]
model_NileNA_stoch$P1inf[1, 1] <- 0
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_NA <- fitSSM(model_NileNA_stoch, kinits, method = "BFGS")
fit_marss_NA <- MARSS(as.vector(NileNA),
  model = mod.nile.stoch,
  inits = inits, method = "BFGS", silent = TRUE
)


