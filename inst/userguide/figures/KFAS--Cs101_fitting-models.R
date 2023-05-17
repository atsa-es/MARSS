###################################################
### code chunk number 3: Cs101_fitting-models
###################################################
model_Nile <- SSModel(Nile ~ SSMtrend(
  degree = 1,
  Q = list(matrix(NA))
),
H = matrix(NA)
)
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_default <- fitSSM(model_Nile, kinits, method = "BFGS")


