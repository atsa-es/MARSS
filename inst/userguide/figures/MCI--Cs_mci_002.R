###################################################
### code chunk number 3: Cs_mci_002
###################################################
covariates <- rbind(
  Temp = fulldat[years, "Temp"],
  TP = fulldat[years, "TP"]
)
# demean the covariates
the.mean <- apply(covariates, 1, mean, na.rm = TRUE)
covariates <- covariates - the.mean


