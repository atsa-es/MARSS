###################################################
### code chunk number 4: Covar_sec2_2_z-score-covar-data
###################################################
covariates <- rbind(
  Temp = fulldat[years, "Temp"],
  TP = fulldat[years, "TP"]
)
# z.score the covariates
covariates <- zscore(covariates)


