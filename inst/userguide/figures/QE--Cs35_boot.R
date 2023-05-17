###################################################
### code chunk number 38: Cs35_boot
###################################################
boot.params <- MARSSboot(kem,
  nboot = 20, output = "parameters", sim = "parametric"
)$boot.params


