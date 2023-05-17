###################################################
### code chunk number 42: Cs39_AICb
###################################################
kem.with.AICb <- MARSSaic(kem,
  output = "AICbp",
  Options = list(nboot = 10, silent = TRUE)
)

print(kem.with.AICb)


