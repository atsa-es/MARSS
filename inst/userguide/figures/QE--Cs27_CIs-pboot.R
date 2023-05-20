###################################################
### code chunk number 30: Cs27_CIs-pboot
###################################################
kem.w.boot.CIs <- MARSSparamCIs(kem, method = "parametric", nboot = 10)
print(kem.w.boot.CIs)


