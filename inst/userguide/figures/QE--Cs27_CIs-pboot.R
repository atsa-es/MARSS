###################################################
### code chunk number 30: Cs27_CIs-pboot
###################################################
kem.w.boot.CIs <- MARSSparamCIs(kem, method = "parametric", nboot = 10)
# nboot should be more like 1000, but set low for example's sake
print(kem.w.boot.CIs)


