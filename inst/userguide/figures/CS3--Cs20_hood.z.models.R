###################################################
### code chunk number 30: Cs20_hood.z.models
###################################################
ZH1 <- factor(c("nc", "nc", "is", "is", "ps",
  "ps", "sc", "sc", "nc", "sc", "is", "ps"))
ZH2 <- factor(c("nc", "nc", "is", "is", "ps",
  "ps", "sc", "sc", "nc", "sc", "is", "hc"))
Z.models.hc <- list(ZH1, ZH2)
names(Z.models.hc) <- c("hood.in.ps", "hood.separate")


