###################################################
### code chunk number 11: Cs03_set.up.Z.models
###################################################
# H1 stock
Z1 <- factor(c("wa.or", "wa.or", rep("ps", 4), "ca", "ca", "wa.or", "wa.or", "bc"))
# H2 coastal+PS
Z2 <- factor(c(rep("coast", 2), rep("ps", 4), rep("coast", 4), "ps"))
# H3 N and S
Z3 <- factor(c(rep("N", 6), "S", "S", "N", "S", "N"))
# H4 North Coast, Inland Strait, Puget Sound, South Coast
Z4 <- factor(c("nc", "nc", "is", "is", "ps", "ps", "sc", "sc", "nc", "sc", "is"))
# H5 panmictic
Z5 <- factor(rep("pan", 11))
# H6 Site
Z6 <- factor(1:11) # site
Z.models <- list(Z1, Z2, Z3, Z4, Z5, Z6)
names(Z.models) <-
  c("stock", "coast+PS", "N-S", "NC+Strait+PS+SC", "panmictic", "site")


