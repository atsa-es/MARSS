###################################################
### code chunk number 3: Cs01_read_in_data
###################################################
data(lakeWAplankton)
# we want lakeWAplanktonTrans, which has been log-transformed
# and the 0s replaced with NAs
plankdat <- lakeWAplanktonTrans
years <- plankdat[, "Year"] >= 1980 & plankdat[, "Year"] < 1990
phytos <- c(
  "Cryptomonas", "Diatoms", "Greens",
  "Unicells", "Other.algae"
)
dat.spp.1980 <- plankdat[years, phytos]


