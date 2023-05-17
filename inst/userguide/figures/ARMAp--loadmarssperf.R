###################################################
### code chunk number 35: loadmarssperf
###################################################
nsim <- 200
TT <- 100
file <- paste("AR2SS", TT, ".RData", sep = "")
if (file %in% dir("./manual_files")) {
  load(paste("./manual_files/", file, sep = ""))
  sims.exist <- TRUE
} else {
  sims.exist <- FALSE
}


