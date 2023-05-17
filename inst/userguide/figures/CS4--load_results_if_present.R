###################################################
### code chunk number 18: load_results_if_present
###################################################
# This is being done to speed up building the user guide
file <- "CS4--model_fits.RData"
if (file %in% dir("./manual_files")) {
  load(paste("./manual_files/", file, sep = ""))
  saved.res <- TRUE
} else {
  saved.res <- FALSE
}


