###################################################
### code chunk number 41: save_results
###################################################
file <- paste(Sys.getenv("TEMP"), "\\CS4--model_fits.RData", sep = "")
save(file = file, list = c("model.data", ls(pattern = "^kemz.")))


