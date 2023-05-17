###################################################
### code chunk number 50: Cs502_plotting
###################################################
plot.type <- ifelse(packageVersion("MARSS") < '3.11.4', "model.ytT", "fitted.ytT")
plot(fit_marss, plot.type = plot.type, pi.int = TRUE)


