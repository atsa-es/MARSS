###################################################
### code chunk number 5: Covar_sec2_3_plot-dat
###################################################
LWA <- ts(cbind(Year = fulldat[years, "Year"], t(dat), t(covariates)), start = c(1965, 1), end = c(1974, 12), freq = 12)
plot.ts(LWA[, c("Greens", "Bluegreens", "Temp", "TP")], main = "", yax.flip = TRUE)


