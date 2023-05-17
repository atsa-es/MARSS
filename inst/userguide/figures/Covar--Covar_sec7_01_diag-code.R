###################################################
### code chunk number 24: Covar_sec7_01_diag-code (eval = FALSE)
###################################################
## for (i in 1:3) {
##   dev.new()
##   modn <- paste("seas.mod", i, sep = ".")
##   for (j in 1:5) {
##     plot.ts(MARSSresiduals(modn, type = "tt1")$model.residuals[j, ],
##       ylab = "Residual", main = phytos[j]
##     )
##     abline(h = 0, lty = "dashed")
##     acf(MARSSresiduals(modn, type = "tt1")$model.residuals[j, ],
##       na.action = na.pass
##     )
##   }
## }


