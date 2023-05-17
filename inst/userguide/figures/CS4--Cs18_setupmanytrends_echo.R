###################################################
### code chunk number 25: Cs18_setupmanytrends_echo (eval = FALSE)
###################################################
## # set new control params
## cntl.list <- list(minit = 200, maxit = 5000, allow.degen = FALSE)
## # set up forms of R matrices
## levels.R <- c(
##   "diagonal and equal",
##   "diagonal and unequal",
##   "equalvarcov",
##   "unconstrained"
## )
## model.data <- data.frame(stringsAsFactors = FALSE)
## # fit lots of models & store results
## # NOTE: this will take a long time to run!
## for (R in levels.R) {
##   for (m in 1:(N.ts - 1)) {
##     dfa.model <- list(A = "zero", R = R, m = m)
##     kemz <- MARSS(dat.z,
##       model = dfa.model, control = cntl.list,
##       form = "dfa", z.score = TRUE
##     )
##     model.data <- rbind(
##       model.data,
##       data.frame(
##         R = R,
##         m = m,
##         logLik = kemz$logLik,
##         K = kemz$num.params,
##         AICc = kemz$AICc,
##         stringsAsFactors = FALSE
##       )
##     )
##     assign(paste("kemz", m, R, sep = "."), kemz)
##   } # end m loop
## } # end R loop


