###################################################
### code chunk number 14: Cs05_run.the.models
###################################################
out.tab <- NULL
fits <- list()
for (i in 1:length(Z.models)) {
  for (Q.model in Q.models) {
    fit.model <- c(list(Z = Z.models[[i]], Q = Q.model), model.constant)
    fit <- MARSS(sealData,
      model = fit.model,
      silent = TRUE, control = list(maxit = 1000)
    )
    out <- data.frame(
      H = names(Z.models)[i], Q = Q.model, U = U.model,
      logLik = fit$logLik, AICc = fit$AICc, num.param = fit$num.params,
      m = length(unique(Z.models[[i]])),
      num.iter = fit$numIter, converged = !fit$convergence,
      stringsAsFactors = FALSE
    )
    out.tab <- rbind(out.tab, out)
    fits <- c(fits, list(fit))
    if (i == 5) next # one m for panmictic so only run 1 Q
  }
}


