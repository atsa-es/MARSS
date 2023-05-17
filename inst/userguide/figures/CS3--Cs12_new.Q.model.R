###################################################
### code chunk number 22: Cs12_new.Q.model
###################################################
for (i in 1:length(Z.models)) {
  if (i == 5) next # don't rerun panmictic
  for (Q.model in c("equalvarcov", "unconstrained")) {
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
  }
}


