###############################################################################################################################################
#  residuals method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is residuals_marxss below
#  residuals_dfa is just residuals_marxss
##############################################################################################################################################
residuals.marssMLE <- function(x, type = c("ytT", "xtT"),
                             form = attr(x[["model"]], "form")[1]) {
  ## Argument checking
  type <- match.arg(type)

  resids.fun <- paste("residuals_", form, sep = "")
  tmp <- try(exists(resids.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    ret <- eval(call(resids.fun, x, type = type))
  } else {
    ret <- paste("No residuals_", form[1], " is available.\n", sep = "")
  }
  ret
}

###############################################################################################################################################
#  residuals method for class marssMLE form marxss
#  returns fitted values, residuals, std err of residuals and std residuals
#  the other forms use this
##############################################################################################################################################
residuals_marxss <- function(x, type) {
  # rotate means to rotate the Z matrix; this is used in DFA
  # but the user is allowed to do this for other cases also
  model <- x[["model"]]
  model.dims <- attr(model, "model.dims")
  data.dims <- model.dims[["y"]]
  nn <- data.dims[1]
  TT <- data.dims[2]
  resids <- MARSSresiduals(x)
  tmp <- apply(resids$var.residuals, 3, function(x) { takediag(x) })
  tmp[tmp<0 & abs(tmp)<sqrt(.Machine$double.eps)] <- 0
  se.resids <- sqrt(tmp)
  model.se.resids <-  se.resids[1:nn,,drop=FALSE]
  model.se.resids[is.na(resids$model.residuals)] <- NA
  if (substr(type, 1, 1) == "y") {
    data.names <- attr(model, "Y.names")
    fit.list <- fitted(x, type = type, interval="none")
    ret <- data.frame(
      .rownames = fit.list$.rownames,
      t = fit.list$t,
      y = fit.list$y,
      .fitted = fit.list$.fitted
    )
    ret <- cbind(ret,
      .resids = vec(t(resids$model.residuals)),
      .sigma = vec(t(model.se.resids)),
      .std.resid = vec(t(resids$std.residuals[1:nn, , drop = FALSE]))
    )
  }
  if (substr(type, 1, 1) == "x") {
    # line up the residuals so that xtT(t) has residuals for xtT(t)-f(xtT(t-1))
    state.names <- attr(model, "X.names")
    state.dims <- attr(model, "model.dims")[["x"]]
    mm <- state.dims[1]
    state.se.resids <- cbind(NA, se.resids[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    state.resids <- cbind(NA, resids$state.residuals[, 1:(TT - 1), drop = FALSE])
    state.std.resids <- cbind(NA, resids$std.residuals[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    state.mar.resids <- cbind(NA, resids$mar.residuals[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    fit.list <- fitted(x, type = type, interval="none")
    ret <- data.frame(
      .rownames = fit.list$.rownames,
      t = fit.list$t,
      xtT = fit.list$xtT,
      .fitted = fit.list$.fitted
    )
    ret <- cbind(ret,
    .resids = vec(t(state.resids)),
    .sigma = vec(t(state.se.resids)),
    .std.resid = vec(t(state.std.resids))
    )
  }
  ret
}

residuals_dfa <- function(x, type) {
  return(residuals_marxss(x, type = type))
}

residuals_marss <- function(x, type) {
  return(residuals_marxss(x, type = type))
}