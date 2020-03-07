###############################################################################################################################################
#  augment method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is augment_marxss below
#  augment_dfa is just augment_marxss
##############################################################################################################################################
augment.marssMLE <- function(x, type = c("ytT", "xtT"),
                             interval = c("none", "confidence", "prediction"),
                             conf.level = 0.95) {
  ## Argument checking
  type <- match.arg(type)
  interval <- match.arg(interval)
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level > 1 || conf.level < 0) {
    stop("augment.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
  }
  form = attr(x[["model"]], "form")[1]
  augment.fun <- paste("augment_", form, sep = "")
  tmp <- try(exists(augment.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    ret <- eval(call(augment.fun, x, type = type, interval = interval, conf.level = conf.level))
  } else {
    ret <- paste("No augment_", form[1], " is available.\n", sep = "")
  }
  ret
}

###############################################################################################################################################
#  augment method for class marssMLE form marxss
#  returns fitted values, residuals, std err of residuals and std residuals
#  the other forms use this
##############################################################################################################################################
augment_marxss <- function(x, type, interval, conf.level) {
  # rotate means to rotate the Z matrix; this is used in DFA
  # but the user is allowed to do this for other cases also
  model <- x[["model"]]
  model.dims <- attr(model, "model.dims")
  data.dims <- model.dims[["y"]]
  nn <- data.dims[1]
  TT <- data.dims[2]
  resids <- residuals(x)
  tmp <- apply(resids$var.residuals, 3, function(x) { takediag(x) })
  tmp[tmp<0 & abs(tmp)<sqrt(.Machine$double.eps)] <- 0
  se.resids <- sqrt(tmp)
  model.se.resids <-  se.resids[1:nn,,drop=FALSE]
  model.se.resids[is.na(resids$model.residuals)] <- NA
  alpha <- 1 - conf.level
  if (substr(type, 1, 1) == "y") {
    data.names <- attr(model, "Y.names")
    fit.list <- fitted(x, type = type, interval=interval, conf.level=conf.level)
    ret <- data.frame(
      .rownames = fit.list$.rownames,
      t = fit.list$t,
      y = fit.list$y,
      .fitted = fit.list$.fitted
    )
    if(interval=="confidence") ret <- cbind(ret,
                   .se.fit = fit.list$.se.fit,
                   .conf.low = fit.list$.conf.low,
                   .conf.up = fit.list$.conf.up)
    if(interval=="prediction") ret <- cbind(ret,
                                            .sd.y = fit.list$.sd.y,
                                            .lwr = fit.list$.lwr,
                                            .upr = fit.list$.upr)
    ret <- cbind(ret,
      .resids = vec(t(resids$model.residuals)),
      .sigma = vec(t(model.se.resids)),
      .std.resid = vec(t(resids$std.residuals[1:nn, , drop = FALSE])),
      .mar.resid = vec(t(resids$mar.residuals[1:nn, , drop = FALSE]))
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
    fit.list <- fitted(x, type = type, interval=interval, conf.level=conf.level)
    ret <- data.frame(
      .rownames = fit.list$.rownames,
      t = fit.list$t,
      xtT = fit.list$xtT,
      .fitted = fit.list$.fitted
    )
    if(interval=="confidence") ret <- cbind(ret,
                                            .se.fit = fit.list$.se.fit,
                                            .conf.low = fit.list$.conf.low,
                                            .conf.up = fit.list$.conf.up)
    if(interval=="prediction") ret <- cbind(ret,
                                            .sd.x = fit.list$.sd.x,
                                            .lwr = fit.list$.lwr,
                                            .upr = fit.list$.upr)
    ret <- cbind(ret,
    .resids = vec(t(state.resids)),
    .sigma = vec(t(state.se.resids)),
    .std.resid = vec(t(state.std.resids)),
    .mar.resid = vec(t(state.mar.resids))
    )
  }
  ret
}

augment_dfa <- function(x, type, interval, conf.level) {
  ret <- augment_marxss(x, type = type, interval = interval, conf.level = conf.level)

  ret
}

augment_marss <- function(x, type, interval, conf.level) {
  return(augment_marxss(x, type = type, interval = interval, conf.level = conf.level))
}