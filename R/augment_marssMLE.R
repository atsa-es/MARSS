###############################################################################################################################################
#  augment method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is augment_marxss below
#  augment_dfa is just augment_marxss
##############################################################################################################################################
augment.marssMLE <- function(x, type = c("ytT", "xtT"),
                             interval = c("none", "confidence", "prediction"),
                             conf.level = 0.95,
                             form = attr(x[["model"]], "form")[1]) {
  ## Argument checking
  type <- match.arg(type)
  interval <- match.arg(interval)
  conditioning <- substr(type, 3, 3)
  if (conditioning != "T") {
    stop("augment.marssMLE: Only conditioning='T' allowed currently.", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level > 1 || conf.level < 0) {
    stop("augment.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
  }
  if (substr(type, 1, 1) == "y") type <- "observations"
  if (substr(type, 1, 1) == "x") type <- "states"
  ## End argument checking

  augment.fun <- paste("augment_", form, sep = "")
  tmp <- try(exists(augment.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    ret <- eval(call(augment.fun, x, type = type, interval = interval, conf.level = conf.level, conditioning = conditioning))
  } else {
    ret <- paste("No augment_", form[1], " is available.\n", sep = "")
  }
  ret
}

###############################################################################################################################################
#  augment method for class marssMLE form marxss
#  returns fitted values, residuals, std err of residuals and std residuals
#  the other forms use this
#  Set up to take other conditioning != T also
##############################################################################################################################################
augment_marxss <- function(x, type, interval, conf.level, conditioning) {
  # rotate means to rotate the Z matrix; this is used in DFA
  # but the user is allowed to do this for other cases also
  model <- x[["model"]]
  resids <- residuals(x)
  tmp <- apply(resids$var.residuals, 3, function(x) { takediag(x) })
  tmp[tmp<0 & abs(tmp)<.Machine$double.eps] <- 0
  se.resids <- sqrt(tmp)
  model.dims <- attr(model, "model.dims")
  data.dims <- model.dims[["y"]]
  nn <- data.dims[1]
  TT <- data.dims[2]
  alpha <- 1 - conf.level
  if (type == "observations") {
    data.names <- attr(model, "Y.names")
    fit.list <- fitted(x, type = "observations", interval=interval, conditioning=conditioning, conf.level=conf.level)
    if(interval=="none") fit.list <- list(.fitted=fit.list)
    ret <- data.frame(
      .rownames = rep(data.names, each = TT),
      t = rep(1:TT, nn),
      y = vec(t(model$data)),
      .fitted = vec(t(fit.list$.fitted))
    )
    if(interval=="confidence") ret <- cbind(ret,
                   .se.fit = vec(t(fit.list$.se.fit)),
                   .conf.low = vec(t(fit.list$.conf.low)),
                   .conf.up = vec(t(fit.list$.conf.up)))
    if(interval=="prediction") ret <- cbind(ret,
                                            .sd.y = vec(t(fit.list$.sd.y)),
                                            .lwr = vec(t(fit.list$.lwr)),
                                            .upr = vec(t(fit.list$.upr)))
    ret <- cbind(ret,
      .resids = vec(t(resids$model.residuals)),
      .sigma = vec(t(se.resids[1:nn, , drop = FALSE])),
      .std.resid = vec(t(resids$std.residuals[1:nn, , drop = FALSE]))
    )
  }
  if (type == "states") {
    # line up the residuals so that xtT(t) has residuals for xtT(t)-f(xtT(t-1))
    state.names <- attr(model, "X.names")
    state.dims <- attr(model, "model.dims")[["x"]]
    mm <- state.dims[1]
    state.se.resids <- cbind(NA, se.resids[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    state.resids <- cbind(NA, resids$state.residuals[, 1:(TT - 1), drop = FALSE])
    state.std.resids <- cbind(NA, resids$std.residuals[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    fit.list <- fitted(x, type = "states", interval=interval, conditioning=conditioning, conf.level=conf.level)
    if(interval=="none") fit.list <- list(.fitted=fit.list)
    ret <- data.frame(
      .rownames = rep(state.names, each = TT),
      t = rep(1:TT, mm),
      xtT = vec(t(x[["states"]])),
      .fitted = vec(t(fit.list$.fitted))
    )
    if(interval=="confidence") ret <- cbind(ret,
                                            .se.fit = vec(t(fit.list$.se.fit)),
                                            .conf.low = vec(t(fit.list$.conf.low)),
                                            .conf.up = vec(t(fit.list$.conf.up)))
    if(interval=="prediction") ret <- cbind(ret,
                                            .sd.x = vec(t(fit.list$.sd.x)),
                                            .lwr = vec(t(fit.list$.lwr)),
                                            .upr = vec(t(fit.list$.upr)))
    ret <- cbind(ret,
    .resids = vec(t(state.resids)),
    .sigma = vec(t(state.se.resids)),
    .std.resid = vec(t(state.std.resids))
    )
  }
  ret
}

augment_dfa <- function(x, type, interval, conf.level, conditioning) {
  ret <- augment_marxss(x, type = type, interval = interval, conf.level = conf.level, conditioning = conditioning)

  ret
}

augment_marss <- function(x, type, interval, conf.level, conditioning) {
  return(augment_marxss(x, type = type, interval = interval, conf.level = conf.level, conditioning = conditioning))
}