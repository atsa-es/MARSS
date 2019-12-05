###############################################################################################################################################
#  augment method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is augment_marxss below
#  augment_dfa uses the optional rotate argument
##############################################################################################################################################
augment.marssMLE <- function(x, type.predict = c("observations", "states"),
                             interval = c("none", "confidence"),
                             conf.level = 0.95,
                             form = attr(x[["model"]], "form")) {
  ## Argument checking
  type.predict <- match.arg(type.predict)
  interval <- match.arg(interval)
  if (!is.numeric(conf.level) || length(conf.level)>1 || conf.level > 1 || conf.level < 0) stop("augment.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
  ## End argument checking

  augment.fun <- paste("augment_", form[1], sep = "")
  tmp <- try(exists(augment.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    ret <- eval(call(augment.fun, x, type.predict = type.predict, interval = interval, conf.level = conf.level))
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
augment_marxss <- function(x, type.predict, interval, conf.level) {
  # rotate means to rotate the Z matrix; this is used in DFA
  # but the user is allowed to do this for other cases also
  model <- x[["model"]]
  resids <- residuals(x)
  se.resids <- sqrt(apply(resids$var.residuals, 3, function(x) {
    takediag(x)
  }))
  data.dims <- attr(model, "model.dims")[["y"]]
  nn <- data.dims[1]
  TT <- data.dims[2]
  alpha <- 1 - conf.level
  if (type.predict == "observations") {
    data.names <- attr(model, "Y.names")
    ret <- data.frame(
      .rownames = rep(data.names, each = TT),
      t = rep(1:TT, nn),
      y = vec(t(model$data)),
      .fitted = vec(t(fitted(x, type = "observations"))),
      .se.fit = vec(t(se.resids[1:nn, , drop = FALSE])),
      .resids = vec(t(resids$model.residuals)),
      .std.resid = vec(t(resids$std.residuals[1:nn, , drop = FALSE]))
    )
    if (interval == "confidence") {
      ret <- cbind(ret,
        .conf.low = qnorm(alpha / 2) * ret$.se.fit + ret$.fitted,
        .conf.up = qnorm(1 - alpha / 2) * ret$.se.fit + ret$.fitted
      )
    }
  }
  if (type.predict == "states") {
    # line up the residuals so that xtT(t) has residuals for xtT(t)-f(xtT(t-1))
    state.names <- attr(model, "X.names")
    state.dims <- attr(model, "model.dims")[["x"]]
    mm <- state.dims[1]
    state.se.resids <- cbind(NA, se.resids[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    state.resids <- cbind(NA, resids$state.residuals[, 1:(TT - 1), drop = FALSE])
    state.std.resids <- cbind(NA, resids$std.residuals[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    ret <- data.frame(
      .rownames = rep(state.names, each = TT),
      t = rep(1:TT, mm),
      xtT = vec(t(x[["states"]])),
      .fitted = vec(t(fitted(x, type = "states"))),
      .se.fit = vec(t(state.se.resids)),
      .resids = vec(t(state.resids)),
      .std.resid = vec(t(state.std.resids))
    )
    if (interval == "confidence") {
      ret <- cbind(ret,
        .conf.low = qnorm(alpha / 2) * ret$.se.fit + ret$.fitted,
        .conf.up = qnorm(1 - alpha / 2) * ret$.se.fit + ret$.fitted
      )
    }
  }
  ret
}

augment_dfa <- function(x, type.predict, interval, conf.level) {
  ret <- augment_marxss(x, type.predict = type.predict, interval = interval, conf.level = conf.level)

  ret
}
