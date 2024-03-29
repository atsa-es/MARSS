###############################################################################################################################################
#  residuals method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is residuals_marxss below
#  residuals_dfa is just residuals_marxss
##############################################################################################################################################
residuals.marssMLE <- function(object, ...,
                               type = c("tt1", "tT", "tt"),
                               standardization = c("Cholesky", "marginal", "Block.Cholesky"),
                               form = attr(object[["model"]], "form")[1],
                               clean = TRUE) {
  ## Argument checking
  type <- match.arg(type)
  standardization <- match.arg(standardization)
  if (is.null(object[["par"]])) {
    stop("residuals.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (object[["convergence"]] == 54) {
    stop("residuals.marssMLE: MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", call. = FALSE)
  }
  

  resids.fun <- paste("residuals_", form, sep = "")
  tmp <- try(exists(resids.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    ret <- try(eval(call(resids.fun, object, type = type, standardization = standardization, clean = clean)), silent = TRUE)
    if (inherits(ret, "try-error")) stop("Stopped in residuals.marssMLE: MARSSresiduals() failed. \n", call. = FALSE)
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
residuals_marxss <- function(x, type, standardization, clean, ...) {
  # rotate means to rotate the Z matrix; this is used in DFA
  # but the user is allowed to do this for other cases also
  model <- x[["model"]]
  model.dims <- attr(model, "model.dims")
  data.dims <- model.dims[["y"]]
  mm <- model.dims[["x"]][1]
  nn <- data.dims[1]
  TT <- data.dims[2]
  resids <- MARSSresiduals(x, type = type, ...)
  tmp <- apply(resids$var.residuals, 3, function(x) {
    takediag(x)
  })
  tmp[tmp < 0 & abs(tmp) < sqrt(.Machine$double.eps)] <- 0
  if (is.null(dim(tmp))) tmp <- matrix(tmp, nrow = 1)
  se.resids <- sqrt(tmp)
  model.se.resids <- se.resids[1:nn, , drop = FALSE]
  model.se.resids[is.na(resids$model.residuals)] <- NA
  # First y
  type1 <- paste0("y", type)
  data.names <- attr(model, "Y.names")
  fit.list <- fitted.marssMLE(x, type = type1, interval = "none")
  ret <- data.frame(
    type = paste0("y", type),
    .rownames = fit.list$.rownames,
    name = "model",
    t = fit.list$t,
    value = fit.list$y,
    .fitted = fit.list$.fitted,
    stringsAsFactors = FALSE
  )
  ret <- cbind(ret,
    .resids = vec(t(resids$model.residuals)),
    .sigma = vec(t(model.se.resids))
  )
  if (standardization == "Cholesky") {
    ret <- cbind(ret, .std.resids = vec(t(resids$std.residuals[1:nn, , drop = FALSE])))
  }
  if (standardization == "marginal") {
    ret <- cbind(ret, .std.resids = vec(t(resids$mar.residuals[1:nn, , drop = FALSE])))
  }
  if (standardization == "Block.Cholesky") {
    ret <- cbind(ret, .std.resids = vec(t(resids$bchol.residuals[1:nn, , drop = FALSE])))
  }

  # Second x
  type1 <- paste0("x", type)
  state.names <- attr(model, "X.names")
  state.dims <- attr(model, "model.dims")[["x"]]
  mm <- state.dims[1]
  state.se.resids <- se.resids[(nn + 1):(nn + mm), , drop = FALSE]
  state.resids <- resids$state.residuals
  if (standardization == "Cholesky") {
    state.std.resids <- resids$std.residuals[(nn + 1):(nn + mm), , drop = FALSE]
  }
  if (standardization == "marginal") {
    state.std.resids <- resids$mar.residuals[(nn + 1):(nn + mm), , drop = FALSE]
  }
  if (standardization == "Block.Cholesky") {
    state.std.resids <- resids$bchol.residuals[(nn + 1):(nn + mm), , drop = FALSE]
  }
  if (type != "tt" && (type == "tT" || !clean)) {
    fit.list <- fitted.marssMLE(x, type = type1, interval = "none")
    # Note .resids(t) = .x(t+1) - .fitted(t+1); State residual time indexing is OFFSET
    # Although I could also offset the value and fitted, I decided not to because it would be confusing.
    # Instead added notes to help file.
    # loc <- rep(c(2:TT,NA), mm)+TT*rep(0:(mm-1), each=TT)
    ret2 <- data.frame(
      type = paste0("x", type),
      .rownames = fit.list$.rownames,
      name = "state",
      t = fit.list$t,
      value = fit.list$.x, # t indexing is offset relative to the state resids
      .fitted = fit.list$.fitted, # t indexing is offset relative to the state resids
      .resids = vec(t(state.resids)),
      .sigma = vec(t(state.se.resids)),
      .std.resids = vec(t(state.std.resids)),
      stringsAsFactors = FALSE
    )
    ret <- rbind(ret, ret2)
  }

  class(ret) <- c("marssResiduals", "data.frame")
  attr(ret, "standardization") <- standardization
  attr(ret, "residual.type") <- type
  attr(ret, "msg") <- resids$msg
  ret
}

residuals_dfa <- function(x, type, standardization, clean, ...) {
  return(residuals_marxss(x, type, standardization, clean))
}

residuals_marss <- function(x, type, standardization, clean, ...) {
  return(residuals_marxss(x, type, standardization, clean))
}
