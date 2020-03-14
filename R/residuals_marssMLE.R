###############################################################################################################################################
#  residuals method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is residuals_marxss below
#  residuals_dfa is just residuals_marxss
##############################################################################################################################################
residuals.marssMLE <- function(x,
                               type=c("smoothations", "innovations"),
                               standardization=c("Cholesky", "marginal"),
                               form = attr(x[["model"]], "form")[1]) {
  ## Argument checking
  type <- match.arg(type)
  if (type == "smoothations") type <- "tT"
  if (type == "innovations") type <- "tt1"
  standardization <- match.arg(standardization)

  resids.fun <- paste("residuals_", form, sep = "")
  tmp <- try(exists(resids.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    ret <- eval(call(resids.fun, x, type = type, standardization = standardization))
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
residuals_marxss <- function(x, type, standardization) {
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
  # First y
  type1 <- paste0("y", type)
    data.names <- attr(model, "Y.names")
    fit.list <- fitted(x, type = type1, interval="none")
    ret <- data.frame(
      .rownames = fit.list$.rownames,
      type = "model",
      t = fit.list$t,
      value = fit.list$y,
      .fitted = fit.list$.fitted
    )
    ret <- cbind(ret,
      .resids = vec(t(resids$model.residuals)),
      .sigma = vec(t(model.se.resids))
      )
    if(standardization=="Cholesky") 
      ret <- cbind(ret, .std.resid = vec(t(resids$std.residuals[1:nn, , drop = FALSE])))
    if(standardization=="marginal") 
      ret <- cbind(ret, .std.resid = vec(t(resids$mar.residuals[1:nn, , drop = FALSE])))
    
  # Second x
    type1 <- paste0("x", type)
    # line up the residuals so that xtT(t) has residuals for xtT(t)-f(xtT(t-1))
    state.names <- attr(model, "X.names")
    state.dims <- attr(model, "model.dims")[["x"]]
    mm <- state.dims[1]
    state.se.resids <- cbind(NA, se.resids[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    state.resids <- cbind(NA, resids$state.residuals[, 1:(TT - 1), drop = FALSE])
    if(standardization=="Cholesky") 
      state.std.resids <- cbind(NA, resids$std.residuals[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    if(standardization=="marginal") 
      state.std.resids <- cbind(NA, resids$mar.residuals[(nn + 1):(nn + mm), 1:(TT - 1), drop = FALSE])
    fit.list <- fitted(x, type = type1, interval="none")
    if (type1=="xtt1") fit.list$value <- fit.list$xtt
    if (type1=="xtT") fit.list$value <- fit.list$xtT
    ret2 <- data.frame(
      .rownames = fit.list$.rownames,
      type = "state",
      t = fit.list$t,
      value = fit.list$value,
      .fitted = fit.list$.fitted,
      .resids = vec(t(state.resids)),
      .sigma = vec(t(state.se.resids)),
      .std.resid = vec(t(state.std.resids))
    )
    ret <- rbind(ret, ret2)
  
  class(ret) <- c("marssResiduals", "tbl_df", "tbl", "data.frame")
  attr(ret, "standardization") <- standardization
  attr(ret, "residual.type") <- type
  ret
}

residuals_dfa <- function(x, type, standardization) {
  return(residuals_marxss(x, type, standardization))
}

residuals_marss <- function(x, type, standardization) {
  return(residuals_marxss(x, type, standardization))
}