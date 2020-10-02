MARSSresiduals.tt <- function(object, method = c("SS"), normalize = FALSE, silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  # These are the residuals and their variance conditioned on the data up to time t
  # state residuals do not exist for this case

  ######################################
  # Set up variables
  MLEobj <- object
  if (missing(fun.kf)) {
    fun.kf <- MLEobj$fun.kf
  } else {
    MLEobj$fun.kf <- fun.kf
    # to ensure that MARSSkf, MARSShatyt and fitted() use the spec'd fun
  }
  method <- match.arg(method)
  model.dims <- attr(MLEobj$marss, "model.dims")
  TT <- model.dims[["x"]][2]
  m <- model.dims[["x"]][1]
  n <- model.dims[["y"]][1]
  y <- MLEobj$marss$data
  # set up holders
  et <- st.et <- mar.st.et <- matrix(NA, n + m, TT)
  model.et <- model.st.et <- model.mar.st.et <- matrix(NA, n, TT)
  model.var.et <- array(0, dim = c(n, n, TT))
  var.et <- array(0, dim = c(n + m, n + m, TT))
  msg <- NULL

  #### list of time-varying parameters
  time.varying <- is.timevarying(MLEobj)

  kf <- MARSSkf(MLEobj)
  Ey <- MARSShatyt(MLEobj, only.kem = FALSE)
  Rt <- parmat(MLEobj, "R", t = 1)$R # returns matrix
  Ht <- parmat(MLEobj, "H", t = 1)$H
  Rt <- Ht %*% tcrossprod(Rt, Ht)
  Zt <- parmat(MLEobj, "Z", t = 1)$Z
  Qtp <- parmat(MLEobj, "Q", t = 2)$Q

  if (method == "SS") {
    # We could set model.et to 0 where no data, but Kt will have a 0 column
    # for any missing y.
    model.et <- Ey$ytt - fitted(MLEobj, type = "ytt", output = "matrix") # model residuals
    et[1:n, ] <- model.et

    cov.et <- matrix(NA, n, m)
    for (t in 1:TT) {
      # model residuals
      if (time.varying$R) Rt <- parmat(MLEobj, "R", t = t)$R # returns matrix
      if (time.varying$H) Ht <- parmat(MLEobj, "H", t = t)$H
      if (time.varying$R || time.varying$H) Rt <- Ht %*% tcrossprod(Rt, Ht)
      if (time.varying$Z) Zt <- parmat(MLEobj, "Z", t = t)$Z
      # compute the variance of the residuals and state.et
      St <- Ey$yxtt[, , t] - tcrossprod(Ey$ytt[, t, drop = FALSE], kf$xtt[, t, drop = FALSE])
      model.var.et[, , t] <- Rt - tcrossprod(Zt %*% kf$Vtt[, , t], Zt) + tcrossprod(St, Zt) + tcrossprod(Zt, St)

      tmpvar.state.et <- matrix(NA, m, m)
      var.et[1:n, , t] <- cbind(model.var.et[, , t], cov.et)
      var.et[(n + 1):(n + m), , t] <- cbind(t(cov.et), tmpvar.state.et)

      if (normalize) {
        Rinv <- psolve(t(pchol(Rt)))
        model.et[, t] <- Rinv %*% model.et[, t]
        et[1:n, t] <- model.et[, t]
        model.var.et[, , t] <- tcrossprod(Rinv %*% model.var.et[, , t], Rinv)
        var.et[1:n, , t] <- cbind(model.var.et[, , t], cov.et)
      }
    }
  }

  # prepare standardized residuals; only model residuals since state residuals don't exist for this case
  for (t in 1:TT) {
    tmpvar <- sub3D(model.var.et, t = t)
    resids <- model.et[, t, drop = FALSE]
    # don't include values for resids if there is no residual (no data)
    # replace NAs with 0s
    is.miss <- is.na(y[, t])
    resids[is.miss] <- 0

    tmpvar[abs(tmpvar) < sqrt(.Machine$double.eps)] <- 0

    # inverse of diagonal of variance matrix for marginal standardization
    tmpvarinv <- try(psolve(makediag(takediag(tmpvar))), silent = TRUE)
    if (inherits(tmpvarinv, "try-error")) {
      model.mar.st.et[, t] <- NA
      msg <- c(msg, paste('MARSSresiduals.tt warning: the diagonal matrix of the variance of the residuals at t =", t, "is not invertible.  NAs returned for mar.residuals at t =", t, "\n'))
    } else { # inverse of diagonal ok
      model.mar.st.et[, t] <- sqrt(tmpvarinv) %*% resids
      model.mar.st.et[is.miss, t] <- NA
    }

    # psolve and pchol deal with 0s on diagonal
    # wrapped in try to prevent crashing if inversion not possible
    tmpchol <- try(pchol(tmpvar), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      model.st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
      next
    }
    # chol is ok
    # chol() returns the upper triangle. We need to lower triangle
    tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
    if (inherits(tmpcholinv, "try-error")) {
      model.st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      next
    }
    # both chol and inverse are ok
    model.st.et[, t] <- tmpcholinv %*% resids
    model.st.et[is.miss, t] <- NA
  }
  st.et[1:n, ] <- model.st.et
  mar.st.et[1:n, ] <- model.mar.st.et

  # Do not need to NA out the state residuals since they are all NA in this case

  # et is the expected value of the residuals conditioned on y(1,t)-the observed data up to time t
  E.obs.v <- et[1:n, , drop = FALSE]
  var.obs.v <- array(0, dim = c(n, n, TT))
  # this will be 0 for observed data
  for (t in 1:TT) var.obs.v[, , t] <- Ey$Ott[, , t] - tcrossprod(Ey$ytt[, t])

  # the observed model residuals are data - E(data), so NA for missing data.
  model.et <- et[1:n, , drop = FALSE]
  model.et[is.na(y)] <- NA
  et[1:n, ] <- model.et

  # add rownames
  Y.names <- attr(MLEobj$model, "Y.names")
  X.names <- attr(MLEobj$model, "X.names")
  rownames(et) <- rownames(st.et) <- rownames(mar.st.et) <- rownames(var.et) <- colnames(var.et) <- c(Y.names, X.names)
  rownames(E.obs.v) <- Y.names
  rownames(var.obs.v) <- colnames(var.obs.v) <- Y.names

  # output any warnings
  if (!is.null(msg) && object[["control"]][["trace"]] >= 0 & !silent) cat("MARSSresiduals.tt reported warnings. See msg element of returned residuals object.\n")

  return(list(
    model.residuals = et[1:n, , drop = FALSE],
    state.residuals = et[(n + 1):(n + m), , drop = FALSE],
    residuals = et,
    var.residuals = var.et,
    std.residuals = st.et,
    mar.residuals = mar.st.et,
    bchol.residuals = st.et,
    E.obs.residuals = E.obs.v,
    var.obs.residuals = var.obs.v,
    msg = msg
  ))
}
