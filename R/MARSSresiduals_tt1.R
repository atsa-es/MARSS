MARSSresiduals.tt1 <- function(object, ..., method=c("SS"), normalize = FALSE) {
  # These are the residuals and their variance conditioned on the data up to time t-1
  MLEobj <- object
  method <- match.arg(method)
  model.dims <- attr(MLEobj$marss, "model.dims")
  TT <- model.dims[["x"]][2]
  m <- model.dims[["x"]][1]
  n <- model.dims[["y"]][1]
  y <- MLEobj$marss$data
  # set up holders
  et <- st.et <- mar.st.et <- matrix(0, n , TT)
  var.et <- array(0, dim = c(n, n, TT))

  #### make a list of time-varying parameters
  time.varying <- list()
  for (elem in attr(MLEobj[["marss"]], "par.names")) {
    if (model.dims[[elem]][3] == 1) {
      time.varying[[elem]] <- FALSE
    } else {
      time.varying[[elem]] <- TRUE
    }
  }

  kf <- MARSSkfss(MLEobj)
  Ey <- MARSShatyt(MLEobj, only.kem=FALSE)
  Rt <- parmat(MLEobj, "R", t = 1)$R # returns matrix
  Ht <- parmat(MLEobj, "H", t = 1)$H
  Rt <- Ht %*% Rt %*% t(Ht)
  Zt <- parmat(MLEobj, "Z", t = 1)$Z

  if (method=="SS") {
    # model.et will be 0 where no data E(y)-modeled(y)
    et <- y - fitted(MLEobj, type = "observations", conditioning = "t-1", output = "matrix") # model residuals

    for (t in 1:TT) {
      # model residuals
      if(time.varying$R) Rt <- parmat(MLEobj, "R", t = t)$R # returns matrix
      if(time.varying$H) Ht <- parmat(MLEobj, "H", t = t)$H
      if(time.varying$R || time.varying$H ) Rt <- Ht %*% Rt %*% t(Ht)
      if(time.varying$Z) Zt <- parmat(MLEobj, "Z", t = t)$Z

      # model.et and state.et defined outside for loop
      var.et[, , t] <- Rt + Zt %*% kf$Vtt1[, , t] %*% t(Zt)

      if (normalize) {
        Rinv <- matrix(0, n, n + m)
        Rinv <- psolve(t(pchol(Rt)))
        et[, t] <- Rinv %*% et[, t]
        var.et[, , t] <- Rinv %*% var.et[, , t] %*% t(Rinv)
      }
    }
  }

  # prepare cholesky standardized residuals
  for (t in 1:TT) {
    tmpvar <- sub3D(var.et, t = t)
    resids <- et[, t, drop = FALSE]
    # don't includ values for resids if there is no residual (no data)
    is.miss <- is.na(y[, t])
    resids[is.miss] <- 0

    tmpvar[abs(tmpvar) < sqrt(.Machine$double.eps)] <- 0

    # psolve and pchol deal with 0s on diagonal
    # wrapped in try to prevent crashing if inversion not possible
    tmpchol <- try(pchol(tmpvar), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      st.et[, t] <- NA
      cat(paste("warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
      next
    }
    tmpcholinv <- try(psolve(tmpchol), silent = TRUE)
    if (inherits(tmpcholinv, "try-error")) {
      st.et[, t] <- NA
      cat(paste("warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, "\n"))
      next
    }
    # inverse of diagonal of variance matrix for marginal standardization
    tmpvarinv <- try(psolve(makediag(takediag(tmpvar))), silent = TRUE)
    if (inherits(tmpvarinv, "try-error")) {
      mar.st.et[, t] <- NA
      cat(paste("warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, "\n"))
      next
    }
    st.et[, t] <- tmpcholinv %*% resids
    st.et[is.miss, t] <- NA
    mar.st.et[, t] <- tmpvarinv %*% resids
    mar.st.et[is.miss, t] <- NA
  }

  # add rownames
  Y.names <- attr(MLEobj$model, "Y.names")
  rownames(et) <- rownames(st.et) <- rownames(var.et) <- colnames(var.et) <- Y.names

  return(list(residuals = et, std.residuals = st.et, mar.residuals = mar.st.et, var.residuals = var.et))
}
