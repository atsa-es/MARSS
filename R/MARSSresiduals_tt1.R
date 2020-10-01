MARSSresiduals.tt1 <- function(object, method=c("SS"), normalize = FALSE, silent=FALSE, fun.kf = c("MARSSkfas", "MARSSkfss") ) {
  # These are the residuals and their variance conditioned on the data up to time t-1

  ######################################
  # Set up variables
  MLEobj <- object
  if ( missing(fun.kf) ){
    fun.kf <- MLEobj$fun.kf
  }else{
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
  et <- st.et <- mar.st.et <- bchol.st.et <- matrix(NA, n+m, TT)
  model.et <- matrix(NA, n, TT)
  model.var.et <- array(0, dim = c(n, n, TT))
  var.et <- array(0, dim = c(n+m, n+m, TT))
  msg <- NULL
  
  #### list of time-varying parameters
  time.varying <- is.timevarying(MLEobj)
  
  kf <- MARSSkf(MLEobj)
  Kt.ok <- TRUE
  if(MLEobj$fun.kf=="MARSSkfas"){
    Kt <- try(MARSSkfss(MLEobj), silent=TRUE)
    if (inherits(Kt, "try-error") || !Kt$ok ){
      Kt <- array(0, dim=c(m, n, TT))
      Kt.ok <- FALSE
    }else{ Kt <- Kt$Kt }
    kf$Kt <- Kt
  }
  Ey <- MARSShatyt(MLEobj, only.kem=FALSE)
  Rt <- parmat(MLEobj, "R", t = 1)$R # returns matrix
  Ht <- parmat(MLEobj, "H", t = 1)$H
  Rt <- Ht %*% tcrossprod(Rt, Ht)
  Zt <- parmat(MLEobj, "Z", t = 1)$Z
  Qtp <- parmat(MLEobj, "Q", t = 2)$Q
  ######################################
  
  ######################################
  # Compute residuals
  if (method=="SS") {
    # We could set model.et to 0 where no data, but Kt will have a 0 column 
    # for any missing y.
    model.et <- Ey$ytt - fitted(MLEobj, type = "ytt1", output = "matrix") # model residuals
    et[1:n,] <- model.et
    
    cov.et <- matrix(0, n, m)
    # get the model residual variance
    for (t in 1:TT) {
      # model residuals
      if(time.varying$R) Rt <- parmat(MLEobj, "R", t = t)$R # returns matrix
      if(time.varying$H) Ht <- parmat(MLEobj, "H", t = t)$H
      if(time.varying$R || time.varying$H ) Rt <- Ht %*% tcrossprod(Rt, Ht)
      if(time.varying$Z) Zt <- parmat(MLEobj, "Z", t = t)$Z
      model.var.et[, , t] <- Rt + tcrossprod(Zt %*% kf$Vtt1[, , t], Zt)
    }
    # Then the states since that needs model var at t+1
    for (t in 1:TT){
      if(t < TT){
        Ktp <- sub3D(kf$Kt, t=t+1)
        et[(n + 1):(n + m), t] <- Ktp %*% model.et[, t+1]
        tmpvar.state.et <- tcrossprod(Ktp %*%model.var.et[, , t+1], Ktp)
      }else {
        tmpvar.state.et <- matrix(0, m, m)
      }
      var.et[1:n, , t] <- cbind(model.var.et[, , t], cov.et)
      var.et[(n + 1):(n + m), , t] <- cbind(t(cov.et), tmpvar.state.et)
      
      if (normalize) {
        Qpinv <- matrix(0, m, n + m)
        Rinv <- matrix(0, n, n + m)
        if (t < TT) {
          if(time.varying$Q) Qtp <- parmat(MLEobj, "Q", t = t + 1)$Q
          Qpinv[, (n + 1):(n + m)] <- psolve(t(pchol(Qtp)))
        }
        Rinv[, 1:n] <- psolve(t(pchol(Rt)))
        RQinv <- rbind(Rinv, Qpinv) # block diag matrix
        et[, t] <- RQinv %*% et[, t]
        var.et[, , t] <- tcrossprod(RQinv %*% var.et[, , t], RQinv)
      }
      
    }
  }
  ######################################
  
  ######################################
  # prepare standardized residuals
  for (t in 1:TT) {
    tmpvar <- sub3D(var.et, t = t)
    resids <- et[, t, drop = FALSE]
    # don't include values for resids if there is no residual (no data)
    # replace NAs with 0s
    is.miss <- c(is.na(y[, t]), rep(FALSE, m))
    resids[is.miss] <- 0
    
    tmpvar[abs(tmpvar) < sqrt(.Machine$double.eps)] <- 0
    
    # Marginal
    # inverse of diagonal of variance matrix for marginal standardization
    tmpvarinv <- try(psolve(makediag(takediag(tmpvar))), silent = TRUE)
    if (inherits(tmpvarinv, "try-error")) {
      mar.st.et[, t] <- NA
      msg <- c(msg, paste('MARSSresiduals.tt1 warning: the diagonal matrix of the variance of the residuals at t =", t, "is not invertible.  NAs returned for mar.residuals at t =", t, "\n'))
    }else{ # inv of diagonal ok, can compute marginal std residuals
      mar.st.et[, t] <- sqrt(tmpvarinv) %*% resids
      mar.st.et[is.miss, t] <- NA
    }
    
    # Block Cholesky
    # psolve and pchol deal with 0s on diagonal
    tmpchol <- try(pchol(tmpvar[(n + 1):(n + m), (n + 1):(n + m), drop = FALSE]), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      bchol.st.et[(n + 1):(n + m), t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tT warning: the variance of the state residuals at t =", t, "is not invertible.  NAs returned for bchol.std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
    } else {
      # chol() returns the upper triangle. We need to lower triangle to t()
      tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
      if (inherits(tmpcholinv, "try-error")) {
        bchol.st.et[(n + 1):(n + m), t] <- NA
        msg <- c(msg, paste("MARSSresiduals.tT warning: the variance of the state residuals at t =", t, "is not invertible.  NAs returned for bchol.std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      } else {
        bchol.st.et[(n + 1):(n + m), t] <- tmpcholinv %*% resids[(n + 1):(n + m), 1, drop = FALSE]
      }
    }
    
    # Cholesky
    # psolve and pchol deal with 0s on diagonal
    tmpchol <- try(pchol(tmpvar), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt1 warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
      next # got to next t since next inversion will fail
    }
    # chol() returns the upper triangle. We need to lower triangle
    tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
    if (inherits(tmpcholinv, "try-error")) {
      st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt1 warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      next # go to next t
    }
    # all ok, can compute Cholesky std.residuals
    st.et[, t] <- tmpcholinv %*% resids
    st.et[is.miss, t] <- NA
  }
  bchol.st.et[1:n, ] <- st.et[1:n, ] # because the upper right block of the lower tri is 0
  ######################################
  
  ######################################
  # the state.residual at the last time step is NA because it is x(T+1) - f(x(T)) and T+1 does not exist.  For the same reason, the var.residuals at TT will have NAs
  et[(n + 1):(n + m), TT] <- NA
  var.et[, (n + 1):(n + m), TT] <- NA
  var.et[(n + 1):(n + m), , TT] <- NA
  st.et[, TT] <- NA
  mar.st.et[(n + 1):(n + m), TT] <- NA  
  ######################################
  
  ######################################
  # et is the expected value of the residuals conditioned on y(1)-the observed data
  E.obs.v <- et[1:n,,drop=FALSE]
  var.obs.v <- array(0, dim = c(n, n, TT))
  # this will be 0 for observed data
  for( t in 1:TT) var.obs.v[,,t] <- Ey$Ott1[,,t] - tcrossprod(Ey$ytt1[,t])
  ######################################
  
  ######################################
  # the observed model residuals are data - E(data), so NA for missing data.
  model.et <- et[1:n, , drop = FALSE]
  model.et[is.na(y)] <- NA
  et[1:n,] <- model.et
  ######################################
  
  ######################################
  # add rownames
  Y.names <- attr(MLEobj$model, "Y.names")
  X.names <- attr(MLEobj$model, "X.names")
  rownames(et) <- rownames(st.et) <- rownames(mar.st.et) <- rownames(bchol.st.et) <- rownames(var.et) <- colnames(var.et) <- c(Y.names, X.names)
  rownames(E.obs.v) <- Y.names
  rownames(var.obs.v) <- colnames(var.obs.v) <- Y.names
  ######################################
  
  ######################################
  # output any warnings
  if(!is.null(msg) && object[["control"]][["trace"]] >= 0 & !silent) cat("MARSSresiduals.tt1 reported warnings. See msg element of returned residuals object.\n")
  ######################################
  
  if( !Kt.ok ){
    et[(n + 1):(n + m),] <- NA
    var.et[(n + 1):(n + m),(n+1):(n+m),] <- NA
    st.et[(n + 1):(n + m),] <- NA
    mar.st.et[(n + 1):(n + m),] <- NA
    bchol.st.et[(n + 1):(n + m),] <- NA
  }
  ret <- list(
    model.residuals = et[1:n, , drop = FALSE], 
    state.residuals = et[(n + 1):(n + m), , drop = FALSE], 
    residuals = et, 
    var.residuals = var.et, 
    std.residuals = st.et, 
    mar.residuals = mar.st.et, 
    bchol.residuals = bchol.st.et, 
    E.obs.residuals = E.obs.v, 
    var.obs.residuals = var.obs.v, 
    msg = msg)
  return(ret)
  
}
