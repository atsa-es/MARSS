###############################################################################################################################################
#  fitted method for class marssMLE.
##############################################################################################################################################
fitted.marssMLE <- function(object, ...,
                            type = c("observations", "states", "y", "x"),
                            conditioning = c("T", "t", "t-1"),
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95,
                            ...) {
  type <- match.arg(type)
  if (type == "y") type <- "observations"
  if (type == "x") type <- "states"
  conditioning <- match.arg(conditioning)
  interval <- match.arg(interval)
  MLEobj <- object
  if (is.null(MLEobj[["par"]])) {
    stop("fitted.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (!missing(level) && interval!="none") warning("fitted.marssMLE: interval='none' thus level is ignored.", call. = FALSE)
  if (!is.numeric(level) && level > 0 && level < 1) stop("fitted.marssMLE: level must be numeric and between 0 and 1.", call. = FALSE)
  alpha <- 1-level
  extras <- list()
  if (!missing(...)) {
    extras <- list(...)
    if ("one.step.ahead" %in% names(extras)) stop("fitted.marssMLE: Use conditioning='t-1' instead of one.step.ahead=TRUE.", call. = FALSE)
  }
  if (conditioning == "t" && type == "states") {
    stop("fitted.marssMLE: Only conditioning = 'T' and 't-1' allowed for type='states'. ")
  }

  # need the model dims in marss form with c in U and d in A
  model.dims <- attr(MLEobj[["marss"]], "model.dims")
  TT <- model.dims[["x"]][2]
  n <- model.dims[["y"]][1]
  m <- model.dims[["x"]][1]

  if (type == "observations") {
    if (conditioning == "T") hatxt <- MLEobj[["states"]]
    if (conditioning == "t-1") hatxt <- MARSSkf(MLEobj)[["xtt1"]]
    if (conditioning == "t") hatxt <- MARSSkf(MLEobj)[["xtt"]]
    if(return.se){
      if (conditioning == "T") hatVt <- MARSSkf(MLEobj)[["VtT"]]
      if (conditioning == "t-1") hatVt <- MARSSkf(MLEobj)[["Vtt1"]]
      if (conditioning == "t") hatVt <- MARSSkf(MLEobj)[["Vtt"]]
    }
    Z.time.varying <- model.dims[["Z"]][3] != 1
    A.time.varying <- model.dims[["A"]][3] != 1
    R.time.varying <- model.dims[["R"]][3] != 1
    
    val <- matrix(NA, n, TT)
    rownames(val) <- attr(MLEobj$marss, "Y.names")
    if(return.se) se <- val

    Zt <- parmat(MLEobj, "Z", t = 1)$Z
    At <- parmat(MLEobj, "A", t = 1)$A
    Rt <- parmat(MLEobj, "R", t = 1)$R
    
    for (t in 1:TT) {
      # parmat returns marss form
      if (Z.time.varying) Zt <- parmat(MLEobj, "Z", t = t)$Z
      if (A.time.varying) At <- parmat(MLEobj, "A", t = t)$A
      val[, t] <- Zt %*% hatxt[, t, drop = FALSE] + At
      if(interval=="confidence") se[, t] <- takediag(Zt %*% tcrossprod( hatVt[, t, drop = FALSE], Zt))
      if(interval=="prediction"){
        if (R.time.varying) Rt <- parmat(MLEobj, "R", t = t)$R
        se[, t] <- takediag(Zt %*% tcrossprod( hatVt[, t, drop = FALSE], Zt)) + Rt
      }
    }
  }

  if (type == "states") {
    if (conditioning == "T") hatxt <- MLEobj[["states"]]
    if (conditioning == "t-1") hatxt <- MARSSkf(MLEobj)[["xtt"]]
    if(return.se){
      if (conditioning == "T") hatVt <- MARSSkf(MLEobj)[["VtT"]]
      if (conditioning == "t-1") hatVt <- MARSSkf(MLEobj)[["Vtt"]]
    }
    
    B.time.varying <- model.dims[["B"]][3] != 1
    U.time.varying <- model.dims[["U"]][3] != 1
    Q.time.varying <- model.dims[["Q"]][3] != 1
    
    val <- matrix(NA, m, TT)
    rownames(val) <- attr(MLEobj$marss, "X.names")
    if(return.se) se <- val
    
    x0 <- coef(MLEobj, type = "matrix")[["x0"]]
    if(return.se) V0 <- coef(MLEobj, type = "matrix")[["V0"]]
    Bt <- parmat(MLEobj, "B", t = 1)[["B"]]
    Ut <- parmat(MLEobj, "U", t = 1)[["U"]]
    Qt <- parmat(MLEobj, "Q", t = 1)[["Q"]]
    if (MLEobj$model$tinitx == 0){
      val[, 1] <- Bt %*% x0 + Ut
      if (interval=="confidence") se[, 1] <- takediag(Bt %*% tcrossprod(V0, Bt))
      if (interval=="prediction") se[, 1] <- takediag(Bt %*% tcrossprod(V0, Bt)) + Qt
    }
    if (MLEobj$model$tinitx == 1){
      val[, 1] <- x0
      if(interval != "none") se[, 1] <- takediag(V0)
    }
    for (t in 2:TT) {
      if (B.time.varying) Bt <- parmat(MLEobj, "B", t = t)[["B"]]
      if (U.time.varying) Ut <- parmat(MLEobj, "U", t = t)[["U"]]
      val[, t] <- Bt %*% hatxt[, t - 1, drop = FALSE] + Ut
      if (interval=="confidence") se[, t] <- takediag(Bt %*% tcrossprod( hatVt[, t-1, drop = FALSE], Bt))
      if (interval=="prediction"){
        if (Q.time.varying) Qt <- parmat(MLEobj, "Q", t = t)[["Q"]]
        se[, 1] <- takediag(Bt %*% tcrossprod(V0, Bt)) + Qt
      } 
    }
  }
  if (interval=="none") return(val)
  if (interval=="confidence"){
    return(list(
      .fitted=val, 
      .se.fit=se,
      .conf.low = val + qnorm(alpha/2) * se,
      .conf.up = val + qnorm(1- alpha/2) * se
      ))
  }
  # Then interval == "prediction"
  if (type == "states"){
    return(list(
      .fitted=val, 
      .sd.x=se, 
      .lwr = val + qnorm(alpha/2) * se,
      .upr = val + qnorm(1- alpha/2) * se
    ))
  }
  if (type == "observations"){
    return(list(
      .fitted=val,
      .sd.y=se,
      .lwr = val + qnorm(alpha/2) * se,
      .upr = val + qnorm(1- alpha/2) * se
    ))
  }
} # end of fitted.marssMLE
