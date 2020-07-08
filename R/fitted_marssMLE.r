###############################################################################################################################################
#  fitted method for class marssMLE.
##############################################################################################################################################
fitted.marssMLE <- function(object, ...,
                            type = c("ytT", "xtT", "ytt", "ytt1", "xtt1"),
                            interval = c("none", "confidence", "prediction"),
                            conf.level = 0.95,
                            output = c("tibble", "matrix")) {
  type <- match.arg(type)
  output <- match.arg(output)
  interval <- match.arg(interval)
  conditioning <- substring(type, 3)
  type <- substr(type, 1, 1)
  MLEobj <- object
  if (is.null(MLEobj[["par"]])) {
    stop("fitted.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (interval != "none" && (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level > 1 || conf.level < 0))
    stop("fitted.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
  alpha <- 1-conf.level
  extras <- list()
  if (!missing(...)) {
    extras <- list(...)
    if ("one.step.ahead" %in% names(extras)) stop("fitted.marssMLE: Use type='ytt1' or 'xtt1' instead of one.step.ahead=TRUE.", call. = FALSE)
  }
  
  # need the model dims in marss form with c in U and d in A
  model.dims <- attr(MLEobj[["marss"]], "model.dims")
  TT <- model.dims[["x"]][2]
  n <- model.dims[["y"]][1]
  m <- model.dims[["x"]][1]
  
  if (type == "y") {
    if (conditioning == "T") hatxt <- MARSSkf(MLEobj)[["xtT"]]
    if (conditioning == "t1") hatxt <- MARSSkf(MLEobj)[["xtt1"]]
    if (conditioning == "t") hatxt <- MARSSkfss(MLEobj)[["xtt"]]
    if(interval!="none"){
      if (conditioning == "T") hatVt <- MARSSkf(MLEobj)[["VtT"]]
      if (conditioning == "t1") hatVt <- MARSSkf(MLEobj)[["Vtt1"]]
      if (conditioning == "t") hatVt <- MARSSkfss(MLEobj)[["Vtt"]]
    }
    Z.time.varying <- model.dims[["Z"]][3] != 1
    A.time.varying <- model.dims[["A"]][3] != 1
    R.time.varying <- model.dims[["R"]][3] != 1
    
    val <- matrix(NA, n, TT)
    rownames(val) <- attr(MLEobj$marss, "Y.names")
    if(interval!="none") se <- val
    
    Zt <- parmat(MLEobj, "Z", t = 1)$Z
    At <- parmat(MLEobj, "A", t = 1)$A
    Rt <- parmat(MLEobj, "R", t = 1)$R
    
    for (t in 1:TT) {
      # parmat returns marss form
      if (Z.time.varying) Zt <- parmat(MLEobj, "Z", t = t)$Z
      if (A.time.varying) At <- parmat(MLEobj, "A", t = t)$A
      val[, t] <- Zt %*% hatxt[, t, drop = FALSE] + At
      if(interval=="confidence") 
        se[, t] <- takediag(Zt %*% tcrossprod( hatVt[, , t], Zt))
      if(interval=="prediction"){
        if (R.time.varying) Rt <- parmat(MLEobj, "R", t = t)$R
        se[, t] <- takediag(Zt %*% tcrossprod( hatVt[, , t], Zt) + Rt)
      }
    }
    
    # Set up output
    if(output=="tibble"){
      data.names <- attr(MLEobj[["model"]], "Y.names")
      data.dims <- attr(MLEobj[["model"]], "model.dims")[["y"]]
      nn <- data.dims[1]
      TT <- data.dims[2]
      ret <- data.frame(
        .rownames = rep(data.names, each = TT),
        t = rep(1:TT, nn),
        y = vec(t(MLEobj[["model"]]$data)),
        stringsAsFactors = FALSE
      )
    }
    
  }
  
  if (type == "x") {
    if (conditioning == "T") hatxt <- MARSSkf(MLEobj)[["xtT"]]
    if (conditioning == "t1") hatxt <- MARSSkfss(MLEobj)[["xtt"]]
    if (interval!="none"){
      if (conditioning == "T") hatVt <- MARSSkf(MLEobj)[["VtT"]]
      if (conditioning == "t1") hatVt <- MARSSkfss(MLEobj)[["Vtt"]]
    }
    
    B.time.varying <- model.dims[["B"]][3] != 1
    U.time.varying <- model.dims[["U"]][3] != 1
    Q.time.varying <- model.dims[["Q"]][3] != 1
    
    val <- matrix(NA, m, TT)
    rownames(val) <- attr(MLEobj[["marss"]], "X.names")
    if(interval!="none") se <- val
    
    x0 <- coef(MLEobj, type = "matrix")[["x0"]]
    if(interval!="none") V0 <- coef(MLEobj, type = "matrix")[["V0"]]
    Bt <- parmat(MLEobj, "B", t = 1)[["B"]]
    Ut <- parmat(MLEobj, "U", t = 1)[["U"]]
    Qt <- parmat(MLEobj, "Q", t = 1)[["Q"]]
    if (MLEobj$model$tinitx == 0){
      val[, 1] <- Bt %*% x0 + Ut
      if (interval=="confidence") se[, 1] <- takediag(Bt %*% tcrossprod(V0, Bt))
      if (interval=="prediction") se[, 1] <- takediag(Bt %*% tcrossprod(V0, Bt) + Qt)
    }
    if (MLEobj[["model"]][["tinitx"]] == 1){
      val[, 1] <- x0
      if(interval != "none") se[, 1] <- takediag(V0)
    }
    for (t in 2:TT) {
      if (B.time.varying) Bt <- parmat(MLEobj, "B", t = t)[["B"]]
      if (U.time.varying) Ut <- parmat(MLEobj, "U", t = t)[["U"]]
      val[, t] <- Bt %*% hatxt[, t - 1, drop = FALSE] + Ut
      if (interval=="confidence") 
        se[, t] <- takediag(Bt %*% tcrossprod( hatVt[, ,t-1], Bt))
      if (interval=="prediction"){
        if (Q.time.varying) Qt <- parmat(MLEobj, "Q", t = t)[["Q"]]
        se[, t] <- takediag(Bt %*% tcrossprod( hatVt[, ,t-1], Bt) + Qt)
      } 
    }
    
    # Set up output
    if(output=="tibble"){
      state.names <- attr(MLEobj[["model"]], "X.names")
      state.dims <- attr(MLEobj[["model"]], "model.dims")[["x"]]
      mm <- state.dims[1]
      TT <- state.dims[2]
      
      if (conditioning == "T")
        ret <- data.frame(
          .rownames = rep(state.names, each = TT),
          t = rep(1:TT, mm),
          xtT = vec(t(hatxt)),
          stringsAsFactors = FALSE
        )
      if (conditioning == "t1")
        ret <- data.frame(
          .rownames = rep(state.names, each = TT),
          t = rep(1:TT, mm),
          xtt = vec(t(hatxt)),
          stringsAsFactors = FALSE
        )
    }
  }
  
  if (interval=="none"){
    if(output=="matrix") return(val)
    retlist = list(.fitted=val)
  }
  if (interval=="confidence"){
    se[se<0 & abs(se)<sqrt(.Machine$double.eps)] <- 0
    se <- sqrt(se) # was not sqrt earlier
    retlist <- list(
      .fitted = val, 
      .se.fit = se,
      .conf.low = val + qnorm(alpha/2) * se,
      .conf.up = val + qnorm(1- alpha/2) * se
    )
  }
  if (interval=="prediction"){
    se[se<0 & abs(se)<sqrt(.Machine$double.eps)] <- 0
    se <- sqrt(se) # was not sqrt earlier
    if (type == "x"){
      retlist <-list(
        .fitted = val, 
        .sd.x = se, 
        .lwr = val + qnorm(alpha/2) * se,
        .upr = val + qnorm(1- alpha/2) * se
      )
    }
    if (type == "y"){
      retlist <- list(
        .fitted=val,
        .sd.y=se,
        .lwr = val + qnorm(alpha/2) * se,
        .upr = val + qnorm(1- alpha/2) * se
      )
    }
  }
  if(output=="matrix") return(retlist)
  return(cbind(ret, as.data.frame(lapply(retlist, function(x){vec(t(x))}))))
} # end of fitted.marssMLE
