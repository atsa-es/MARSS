###############################################################################################################################################
#  fitted method for class marssMLE.
#  returns the fitted value of y conditioned on all the data or data up to t-1 if one.step.ahead=TRUE
##############################################################################################################################################
fitted.marssMLE <- function(object, ...,
                            type = c("observations", "states", "y", "x"),
                            conditioning = c("T", "t", "t-1")) {
  type <- match.arg(type)
  if (type == "y") type <- "observations"
  if (type == "x") type <- "states"
  conditioning <- match.arg(conditioning)
  MLEobj <- object
  if (is.null(MLEobj[["par"]])) {
    stop("fitted.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.")
  }
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
    Z.time.varying <- model.dims[["Z"]][3] != 1
    A.time.varying <- model.dims[["A"]][3] != 1

    val <- matrix(NA, n, TT)
    rownames(val) <- attr(MLEobj$marss, "Y.names")

    Zt <- parmat(MLEobj, "Z", t = 1)$Z
    At <- parmat(MLEobj, "A", t = 1)$A

    for (t in 1:TT) {
      # parmat returns marss form
      if (Z.time.varying) Zt <- parmat(MLEobj, "Z", t = t)$Z
      if (A.time.varying) At <- parmat(MLEobj, "A", t = t)$A
      val[, t] <- Zt %*% hatxt[, t, drop = FALSE] + At
    }
  }

  if (type == "states") {
    if (conditioning == "T") hatxt <- MLEobj[["states"]]
    if (conditioning == "t-1") hatxt <- MARSSkf(MLEobj)[["xtt"]]

    B.time.varying <- model.dims[["B"]][3] != 1
    U.time.varying <- model.dims[["U"]][3] != 1

    val <- matrix(NA, m, TT)
    rownames(val) <- attr(MLEobj$marss, "X.names")

    x0 <- coef(MLEobj, type = "matrix")[["x0"]]
    Bt <- parmat(MLEobj, "B", t = 1)[["B"]]
    Ut <- parmat(MLEobj, "U", t = 1)[["U"]]
    if (MLEobj$model$tinitx == 0) val[, 1] <- Bt %*% x0 + Ut
    if (MLEobj$model$tinitx == 1) val[, 1] <- x0
    for (t in 2:TT) {
      if (B.time.varying) Bt <- parmat(MLEobj, "B", t = t)[["B"]]
      if (U.time.varying) Ut <- parmat(MLEobj, "U", t = t)[["U"]]
      val[, t] <- Bt %*% hatxt[, t - 1, drop = FALSE] + Ut
    }
  }
  val
} # end of fitted.marssMLE
