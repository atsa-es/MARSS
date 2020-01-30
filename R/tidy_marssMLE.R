###############################################################################################################################################
#  tidy method for class marssMLE
##############################################################################################################################################
tidy.marssMLE <- function(x, type = c("parameters", "xtT", "fitted.ytT", "ytT"),
                          conf.int = TRUE, conf.level = 0.95,
                          form = attr(x[["model"]], "form")[1], ...) {
  ## Argument checking
  type <- match.arg(type)
  conditioning <- "T" # MARSShatyt missing needed var.ytt1 and var.Eytt1 for other cases
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level > 1 || conf.level < 0) stop("tidy.marssMLE: conf.level must be between 0 and 1.", call. = FALSE)
  if (!(conf.int %in% c(TRUE, FALSE))) stop("tidy.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)
  if (type == "xtT") type <- "x"
  if (type == "ytT") type <- "y"
  if (type == "fitted.ytT") type <- "fitted.y"
  if (!(form %in% c("marss", "marxss", "dfa"))) stop("tidy.marssMLE: Allowed forms are marss, marxss, and dfa.\n", call. = FALSE)
  if (length(form) != 1) stop("tidy.marssMLE: Please enter one form from marss, marxss, and dfa.\n", call. = FALSE)
  ## End Argument checking

  alpha <- 1 - conf.level
  csuffix <- switch(conditioning, T = "tT", `t-1` = "tt1", t = "tt")
  extras <- list()

  rerun.MARSSparamCIs <- FALSE
  model.has.cis <- all(c("par.se", "par.lowCI", "par.upCI") %in% names(x))
  if (conf.int & type == "parameters") rerun.MARSSparamCIs <- ifelse(model.has.cis, FALSE, TRUE)
  if (!missing(...)) {
    extras <- list(...)
    if (!all(names(extras) %in% c("rotate", "method", "hessian.fun", "nboot"))) stop("Unknown extra argument. See ?tidy.marssMLE for allowed arguments.\n")
  }

  if (type == "parameters" && conf.int && (!missing(...) || !missing(conf.level))) {
    if (model.has.cis) warning("tidy.marssMLE: Your marssMLE object has par.se and CIs, but you have passed in arguments for calculating CIs.  MARSSparamCIs() will be re-run with these values.\n")
    rerun.MARSSparamCIs <- TRUE
  }

  # set rotate
  rotate <- FALSE
  if ("rotate" %in% names(extras)) {
    if (form != "dfa") stop("tidy.marssMLE: rotate only makes sense if form='dfa'.\n  Pass in form='dfa' if your model is a DFA model, but the form\n attribute is not set (because you set up your DFA model manually). \n", call. = FALSE)
    rotate <- extras[["rotate"]]
    if (!(rotate %in% c(TRUE, FALSE))) stop("tidy.marssMLE: rotate must be TRUE/FALSE. \n", call. = FALSE)
    if( rotate && attr(x[["model"]], "model.dims")[["Z"]][3]!=1 ) stop("tidy.marssMLE: if rotate = TRUE, Z must be time-constant. \n", call. = FALSE)
  }

  if (type == "parameters") {
    ests <- coef(x, type = "vector")
    if (length(ests) == 0) stop("tidy.marssMLE: No estimated parameters in your fitted model.\n", call. = FALSE)
    if (form == "dfa" && rotate && length(x[["par"]][["Z"]]) != 0) {
      stop("tidy.marssMLE: You are requesting the parameters for a DFA \n and requested that the Z matrix be rotated. You need to do the rotation yourself.  See ?tidy.marssMLE for the code.\n", call. = FALSE)
    } else {
      ret <- data.frame(
        term = names(ests),
        estimate = ests
      )
      if (conf.int) {
        if (rerun.MARSSparamCIs) x <- MARSSparamCIs(x, alpha = alpha, ...)
        ret <- cbind(ret,
          std.error = coef(x, type = "vector", what = "par.se"),
          conf.low = coef(x, type = "vector", what = "par.lowCI"),
          conf.high = coef(x, type = "vector", what = "par.upCI")
        )
      }
    }
    rownames(ret) <- NULL
  }
  if (type == "x") {
    xtype <- paste0(type, csuffix)
    model <- x[["model"]]
    state.names <- attr(model, "X.names")
    state.dims <- attr(model, "model.dims")[["x"]]
    mm <- state.dims[1]
    TT <- state.dims[2]
    if(csuffix=="tt"){ 
      kfss <- MARSSkfss(x)
    }else{
      kfss <- MARSSkfas(x)
    }
    states <- kfss[[xtype]]
    vtype <- str_replace(xtype, "x", "V")
    states.se <- apply(kfss[[vtype]], 3, function(x) {
      takediag(x)
    })
    states.se[states.se < 0] <- NA
    states.se <- sqrt(states.se)
    if (mm == 1) states.se <- matrix(states.se, 1, TT)

    # if user specified rotate, 
    # I specified that Z (in marxss form) must be time-constant
    if (form == "dfa" && rotate && length(x[["par"]][["Z"]]) != 0) {
      Z.est <- coef(x, type = "matrix")[["Z"]]
      H <- 1
      if (ncol(Z.est) > 1) {
        H <- solve(varimax(Z.est)[["rotmat"]])
        states <- H %*% states # rotated states
        states.var <- kfss[[vtype]]
        for (t in 1:TT) {
          states.se[, t] <- sqrt(takediag(H %*% states.var[, , t] %*% t(H)))
        }
      }
    }
    ret <- data.frame(
      .rownames = rep(state.names, each = TT),
      t = rep(1:TT, mm),
      estimate = vec(t(states)),
      std.error = vec(t(states.se))
    )
    if (conf.int) {
      conf.low <- qnorm(alpha / 2) * ret$std.error + ret$estimate
      conf.up <- qnorm(1 - alpha / 2) * ret$std.error + ret$estimate
      ret <- cbind(ret,
        conf.low = conf.low,
        conf.high = conf.up
      )
    }
    rownames(ret) <- NULL
  }
  if (type == "y") {
    ytype <- paste0(type, csuffix)
    model <- x[["model"]]
    Y.names <- attr(model, "Y.names")
    Y.dims <- attr(model, "model.dims")[["y"]]
    nn <- Y.dims[1]
    TT <- Y.dims[2]
    hatyt <- MARSShatyt(x, only.kem=FALSE)
    Ey <- hatyt[[ytype]]
    vtype <- str_replace(ytype, "y", "var.y")
    y.var <- hatyt[[vtype]]
    vtype <- str_replace(ytype, "y", "var.Ey")
    Ey.var <- hatyt[[vtype]]
    Ey.se <- apply(Ey.var, 3, function(x) { takediag(x) })
    Ey.se[Ey.se < 0] <- NA
    Ey.se <- sqrt(Ey.se)
    y.sd <- apply(y.var, 3, function(x) { takediag(x) })
    y.sd[y.sd < 0] <- NA
    y.sd <- sqrt(y.sd)
    if (nn == 1) Ey.se <- matrix(Ey.se, 1, TT)
    if (nn == 1) y.sd <- matrix(y.sd, 1, TT)
    
    ret <- data.frame(
      .rownames = rep(Y.names, each = TT),
      t = rep(1:TT, nn),
      y = vec(t(x[["model"]][["data"]])),
      estimate = vec(t(Ey)),
      std.error = vec(t(Ey.se)),
      std.dev = vec(t(y.sd))
    )
    if (conf.int) {
      ret <- cbind(ret,
        conf.low = qnorm(alpha / 2) * ret$std.error + ret$estimate,
        conf.high = qnorm(1 - alpha / 2) * ret$std.error + ret$estimate,
        pred.low = qnorm(alpha / 2) * ret$std.dev + ret$estimate,
        pred.high = qnorm(1 - alpha / 2) * ret$std.dev + ret$estimate
      )
    }
    rownames(ret) <- NULL
  }
  if (type == "fitted.y") {
    Y.names <- attr(x[["model"]], "Y.names")
    Y.dims <- attr(x[["model"]], "model.dims")[["y"]]
    nn <- Y.dims[1]
    TT <- Y.dims[2]
    fit.y.conf <- fitted.marssMLE(x, type="ytT", 
                             interval="confidence", 
                             conf.level=conf.level)
    fit.y.pred <- fitted.marssMLE(x, type="ytT", 
                             interval="prediction", 
                             conf.level=conf.level)

    ret <- fit.y.conf[,c(".rownames","t","y")]
    ret <- cbind( ret,
      estimate = fit.y.conf$.fitted,
      std.error = fit.y.conf$.se.fit,
      std.dev = fit.y.pred$.sd.y
    )
    if (conf.int) {
      ret <- cbind(ret,
                   conf.low = qnorm(alpha / 2) * ret$std.error + ret$estimate,
                   conf.high = qnorm(1 - alpha / 2) * ret$std.error + ret$estimate,
                   pred.low = qnorm(alpha / 2) * ret$std.dev + ret$estimate,
                   pred.high = qnorm(1 - alpha / 2) * ret$std.dev + ret$estimate
      )
    }
    rownames(ret) <- NULL
  }
  ret
}
