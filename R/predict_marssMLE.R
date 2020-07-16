######################################################################################################  predict method for class marssMLE. Prediction intervals
##################################################################################
predict.marssMLE <- function(object, n.ahead = 0,
                             conf.level = c(0.80, 0.95),
                             type = c("ytT", "xtT"),
                             newdata = list(t = NULL, y = NULL, c = NULL, d = NULL),
                             interval = c("prediction", "confidence", "none"),
                             fun.kf = c("MARSSkfas", "MARSSkfss"),
                             x0 = "reestimate",
                             ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  fun.kf <- match.arg(fun.kf)
  if (object[["fun.kf"]] != fun.kf) message(paste0(fun.kf, "is being used for prediction. This is different than fun.kf in the marssMLE object."))
  MODELobj <- object[["model"]]
  model.dims <- attr(MODELobj, "model.dims")
  TT <- model.dims[["y"]][2]
  nx <- switch(type,
    ytT = model.dims[["y"]][1],
    xtT = model.dims[["x"]][1]
  )

  if (is.null(object[["par"]])) {
    stop("predict.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (interval == "none") conf.level <- c()
  if (length(conf.level) > 0 && (!is.numeric(conf.level) || any(conf.level > 1) || any(conf.level < 0))) {
    stop("predict.marssMLE: conf.level must be between 0 and 1.", call. = FALSE)
  }
  if (length(conf.level) == 0) interval <- "none"
  if (!missing(n.ahead) && (!is.numeric(n.ahead) || length(n.ahead) != 1 
                            || n.ahead < 0 || (n.ahead %% 1) != 0)) {
    stop("predict.marssMLE: n.ahead must be an integer >= 0.", call. = FALSE)
  }
  missingx0 <- TRUE
  is.matrixx0 <- FALSE
  if (!missing(x0)) {
    missingx0 <- FALSE
    if (inherits(x0, "matrix")) is.matrixx0 <- TRUE
  }

  # test x0 structure
  if (!(identical(x0, "reestimate") |
    identical(x0, "use.model") |
    inherits(x0, "matrix"))) {
    stop("predict.marssMLE: x0 must be 'reestimate', 'use.model' or a new x0 matrix to use as the x0 value.", call. = FALSE)
  }
  dimx0 <- attr(object$model, "model.dims")$x0[1:2]
  if (inherits(x0, "matrix") && !identical(dim(x0), dimx0)) {
    stop(paste0("predict.marssMLE: if x0 is a matrix, it must be the same dimensions as the model x0 (", dimx0[1], " x ", dimx0[2], ")."), call. = FALSE)
  }
  # Use the x0 estimated in the fitted model (object)
  if (identical(x0, "use.model")) x0list <- list(x0 = coef(object, type = "matrix")[["x0"]], tinitx = object[["model"]][["tinitx"]])
  # Set up x0 to be the same form (sharing, independent, etc) as in the model call
  if (identical(x0, "reestimate")) x0list <- list(x0 = object[["call"]][["model"]][["x0"]], tinitx = object[["call"]][["model"]][["tinitx"]])
  if (inherits(x0, "matrix")) x0list <- list(x0 = x0, tinitx = object[["model"]][["tinitx"]])
  extras <- list()
  if (!missing(...)) {
    extras <- list(...)
  }

  if (!missing(n.ahead) && n.ahead > 0) {
    if (!missingx0) message("predict.marssMLE: x0 was passed in. This is ignored since n.ahead > 0 (forecast). x at time T is used for x0.")

    outlist <- forecast.marssMLE(object,
      h = n.ahead, conf.level = conf.level,
      interval = interval,
      type = type, newdata = newdata,
      fun.kf = fun.kf, ...
    )
    
    return(outlist)
  }

  # h not passed and/or h=0

  # check newdata for errors
  if (!is.null(newdata[["t"]])) {
    if (!is.numeric(newdata[["t"]]) || !is.vector(newdata[["t"]])) {
      stop("predict.marssMLE: t in newdata must be a numeric vector.", call. = FALSE)
    }
    if (!all(diff(newdata[["t"]]) == 1) || newdata[["t"]][1] < 1) {
      stop("predict.marssMLE: t in newdata must be a positive ordered sequence of integers, one time step apart (like 1,2,3...).", call. = FALSE)
    }
  }
  nonewdata <- all(unlist(lapply(newdata[c("y", "c", "d")], is.null)))
  isnullnewdata <- all(unlist(lapply(newdata[c("t", "y", "c", "d")], is.null)))
  if (nonewdata && !is.null(newdata[["t"]])) {
    if (length(newdata[["t"]]) > TT || newdata[["t"]][length(newdata[["t"]])] > TT) {
      stop("predict.marssMLE: if no y, c or d in newdata, t must be within the original time steps.", call. = FALSE)
    } else {
      newdata[["y"]] <- object[["model"]][["data"]][, newdata[["t"]], drop = FALSE]
      message(paste0("predict.marssMLE(): prediction is conditioned on the model data t = ", newdata[["t"]][1], " to ", max(newdata[["t"]]), "."))
    }
  }
  if (isnullnewdata) newdata[["t"]] <- 1:TT

  # h=0 if here
  if (isnullnewdata) { # Use original data if no new data
    if (is.matrixx0) { # Need to reestimate states and ytT
      new.MODELlist <- coef(object, type = "matrix", form = "marxss")
      new.MODELlist$c <- object$call$model$c
      new.MODELlist$d <- object$call$model$d
      new.MODELlist$x0 <- x0
      newMLEobj <- MARSS(object$call$data,
        model = new.MODELlist, silent = TRUE,
        method = object[["method"]],
        fun.kf = object$call$fun.kf, form = "marxss"
      )
    } else {
      newMLEobj <- object
    }
  } else {
    # newdata was passed in. Need to make newMLEobj
    if (!is.list(newdata) || !all(names(newdata) %in% c("t", "y", "c", "d"))) {
      stop("predict.marssMLE: newdata must be list with only t, y, c and/or d.", call. = FALSE)
    }

    # We need the model in marxss form
    new.MODELlist <- coef(object, type = "matrix", form = "marxss")

    # Does the model contain covariates?
    isxreg <- list(c = TRUE, d = TRUE, y = TRUE)
    for (elem in c("c", "d")) {
      tmp <- coef(object, type = "matrix", form = "marxss")[[elem]]
      if (dim(tmp)[1] == 1 && dim(tmp)[2] == 1 && tmp == 0) {
        isxreg[[elem]] <- FALSE
      } else {
        dim(new.MODELlist[[elem]]) <- model.dims[[elem]][c(1, 3)]
      }
    }

    # if y is null, then c or d was passed in
    if (identical(newdata[["y"]], "model")) {
      newdata[["y"]] <- object[["model"]][["data"]]
      tmp <- unlist(lapply(newdata[c("c", "d")], ncol))
      ncol.cd <- unname(tmp[1])
      if (is.null(newdata[["t"]]) && ncol.cd <= TT) newdata[["t"]] <- 1:ncol.cd
      if (!is.null(newdata[["t"]])) {
        if (max(newdata[["t"]]) > TT) stop("predict.marssMLE: newdata t must be within the original data if using model data for y.", call. = FALSE)
        newdata[["y"]] <- newdata[["y"]][, newdata[["t"]], drop = FALSE]
        message(paste("predict.marssMLE(): prediction is conditioned on the model data t =", newdata[["t"]][1], "to", max(newdata[["t"]]), "."))
      }
    }
    if (identical(newdata[["y"]], "none") || is.null(newdata[["y"]])) {
      tmp <- unlist(lapply(newdata[c("c", "d")], ncol))
      ncol.cd <- unname(tmp[1])
      newdata[["y"]] <- matrix(as.numeric(NA), model.dims[["y"]][1], ncol.cd)
      message("predict.marssMLE(): prediction is not conditioned on any data, only c or d covariates.")
    }
    for (elem in c("y", "c", "d")) {
      if (!isxreg[[elem]]) {
        if (!is.null(newdata[[elem]])) {
          message(paste0("predict.marssMLE(): model does not include ", elem, ". ", elem, " in newdata is being ignored."))
        }
      } else {
        if (is.null(newdata[[elem]])) {
          stop(paste0("predict.marssMLE(): model includes ", elem, ". ", elem, " must be in newdata."), call. = FALSE)
        }
        if (!is.numeric(newdata[[elem]])) {
          stop("predict.marssMLE: y, c, and d in newdata must be numeric (use class() and is.numeric() to test what you are passing in).", call. = FALSE)
        }
        if (is.vector(newdata[[elem]])) newdata[[elem]] <- matrix(newdata[[elem]], nrow = 1)
        if (!is.matrix(newdata[[elem]])) stop(paste0("predict.marssMLE(): newdata ", elem, " must be a matrix with ", model.dims[[elem]][1], " rows."), call. = FALSE)
        if (dim(newdata[[elem]])[1] != model.dims[[elem]][1]) stop(paste0("predict.marssMLE(): model ", elem, " has ", model.dims[[elem]][1], " rows.", elem, " in newdata does not."), call. = FALSE)

        if (elem != "y") new.MODELlist[[elem]] <- newdata[[elem]]
        if (elem == "y") new.data <- newdata[[elem]]
      }
    }
    # check that ncols match
    tmp <- unlist(lapply(newdata[c("y", "c", "d")], ncol))
    ncol.newdata <- unname(tmp[1])
    if (!all(tmp == ncol.newdata)) stop("predict.marssMLE(): y, c, and d in newdata must all have the same number of columns.", call. = FALSE)
    # check that t length matches ncols
    if (!is.null(newdata[["t"]]) && length(newdata[["t"]]) != ncol.newdata) {
      stop("predict.marssMLE(): t in newdata must be the same length as the number of columns in y, c and d.", call. = FALSE)
    }

    # use x0 in model if no data and user didn't specify not to use x0
    nodata <- (all(is.na(newdata[["y"]])) || identical(newdata[["y"]], "none"))
    if (missingx0 && nodata) {
      x0list <- list(x0 = coef(object, type = "matrix")[["x0"]], tinitx = object[["model"]][["tinitx"]])
      message("predict.marssMLE(): x0 and tinitx from model are being used for prediction.")
    }
    if (!missingx0 && identical(x0, "reestimate") && nodata) {
      stop("predict.marssMLE(): to reestimate x0 (x0='reestimate'), data (y in newdata) are required.", call. = FALSE)
    }

    # check if parameters are time-varying. If so, t used to specify which parameters to use
    for (elem in names(new.MODELlist)) {
      if (elem %in% c("c", "d", "x0", "V0")) next
      if (model.dims[[elem]][3] == TT) {
        if (is.null(newdata[["t"]])) stop("predict.marssMLE(): the model has time-varying parameters. t in newdata needed in order to specify which parameters to use.", call. = FALSE)
        nend <- newdata[["t"]][length(newdata[["t"]])]
        nstart <- newdata[["t"]][1]
        tmp <- array(NA, dim = c(model.dims[[elem]][1:2], ncol.newdata))
        if (nstart <= TT && nend > TT) {
          tmp[, , 1:(TT - nstart + 1)] <- new.MODELlist[[elem]][, , nstart:TT, drop = FALSE]
          tmp[, , (TT - nstart + 2):ncol.newdata] <- new.MODELlist[[elem]][, , TT, drop = FALSE]
          message(paste0(elem, " is time-varying. The value at t = ", TT, " is used for any newdata t past the original data."))
        }
        if (nstart > TT) {
          tmp[, , 1:ncol.newdata] <- new.MODELlist[[elem]][, , TT, drop = FALSE]
          message(paste0(elem, " is time-varying. The value at ", TT, " is used since newdata t past the original data."))
        }
        if (nend <= TT) {
          tmp[, , 1:ncol.newdata] <- new.MODELlist[[elem]][, , newdata[["t"]], drop = FALSE]
        }
        new.MODELlist[[elem]] <- tmp
      }
    }
    # Passed all checks. Can set t now if still null
    if (is.null(newdata[["t"]])) newdata[["t"]] <- 1:ncol.newdata

    # The x0 values based on use.initial.values is set at top
    new.MODELlist[["tinitx"]] <- x0list[["tinitx"]]
    new.MODELlist[["x0"]] <- x0list[["x0"]]
    newMLEobj <- MARSS(newdata[["y"]],
      model = new.MODELlist, silent = TRUE,
      method = object[["method"]], control = object$call$control,
      inits = object$call$inits, 
      fun.kf = object$call$fun.kf, form = "marxss"
    )
  } # end setting up newMLEobj

  if (type == "ytT") {
    if (length(conf.level) == 0) interval <- "none"
    cols <- switch(interval,
      prediction = c(".rownames", "t", "y", ".fitted", ".sd.y", ".lwr", ".upr"),
      none = c(".rownames", "t", "y", ".fitted"),
      confidence = c(".rownames", "t", "y", ".fitted", ".se.fit", ".conf.low", ".conf.up")
    )
    ret <- fitted.marssMLE(newMLEobj, type = type, interval = interval, conf.level = conf.level[1])[cols]
    colnames(ret)[which(colnames(ret) == ".fitted")] <- "estimate"
    colnames(ret)[which(colnames(ret) == ".sd.y")] <- "se"
    colnames(ret)[which(colnames(ret) == ".se.fit")] <- "se"
    colnames(ret)[which(colnames(ret) == ".lwr")] <- paste("Lo", 100 * conf.level[1])
    colnames(ret)[which(colnames(ret) == ".upr")] <- paste("Hi", 100 * conf.level[1])
    colnames(ret)[which(colnames(ret) == ".conf.low")] <- paste("Lo", 100 * conf.level[1])
    colnames(ret)[which(colnames(ret) == ".conf.up")] <- paste("Hi", 100 * conf.level[1])
    if (interval != "none" && length(conf.level) > 1) {
      for (i in 2:length(conf.level)) {
        cols <- switch(interval,
          prediction = c(".lwr", ".upr"),
          confidence = c(".conf.low", ".conf.up")
        )
        tmp <- fitted.marssMLE(newMLEobj, type = type, interval = interval, conf.level = conf.level[i])[cols]
        colnames(tmp) <- paste(c("Lo", "Hi"), 100 * conf.level[i])
        ret <- cbind(ret, tmp)
      }
    }
  }
  if (type == "xtT") {
    conf.int <- TRUE
    if (length(conf.level) == 0 || interval == "none") conf.int <- FALSE
    ret <- tidy.marssMLE(newMLEobj, type = type, conf.int = conf.int, conf.level = conf.level[1])
    colnames(ret)[which(colnames(ret) == "std.error")] <- "se"
    colnames(ret)[which(colnames(ret) == "conf.low")] <- paste("Lo", 100 * conf.level[1])
    colnames(ret)[which(colnames(ret) == "conf.high")] <- paste("Hi", 100 * conf.level[1])
    if (conf.int && length(conf.level) > 1) {
      for (i in 2:length(conf.level)) {
        tmp <- tidy.marssMLE(newMLEobj, type = type, conf.int = conf.int, conf.level = conf.level[i])[c("conf.low", "conf.high")]
        colnames(tmp) <- paste(c("Lo", "Hi"), 100 * conf.level[i])
        ret <- cbind(ret, tmp)
      }
    }
  }

  # set t in ret with t in newdata
  ret$t <- rep(newdata[["t"]], nx)

  # Set up output
  outlist <- list(
    method = c("MARSS", object[["method"]]),
    model = object,
    level = 100 * conf.level,
    type = type,
    pred = ret,
    t = newdata[["t"]],
    h = n.ahead,
    n.ahead = n.ahead,
    x0 = x0list[["x0"]],
    tinitx = x0list[["tinitx"]],
    newdata = newdata
  )

  class(outlist) <- "marssPredict"

  return(outlist)
}
