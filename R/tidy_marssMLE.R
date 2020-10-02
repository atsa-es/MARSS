###############################################################################################################################################
#  tidy method for class marssMLE
##############################################################################################################################################
tidy.marssMLE <-
  function(x,
           conf.int = TRUE,
           conf.level = 0.95,
           ...) {
    ## Argument checking
    if (conf.int &&
      (!is.numeric(conf.level) ||
        length(conf.level) != 1 ||
        conf.level > 1 ||
        conf.level < 0)) {
      stop("tidy.marssMLE: conf.level must be between 0 and 1.", call. = FALSE)
    }
    if (!(conf.int %in% c(TRUE, FALSE))) {
      stop("tidy.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)
    }
    ## End Argument checking

    alpha <- 1 - conf.level

    extras <- list()

    rerun.MARSSparamCIs <- FALSE
    model.has.cis <-
      all(c("par.se", "par.lowCI", "par.upCI") %in% names(x))
    if (conf.int) {
      rerun.MARSSparamCIs <- ifelse(model.has.cis, FALSE, TRUE)
    }
    if (!missing(...)) {
      extras <- list(...)
      if (!all(names(extras) %in% c("method", "hessian.fun", "nboot"))) {
        stop("Unknown extra argument. See ?tidy.marssMLE for allowed arguments.\n")
      }
    }

    if (conf.int && (!missing(...) || !missing(conf.level))) {
      if (model.has.cis) {
        warning(
          "tidy.marssMLE: Your marssMLE object has par.se and CIs, but you have passed in arguments for calculating CIs.  MARSSparamCIs() will be re-run with these values.\n"
        )
      }
      rerun.MARSSparamCIs <- TRUE
    }

    ests <- coef(x, type = "vector")

    ret <- data.frame(
      term = names(ests),
      estimate = ests,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      if (rerun.MARSSparamCIs) {
        x <- MARSSparamCIs(x, alpha = alpha, ...)
      }
      ret <- cbind(
        ret,
        std.error = coef(x, type = "vector", what = "par.se"),
        conf.low = coef(x, type = "vector", what = "par.lowCI"),
        conf.up = coef(x, type = "vector", what = "par.upCI")
      )
    }

    rownames(ret) <- NULL

    ret
  }
