plot.marssMLE <-
  function(x,
           plot.type = c(
             "model.ytT", "xtT",
             "model.resids.ytt1", "qqplot.model.resids.ytt1", "acf.model.resids.ytt1",
             "std.model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
             "model.resids.ytT", "qqplot.model.resids.ytT",
             "std.model.resids.ytT", "qqplot.std.model.resids.ytT",
             "model.resids.ytt", "qqplot.model.resids.ytt",
             "std.model.resids.ytt", "qqplot.std.model.resids.ytt",
             "state.resids.xtT", "qqplot.state.resids.xtT",
             "std.state.resids.xtT", "qqplot.std.state.resids.xtT",
             "ytT", "residuals", "all"
           ),
           form = c("marxss", "marss", "dfa"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE, pi.int = FALSE,
           plot.par = list(),
           ..., silent = FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for autoplot.marssMLE. Please install it.", call. = FALSE)
    }
    if (!inherits(x, "marssMLE")) {
      stop("autoplot.marssMLE: x must be class marssMLE.", call. = FALSE)
    }

    # Argument checks: plot.type
    if (missing(plot.type)) {
      plot.type <- c(
        "model.ytT", "xtT",
        "model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
        "std.model.resids.ytT",
        "std.state.resids.xtT", "qqplot.std.state.resids.xtT"
      )
    }
    plot.type <- match.arg(plot.type, several.ok = TRUE)
    if (identical(plot.type, "residuals")) {
      plot.type <- c(
        "model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
        "std.model.resids.ytT", "std.state.resids.xtT", "qqplot.std.state.resids.xtT"
      )
    }
    if (identical(plot.type, "all")) {
      plot.type <- eval(formals()$plot.type)
      plot.type <- plot.type[!(plot.type %in% c("residuals", "all"))]
    }

    if (!is.numeric(conf.level) || length(conf.level) > 1 || conf.level > 1 || conf.level < 0) stop("autoplot.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
    if (!(conf.int %in% c(TRUE, FALSE))) stop("autoplot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)

    if (missing(form)) {
      model_form <- attr(x[["model"]], "form")[1]
    } else {
      model_form <- match.arg(form)
    }

    # Argument checks: plotpar
    plotpar <- list(
      point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
      line.col = "black", line.size = 1, line.linetype = "solid",
      ci.fill = "grey70", ci.col = "grey70", ci.linetype = "solid",
      ci.linesize = 0, ci.alpha = 0.6
    )
    if (!is.list(plot.par)) stop("autoplot.marssMLE: plot.par must be a list.", call. = FALSE)
    if (!missing(plot.par)) {
      if (!all(names(plot.par) %in% names(plotpar))) {
        stop(paste0("autoplot.marssMLE: Allowed plot.par names are ", paste(names(plotpar), collapse = ", "), ".\n"), call. = FALSE)
      } else {
        for (i in names(plot.par)) {
          plotpar[[i]] <- plot.par[[i]]
        }
      }
    }

    extras <- list()
    if (!missing(...)) {
      extras <- list(...)
      allowednames <- c("rotate", "method", "hessian.fun", "nboot")
      bad.names <- names(extras)[!(names(extras) %in% allowednames)]
      if (!all(names(extras) %in% allowednames)) stop(paste("autoplot.marssMLE:", paste(bad.names, collapse = " "), "is/are unknown argument(s). See ?autoplot.marssMLE for allowed arguments.\n"), call. = FALSE)
      if (model_form != "dfa" & "rotate" %in% names(extras)) {
        cat("autoplot.marssMLE: 'rotate' argument is ignored if form!='dfa'\n Pass in form='dfa' if your model is a DFA model, but the form \n attribute is not set (because you set up your DFA model manually).\n\n")
        rotate <- FALSE
      }
    }
    # End Argument checks

    alpha <- 1 - conf.level

    # If user requests any residuals plots, set up the residuals data frames
    resids <- c()
    if (any(stringr::str_detect(plot.type, "tt1"))) {
      resids <- residuals.marssMLE(x, type = "tt1", standardization = "Cholesky")
    }
    if (any(stringr::str_detect(plot.type, "tt") & !stringr::str_detect(plot.type, "tt1"))) {
      resids <- rbind(resids, residuals.marssMLE(x, type = "tt", standardization = "Cholesky"))
    }
    if (any(stringr::str_detect(plot.type, "tT"))) {
      resids <- rbind(resids, residuals.marssMLE(x, type = "tT", standardization = "Cholesky"))
    }

    if ("model.ytT" %in% plot.type) {
      # make plot of observations
      df <- fitted.marssMLE(x, type = "ytT", interval = "confidence", level = conf.level)
      df$ymin <- df$.conf.low
      df$ymax <- df$.conf.up
      if (pi.int) {
        df2 <- fitted.marssMLE(x, type = "ytT", interval = "prediction", level = conf.level)
        df$ymin.pi <- df2$.lwr
        df$ymax.pi <- df2$.upr
      }
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        if (conf.int) tit <- paste(tit, "+ CI")
        if (pi.int) tit <- paste(tit, "+ PI (dashed)")
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.fitted, y, ymin, ymax, na.rm = TRUE), max(.fitted, y, ymin, ymax, na.rm = TRUE))
          plot(t, .fitted, type = "l", xlab = "", ylab = "Estimate", ylim = ylims)
          title(tit)
          if (conf.int) polygon(c(t, rev(t)), c(ymin, rev(ymax)), col = plotpar$ci.col, border = plotpar$ci.border)
          if (decorate) points(t, y, col = plotpar$point.col, pch = plotpar$point.pch, cex = plotpar$point.size)
          lines(t, .fitted, col = plotpar$line.col, lwd = plotpar$line.lwd)
          if (pi.int) {
            lines(t, ymin.pi, col = "black", lwd = 1, lty = 2)
            lines(t, ymax.pi, col = "black", lwd = 1, lty = 2)
          }
          box()
        })
      }
      plot.type <- plot.type[plot.type != "model.ytT"]
      if (!silent) cat(paste("plot type = \"model.ytT\" Observations with fitted values\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("xtT" %in% plot.type) {
      # make plot of states and CIs

      if ("rotate" %in% names(extras)) {
        rotate <- extras[["rotate"]]
        if (!(rotate %in% c(TRUE, FALSE))) stop("plot.marssMLE: rotate must be TRUE/FALSE. \n")
      } else {
        rotate <- FALSE
      }

      states <- tsSmooth.marssMLE(x, type = "xtT", interval = ifelse(conf.int, "confidence", "none"), level = conf.level, ...)
      if (model_form == "dfa") {
        if (rotate) {
          rottext <- "rotated"
        } else {
          rottext <- ""
        }
        states$.rownames <- paste("DFA", rottext, "trend", states$.rownames)
      } else {
        states$.rownames <- paste0("State ", states$.rownames)
      }

      nX <- min(9, attr(x$model, "model.dims")$x[1])
      plot.nrow <- round(sqrt(nX))
      plot.ncol <- ceiling(nX / plot.nrow)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(states$.rownames)) {
        with(subset(states, states$.rownames == plt), {
          ylims <- c(min(.estimate, .conf.low, na.rm = TRUE), max(.estimate, .conf.up, na.rm = TRUE))
          plot(t, .estimate, type = "l", xlab = "", ylab = "Estimate", ylim = ylims)
          title(plt)
          if (conf.int) polygon(c(t, rev(t)), c(.conf.low, rev(.conf.up)), col = plotpar$ci.col, border = plotpar$ci.border)
          lines(t, .estimate)
          box()
        })
      }
      plot.type <- plot.type[plot.type != "xtT"]
      if (!silent) cat(paste("plot type = \"xtT\" Estimated states\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    resids.vs.time.plots <- paste0(rep(c("model.resids", "std.model.resids"), each = 3), c(".ytT", ".ytt", ".ytt1"))
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1] # ytT, ytt or ytt1
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
      df$.resids[is.na(df$value)] <- NA
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.resids, na.rm = TRUE), max(.resids, na.rm = TRUE))
          plot(t, .resids,
            type = "p", xlab = "",
            ylab = "", ylim = ylims,
            col = plotpar$point.col, pch = plotpar$point.pch,
            cex = plotpar$point.size
          )
          if (decorate) {
            lo <- predict(loess(.resids ~ t), newdata = data.frame(t = t), se = TRUE)
            lo.t <- names(lo$fit)
            if (stringr::str_detect(i, "std")) {
              sigma <- 1
            } else {
              sigma <- .sigma
              sigma[is.na(value)] <- 0
            }
            ymin <- qnorm(alpha / 2) * sigma
            ymax <- -qnorm(alpha / 2) * sigma
            ylims <- c(min(.resids, ymin, na.rm = TRUE), max(.resids, ymax, na.rm = TRUE))
            lines(t, lo$fit, col = plotpar$line.col, lwd = plotpar$line.lwd)
          }
          title(plt)
          box()
          abline(h = 0, lty = 3)
        })
        mtext(ifelse(stringr::str_detect(i, "std"), "Standardized observation residuals, y - E[y]", "Observation residuals, y - E[y]"), side = 2, outer = TRUE, line = -1)
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    resids.vs.time.plots <- c("state.resids.xtT", "std.state.resids.xtT")
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      # make plot of process residuals; set form='marxss' to get process resids
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1]
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
      df$.rownames <- paste0("State ", df$.rownames)
      nX <- min(9, attr(x$model, "model.dims")$x[1])
      plot.nrow <- round(sqrt(nX))
      plot.ncol <- ceiling(nX / plot.nrow)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.resids, na.rm = TRUE), max(.resids, na.rm = TRUE))
          plot(t, .resids,
            type = "p", xlab = "",
            ylab = "", ylim = ylims,
            col = plotpar$point.col, pch = plotpar$point.pch,
            cex = plotpar$point.size
          )
          if (decorate) {
            lo <- predict(loess(.resids ~ t), newdata = data.frame(t = t), se = TRUE)
            lo.t <- names(lo$fit)
            sigma <- .sigma
            ymin <- qnorm(alpha / 2) * sigma
            ymax <- -qnorm(alpha / 2) * sigma
            ylims <- c(min(.resids, ymin, na.rm = TRUE), max(.resids, ymax, na.rm = TRUE))
            lines(t, lo$fit, col = plotpar$line.col, lwd = plotpar$line.lwd)
          }

          title(plt)
          box()
          abline(h = 0, lty = 3)
        })
        mtext(ifelse(stringr::str_detect(i, "std"), "Standardized state residuals, xtT - E[x]", "State residuals, xtT - E[x]"),
          side = 2, outer = TRUE, line = -1
        )
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }


    # QQplots for normality
    slp <- function(yy) {
      y <- quantile(yy[!is.na(yy)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y) / diff(x)
      return(slope)
    }
    int <- function(yy) {
      y <- quantile(yy[!is.na(yy)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y) / diff(x)
      int <- y[1L] - slope * x[1L]
      return(int)
    }

    qqplot.plots <- plot.type[stringr::str_detect(plot.type, "qqplot")]
    for (i in qqplot.plots[qqplot.plots %in% plot.type]) {
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
      slope <- tapply(df$.resids, df$.rownames, slp)
      intercept <- tapply(df$.resids, df$.rownames, int)
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          qqnorm(.resids, main = plt)
          abline(a = intercept[plt], b = slope[plt], col = plotpar$line.col, lwd = plotpar$line.lwd)
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("ytT" %in% plot.type) {
      # make plot of expected value of y
      df <- tsSmooth.marssMLE(x, type = "ytT", interval = ifelse(conf.int, "confidence", "none"))
      if (conf.int) {
        df$ymin <- df$.conf.low
        df$ymax <- df$.conf.up
      }
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        if (conf.int) tit <- paste(tit, "+ CI")
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.estimate, y, ymin, ymax, na.rm = TRUE), max(.estimate, y, ymin, ymax, na.rm = TRUE))
          plot(t, .estimate, type = "l", xlab = "", ylab = "Estimate", ylim = ylims)
          title(tit)
          if (conf.int) polygon(c(t, rev(t)), c(ymin, rev(ymax)), col = plotpar$ci.col, border = plotpar$ci.border)
          points(t, y, col = plotpar$point.col, pch = plotpar$point.pch, cex = plotpar$point.size)
          lines(t, .estimate, col = plotpar$line.col, lwd = plotpar$line.lwd)
          box()
        })
      }
      plot.type <- plot.type[plot.type != "ytT"]
      if (!silent) cat(paste("plot type = \"ytT\" Expected value of Y conditioned on data\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    acf.plots <- plot.type[stringr::str_detect(plot.type, "acf")]
    for (i in acf.plots[acf.plots %in% plot.type]) {
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        with(subset(df, df$.rownames == plt), {
          stats::acf(.resids, na.action = na.pass, main = "")
          title(tit)
          box()
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }
  }
