plot.marssMLE <-
  function(x,
           plot.type = c("model.ytT", "xtT", "model.resids", "state.resids", "qqplot.model.resids", "qqplot.state.resids", "ytT", "acf.model.resids", "acf.state.resids"),
           form = c("marxss", "marss", "dfa"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE, pi.int = FALSE,
           plot.par = list(), ...) {

    # Argument checks
    plot.type <- match.arg(plot.type, several.ok = TRUE)
    old.plot.type <- c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqplot", "state.residuals.qqplot", "expected.value.observations")
    new.plot.type <- c("model.ytT", "xtT", "model.resids", "state.resids", "qqplot.model.resids", "qqplot.state.resids", "ytT")
    for (i in 1:NROW(old.plot.type)) if (old.plot.type[i] %in% plot.type) plot.type[plot.type == old.plot.type[i]] <- new.plot.type[i]
    if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level > 1 || conf.level < 0) stop("plot.marssMLE: conf.level must be between 0 and 1.", call. = FALSE)
    if (!(conf.int %in% c(TRUE, FALSE))) stop("plot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)

    # run residuals if needed
    if (any(str_detect(plot.type, "acf"))) {
      tt1.resids <- residuals.marssMLE(x, type = "tt1", standardization = "Cholesky")
    }
    if (any(str_detect(plot.type, "resids"))) {
      tT.resids <- residuals.marssMLE(x, type = "tT", standardization = "Cholesky")
    }

    if (missing(form)) {
      model_form <- attr(x[["model"]], "form")[1]
    } else {
      model_form <- match.arg(form)
    }

    plotpar <- list(
      point.pch = 20, point.col = "blue", point.fill = "blue", point.size = 1,
      line.col = "black", line.size = 1, line.linetype = "solid",
      ci.fill = "grey70", ci.col = "grey70", ci.border = FALSE,
      ci.linesize = 0, ci.alpha = 0.6
    )
    if (!is.list(plot.par)) stop("plot.marssMLE: plot.par must be a list.", call. = FALSE)
    if (!missing(plot.par)) {
      if (!all(names(plot.par) %in% names(plotpar))) {
        stop(paste0("plot.marssMLE: Allowed plot.par names are ", paste(names(plotpar), collapse = ", "), ".\n"), call. = FALSE)
      } else {
        for (i in names(plot.par)) {
          plotpar[[i]] <- plot.par[[i]]
        }
      }
    }

    extras <- list()
    if (!missing(...)) {
      extras <- list(...)
      allowednames <- c("rotate")
      bad.names <- names(extras)[!(names(extras) %in% allowednames)]
      if (!all(names(extras) %in% allowednames)) stop(paste("plot.marssMLE:", paste(bad.names, collapse = " "), "is/are unknown argument(s). See ?fitted.marssMLE for allowed arguments.\n"), call. = FALSE)
      if (model_form != "dfa" & "rotate" %in% names(extras)) {
        cat("plot.marssMLE: 'rotate' argument is ignored if form!='dfa'\n Pass in form='dfa' if your model is a DFA model, but the form \n attribute is not set (because you set up your DFA model manually).\n\n")
        rotate <- FALSE
      }
    }
    # End Argument checks

    alpha <- 1 - conf.level

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
      cat(paste("plot type = \"model.ytT\" Observations with fitted values\n"))
      if (length(plot.type) != 0) {
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
      cat(paste("plot type = \"xtT\" Estimated states\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("model.resids" %in% plot.type) {
      # make plot of observation residuals
      df <- subset(tT.resids, tT.resids$name == "model")
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
            sigma <- .sigma
            sigma[is.na(value)] <- 0
            ymin <- qnorm(alpha / 2) * sigma
            ymax <- -qnorm(alpha / 2) * sigma
            ylims <- c(min(.resids, ymin, na.rm = TRUE), max(.resids, ymax, na.rm = TRUE))
            lines(t, lo$fit, col = plotpar$line.col, lwd = plotpar$line.lwd)
          }
          title(plt)
          box()
          abline(h = 0, lty = 3)
        })
        mtext("Observation residuals, y - E[y]", side = 2, outer = TRUE, line = -1)
      }
      plot.type <- plot.type[plot.type != "model.resids"]
      cat(paste("plot type = \"model.resids\" Model residuals\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("std.model.resids" %in% plot.type) {
      # make plot of standardized observation residuals
      df <- subset(tT.resids, tT.resids$name == "model")
      df$.std.resid[is.na(df$value)] <- NA
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.resids, na.rm = TRUE), max(.resids, na.rm = TRUE))
          plot(t, .std.resid,
            type = "p", xlab = "",
            ylab = "", ylim = ylims,
            col = plotpar$point.col, pch = plotpar$point.pch,
            cex = plotpar$point.size
          )
          title(plt)
          box()
          abline(h = 0, lty = 1)
          abline(h = -2, lty = 3, col = "blue")
          abline(h = 2, lty = 3, col = "blue")
        })
        mtext("Observation standardized residuals", side = 2, outer = TRUE, line = -1)
      }
      plot.type <- plot.type[plot.type != "std.model.resids"]
      cat(paste("plot type = \"std.model.resids\" Standardized model residuals\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("state.resids" %in% plot.type) {
      # make plot of process residuals; set form='marxss' to get process resids
      df <- subset(tT.resids, tT.resids$name == "state")
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
        mtext("State residuals, xtT - E[x]", side = 2, outer = TRUE, line = -1)
      }
      plot.type <- plot.type[plot.type != "state.resids"]
      cat(paste("plot type = \"state.resids\" State residuals\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

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

    if ("qqplot.model.resids" %in% plot.type) {
      # make plot of observation residuals
      df <- subset(tT.resids, tT.resids$name == "model")
      slope <- tapply(df$.std.resid, df$.rownames, slp)
      intercept <- tapply(df$.std.resid, df$.rownames, int)
      nY <- min(9, attr(x$model, "model.dims")$y[1])
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          qqnorm(.std.resid, main = plt)
          abline(a = intercept[plt], b = slope[plt], col = plotpar$line.col, lwd = plotpar$line.lwd)
        })
      }
      plot.type <- plot.type[plot.type != "qqplot.model.resids"]
      cat(paste("plot type = \"qqplot.model.resids\" QQ plot of model standardized smoothed residuals\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("qqplot.state.resids" %in% plot.type) {
      # make qqplot of state residuals
      df <- subset(tT.resids, tT.resids$name == "state")
      df$.rownames <- paste0("State ", df$.rownames)
      slope <- tapply(df$.std.resid, df$.rownames, slp)
      intercept <- tapply(df$.std.resid, df$.rownames, int)
      nX <- min(9, attr(x$model, "model.dims")$x[1])
      plot.nrow <- round(sqrt(nX))
      plot.ncol <- ceiling(nX / plot.nrow)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          qqnorm(.std.resid, main = plt)
          abline(a = intercept[plt], b = slope[plt], col = plotpar$line.col, lwd = plotpar$line.lwd)
        })
      }
      plot.type <- plot.type[plot.type != "qqplot.state.resids"]
      cat(paste("plot type = \"qqplot.state.resids\" QQ plot of state standardized smoothed residuals\n"))
      if (length(plot.type) != 0) {
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
          lines(t, estimate, col = plotpar$line.col, lwd = plotpar$line.lwd)
          box()
        })
      }
      plot.type <- plot.type[plot.type != "ytT"]
      cat(paste("plot type = \"ytT\" Expected value of Y conditioned on data\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("acf.model.resids" %in% plot.type) {
      df <- subset(tt1.resids, tt1.resids$name == "model")
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
      plot.type <- plot.type[plot.type != "acf.model.resids"]
      cat(paste("plot type = \"acf.model.resids\" Model innovations residuals ACF\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    if ("acf.state.resids" %in% plot.type) {
      df <- subset(tt1.resids, tt1.resids$name == "state")
      df$.rownames <- paste0("State ", df$.rownames)
      nX <- min(9, attr(x$model, "model.dims")$x[1])
      plot.nrow <- round(sqrt(nX))
      plot.ncol <- ceiling(nX / plot.nrow)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        with(subset(df, df$.rownames == plt), {
          stats::acf(.resids, na.action = na.pass, main = "")
          title(tit)
          box()
        })
      }
      plot.type <- plot.type[plot.type != "acf.state.resids"]
      cat(paste("plot type = \"acf.state.resids\" State one-step ahead residuals ACF\n"))
      if (length(plot.type) != 0) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }
  }
