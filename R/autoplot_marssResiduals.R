autoplot.marssResiduals <-
  function(x,
           plot.type = c("model.resids", "state.resids", "qqplot.model.resids", "qqplot.state.resids", "acf.model.resids"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE,
           plot.par = list(),
           silent = FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for autoplot.marssMLE. Please install it.", call. = FALSE)
    }
    if (!inherits(x, "marssResiduals")) {
      stop("autoplot.marssResiduals: x must be class marssResiduals.", call. = FALSE)
    }
    plot.type <- match.arg(plot.type, several.ok = TRUE)
    if (!("model" %in% x$name)) plot.type <- plot.type[!stringr::str_detect(plot.type, "model")]
    if (!("state" %in% x$name)) plot.type <- plot.type[!stringr::str_detect(plot.type, "state")]
    if (length(plot.type) == 0) {
      message("Nothing to plot. Either your MARSSresiduals object does not include model (ytt1) or state (xtT) residuals or you have passed in the wrong plot.type, i.e. model residual plots when your MARRSSresiduals object only includes state residuals (xtT).")
      return()
    }

    # Argument checks
    if (!is.numeric(conf.level) || length(conf.level) > 1 || conf.level > 1 || conf.level < 0) {
      stop("autoplot.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
    }
    if (!(conf.int %in% c(TRUE, FALSE))) stop("autoplot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)

    # Various warnings about residuals properties

    if (attr(x, "standardization") != "Cholesky") {
      if (!silent) message("Note: the residuals standardization is not Cholesky (the default for the residuals.marssMLE() function). See discussion of residual standardization and correlation in state-space residuals in ?MARSSresiduals.\n\n")
    }

    if ("model" %in% x$name && "ytT" %in% x$type) {
      if (!silent) message("Note: the model residuals with type equal to 'ytT' are smoothation residuals. The ACF will be shown for these but note that smoothation model residuals are not temporally independent. If you are checking model residuals for temporal independence, use innovation residuals (type='tt1' which is the default for residuals.marssMLE().\n\n")
    }

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

    alpha <- 1 - conf.level
    plts <- list()

    if ("model.resids" %in% plot.type) {
      # make plot of observation residuals
      resids.types <- unique(subset(x, x$name == "model")$type)
      for (resids.type in resids.types) {
        df <- subset(x, x$name == "model" & x$type == resids.type)
        df$.rownames <- factor(df$.rownames) # drop levels
        p1 <- ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$value)), ], ggplot2::aes_(~t, ~.resids)) +
          ggplot2::geom_point(
            shape = plotpar$point.pch, fill = plotpar$point.fill,
            col = plotpar$point.col, size = plotpar$point.size
          ) +
          ggplot2::xlab("Time") +
          ggplot2::ylab("Observation residuals, y - E[y]") +
          ggplot2::facet_wrap(~.rownames, scale = "free_y") +
          ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
          ggplot2::ggtitle("Model residuals")
        if (decorate) {
          p1 <- p1 + ggplot2::stat_smooth(formula = y ~ x, method = "loess", se = FALSE, na.rm = TRUE)
          df$sigma <- df$.sigma
          df$sigma[is.na(df$value)] <- 0
          df$ymin.resid <- qnorm(alpha / 2) * df$sigma
          df$ymax.resid <- -qnorm(alpha / 2) * df$sigma
          p1 <- p1 +
            ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin.resid, ymax = ~ymax.resid), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
        }
        plts[[paste0("model.resids.", resids.type)]] <- p1
      }
      if (identical(plot.type, "model.resids") && length(resids.types) == 1) {
        return(p1)
      }
    }

    if ("state.resids" %in% plot.type) {
      # make plot of process residuals; set form='marxss' to get process resids
      df <- subset(x, x$name == "state")
      df$.rownames <- factor(df$.rownames) # drop levels
      df$.rownames <- paste0("State ", df$.rownames)
      p1 <- ggplot2::ggplot(df[!is.na(df$.resids), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size
        ) +
        ggplot2::xlab("Time") +
        ggplot2::ylab("State residuals, xtT - E[x]") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle("State residuals")
      if (decorate) {
        p1 <- p1 + ggplot2::stat_smooth(formula = y ~ x, method = "loess", se = FALSE, na.rm = TRUE)
        df$ymin.resid <- qnorm(alpha / 2) * df$.sigma
        df$ymax.resid <- -qnorm(alpha / 2) * df$.sigma
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin.resid, ymax = ~ymax.resid), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      plts[["state.resids"]] <- p1
      if (identical(plot.type, "state.resids")) {
        return(p1)
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
      resids.types <- unique(subset(x, x$name == "model")$type)
      for (resids.type in resids.types) {
        df <- subset(x, x$name == "model" & x$type == resids.type)
        df$.rownames <- factor(df$.rownames) # drop levels
        slope <- tapply(df$.std.resids, df$.rownames, slp)
        intercept <- tapply(df$.std.resids, df$.rownames, int)
        abline.dat <- data.frame(
          .rownames = names(slope),
          slope = slope,
          intercept = intercept,
          stringsAsFactors = FALSE
        )
        p1 <- ggplot2::ggplot(df) +
          ggplot2::geom_qq(ggplot2::aes_(sample = ~.std.resids), na.rm = TRUE) +
          ggplot2::xlab("Theoretical quantiles") +
          ggplot2::ylab("Cholesky standardized model smoothed residuals") +
          ggplot2::facet_wrap(~.rownames, scale = "free_y") +
          ggplot2::ggtitle("Model residuals normality test")
        if (decorate) p1 <- p1 + ggplot2::geom_abline(data = abline.dat, ggplot2::aes_(slope = ~slope, intercept = ~intercept), color = "blue")
        plts[[paste0("qqplot.model.resids.", resids.type)]] <- p1
      }
      if (identical(plot.type, "qqplot.model.resids") && length(resids.types) == 1) {
        return(p1)
      }
    }

    if ("qqplot.state.resids" %in% plot.type) {
      # make qqplot of state residuals
      df <- subset(x, x$name == "state")
      df$.rownames <- factor(df$.rownames) # drop levels
      df$.rownames <- paste0("State ", df$.rownames)
      slope <- tapply(df$.std.resids, df$.rownames, slp)
      intercept <- tapply(df$.std.resids, df$.rownames, int)
      abline.dat <- data.frame(
        .rownames = names(slope),
        slope = slope,
        intercept = intercept,
        stringsAsFactors = FALSE
      )
      p1 <- ggplot2::ggplot(df) +
        ggplot2::geom_qq(ggplot2::aes_(sample = ~.std.resids), na.rm = TRUE) +
        ggplot2::xlab("Theoretical quantiles") +
        ggplot2::ylab("Cholesky standardized state smoothed residuals") +
        ggplot2::facet_wrap(~.rownames, scales = "free_y") +
        ggplot2::ggtitle("State residuals normality test")
      if (decorate) p1 <- p1 + ggplot2::geom_abline(data = abline.dat, ggplot2::aes_(slope = ~slope, intercept = ~intercept), color = "blue")
      plts[["qqplot.state.resids"]] <- p1
      if (identical(plot.type, "qqplot.state.resids")) {
        return(p1)
      }
    }

    # ACF functions
    acffun <- function(x, y) {
      bacf <- acf(x, plot = FALSE, na.action = na.pass)
      bacfdf <- with(bacf, data.frame(lag, acf))
      return(bacfdf)
    }
    acfci <- function(x) {
      ciline <- qnorm((1 - conf.level) / 2) / sqrt(length(x))
      return(ciline)
    }
    if ("acf.model.resids" %in% plot.type) {
      resids.types <- unique(subset(x, x$name == "model")$type)
      for (resids.type in resids.types) {
        df <- subset(x, x$name == "model" & x$type == resids.type)
        df$.rownames <- factor(df$.rownames) # drop state levels
        acfdf <- tapply(df$.resids, df$.rownames, acffun)
        fun <- function(x, y) {
          data.frame(.rownames = y, lag = x$lag, acf = x$acf, stringsAsFactors = FALSE)
        }
        acfdf <- mapply(fun, acfdf, names(acfdf), SIMPLIFY = FALSE)
        acf.dat <- data.frame(
          .rownames = unlist(lapply(acfdf, function(x) {
            x$.rownames
          })),
          lag = unlist(lapply(acfdf, function(x) {
            x$lag
          })),
          acf = unlist(lapply(acfdf, function(x) {
            x$acf
          })),
          stringsAsFactors = FALSE
        )

        cidf <- tapply(df$.resids, df$.rownames, acfci)
        ci.dat <- data.frame(.rownames = names(cidf), ci = cidf, stringsAsFactors = FALSE)

        p1 <- ggplot2::ggplot(acf.dat, mapping = ggplot2::aes(x = lag, y = acf)) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
          ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0), na.rm = TRUE) +
          ggplot2::xlab("Lag") +
          ggplot2::ylab("ACF") +
          ggplot2::facet_wrap(~.rownames, scales = "free_y") +
          ggplot2::ggtitle("Model innovations residuals ACF")
        p1 <- p1 +
          ggplot2::geom_hline(data = ci.dat, ggplot2::aes_(yintercept = ~ci), color = "blue", linetype = 2) +
          ggplot2::geom_hline(data = ci.dat, ggplot2::aes_(yintercept = ~ -ci), color = "blue", linetype = 2)
        plts[[paste0("acf.model.resids.", resids.type)]] <- p1
      }
      if (identical(plot.type, "acf.model.resids") && length(resids.types) == 1) {
        return(p1)
      }
    }
    last.plot <- names(plts)[length(names(plts))]
    for (i in names(plts)) {
      print(plts[[i]])
      if (!silent) cat(paste("plot.type =", i, "\n"))
      if (i != last.plot && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      } else {
        if (!silent) cat("Finished plots.\n")
      }
    }
  }
