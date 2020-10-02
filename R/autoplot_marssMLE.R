autoplot.marssMLE <-
  function(x,
           plot.type = c("model.ytT", "xtT", "model.resids", "state.resids", "qqplot.model.resids", "qqplot.state.resids", "ytT", "acf.model.resids"),
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
    # Argument checks
    plot.type <- match.arg(plot.type, several.ok = TRUE)
    old.plot.type <- c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqplot", "state.residuals.qqplot", "expected.value.observations")
    new.plot.type <- c("model.ytT", "xtT", "model.resids", "state.resids", "qqplot.model.resids", "qqplot.state.resids", "ytT")
    for (i in 1:NROW(old.plot.type)) if (old.plot.type[i] %in% plot.type) plot.type[plot.type == old.plot.type[i]] <- new.plot.type[i]
    if (!is.numeric(conf.level) || length(conf.level) > 1 || conf.level > 1 || conf.level < 0) stop("autoplot.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
    if (!(conf.int %in% c(TRUE, FALSE))) stop("autoplot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)

    if (missing(form)) {
      model_form <- attr(x[["model"]], "form")[1]
    } else {
      model_form <- match.arg(form)
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
    plts <- list()

    if (any(str_detect(plot.type, "acf"))) {
      tt1.resids <- residuals.marssMLE(x, type = "tt1", standardization = "Cholesky")
    }
    if (any(str_detect(plot.type, "resids"))) {
      tT.resids <- residuals.marssMLE(x, type = "tT", standardization = "Cholesky")
    }
    
    if ("xtT" %in% plot.type) {
      # make plot of states and CIs

      if ("rotate" %in% names(extras)) {
        rotate <- extras[["rotate"]]
        if (!(rotate %in% c(TRUE, FALSE))) stop("autoplot.marssMLE: rotate must be TRUE/FALSE. \n")
      } else {
        rotate <- FALSE
      }

      states <- tsSmooth.marssMLE(x, type = "xtT", ifelse(conf.int, "confidence", "none"), level = conf.level, ...)
      if (model_form == "dfa") {
        if (rotate) {
          rottext <- "rotated"
        } else {
          rottext <- ""
        }
        states$term <- paste("DFA", rottext, "trend", states$term)
      } else {
        states$term <- paste0("State ", states$term)
      }
      p1 <- ggplot2::ggplot(data = states, ggplot2::aes_(~t, ~.estimate))
      if (conf.int) {
        p1 <- p1 + ggplot2::geom_ribbon(data = states, ggplot2::aes_(ymin = ~.conf.low, ymax = ~.conf.up), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle("States")
      plts[["xtT"]] <- p1
      if (identical(plot.type, "xtT")) {
        return(p1)
      }
    }

    if ("model.ytT" %in% plot.type) {
      # make plot of observations
      tit <- "Model fitted Y"
      if (conf.int) tit <- paste(tit, "+ CI")
      if (pi.int) tit <- paste(tit, "+ PI (dashed)")
      df <- fitted.marssMLE(x, type = "ytT", interval = "confidence", level = conf.level)
      df$ymin <- df$.conf.low
      df$ymax <- df$.conf.up
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~.fitted))
      if (conf.int) {
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin, ymax = ~ymax), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      if (pi.int) {
        df2 <- fitted.marssMLE(x, type = "ytT", interval = "prediction", level = conf.level)
        df$ymin.pi <- df2$.lwr
        df$ymax.pi <- df2$.upr
        p1 <- p1 + ggplot2::geom_line(data = df, ggplot2::aes_(~t, ~ymin.pi), linetype = "dashed")
        p1 <- p1 + ggplot2::geom_line(data = df, ggplot2::aes_(~t, ~ymax.pi), linetype = "dashed")
      }
      # Add data points
      if (decorate) {
        p1 <- p1 + ggplot2::geom_point(
          data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y),
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size, na.rm = TRUE
        )
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle(tit)
      plts[["model.ytT"]] <- p1
      if (identical(plot.type, "model.ytT")) {
        return(p1)
      }
    }

    if ("ytT" %in% plot.type) {
      # make plot of expected value of Y condtioned on y(1)
      df <- tsSmooth.marssMLE(x, type = "ytT", ifelse(conf.int, "confidence", "none"), level=conf.level)
      if(conf.int){
        df$ymin <- df$.conf.low
        df$ymax <- df$.conf.up
      }
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~.estimate)) +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size)
      if (conf.int) {
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin, ymax = ~ymax), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_point(
          data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y),
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size
        ) +
        ggplot2::ggtitle("Expected value of Y conditioned on data")
      plts[["ytT"]] <- p1
      if (identical(plot.type, "ytT")) {
        return(p1)
      }
    }

    if ("model.resids" %in% plot.type) {
      # make plot of observation residuals
      df <- subset(tT.resids, tT.resids$name == "model")
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
      plts[["model.resids"]] <- p1
      if (identical(plot.type, "model.resids")) {
        return(p1)
      }
    }

    if ("state.resids" %in% plot.type) {
      # make plot of process residuals; set form='marxss' to get process resids
      df <- subset(tT.resids, tT.resids$name == "state")
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
      # make plot of observation residuals
      df <- subset(tT.resids, tT.resids$name == "model")
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
      plts[["qqplot.model.resids"]] <- p1
      if (identical(plot.type, "qqplot.model.resids")) {
        return(p1)
      }
    }

    if ("qqplot.state.resids" %in% plot.type) {
      # make qqplot of state residuals
      df <- subset(tT.resids, tT.resids$name == "state")
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
      df <- subset(tt1.resids, tt1.resids$name == "model")
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
      plts[["acf.model.resids"]] <- p1
      if (identical(plot.type, "acf.model.resids")) {
        return(p1)
      }
    }
    for (i in plot.type) {
      print(plts[[i]])
      cat(paste("plot.type =", i, "\n"))
      if (i != plot.type[length(plot.type)] && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      } else {
        if(!silent) cat("Finished plots.\n")
      }
    }
  }
