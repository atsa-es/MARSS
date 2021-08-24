autoplot.marssMLE <-
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
    if (identical(plot.type,"all")){
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
    plts <- list()

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
      df <- tsSmooth.marssMLE(x, type = "ytT", ifelse(conf.int, "confidence", "none"), level = conf.level)
      if (conf.int) {
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

    resids.vs.time.plots <- paste0(rep(c("model.resids", "std.model.resids"), each = 3), c(".ytT", ".ytt", ".ytt1"))
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1] # ytT, ytt or ytt1
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
      df$.rownames <- factor(df$.rownames) # drop levels
      p1 <- ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$value)), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size
        ) +
        ggplot2::xlab("Time") +
        ggplot2::ylab(ifelse(stringr::str_detect(i, "std"), "Standardized observation residuals, y - E[y]", "Observation residuals, y - E[y]")) +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle(ifelse(stringr::str_detect(i, "std"), "Standardized model residuals", "Model residuals"))
      if (decorate) {
        p1 <- p1 + ggplot2::stat_smooth(formula = y ~ x, method = "loess", se = FALSE, na.rm = TRUE)
        if (stringr::str_detect(i, "std")) {
          df$sigma <- 1
          if (i == "std.model.resids.ytt1") note <- "Cholesky standardized innovations residuals. Use standardized model smoothation model residuals ('std.model.resids.ytT') for outlier detection."
          if (i == "std.model.resids.ytt") note <- "Cholesky standardized ytt residuals. Use standardized model smoothation model residuals ('std.model.resids.ytT') for outlier detection."
          if (i == "std.model.resids.ytT") note <- "Cholesky standardized model smoothation residuals. Residuals outside the +/- 2 limits are potential outliers."
        } else {
          df$sigma <- df$.sigma
          df$sigma[is.na(df$value)] <- 0
          if (i == "model.resids.ytt1") note <- "Innovations (one-step ahead) residuals. (1-alpha) fraction of residuals should fall within the CIs."
          if (i == "model.resids.ytt") note <- "Model residuals conditioned on data up to t. (1-alpha) fraction of residuals should fall within the CIs."
          if (i == "model.resids.ytT") note <- "Model smoothation residuals (y - E[y|all data]). Note, standardized smoothation model residuals are the more typical outlier diagnostic."
        }
        df$ymin.resid <- stats::qnorm(alpha / 2) * df$sigma
        df$ymax.resid <- -stats::qnorm(alpha / 2) * df$sigma
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin.resid, ymax = ~ymax.resid), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize) +
          ggplot2::labs(caption = note) + ggplot2::theme(plot.caption = element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        return(p1)
      }
    }

    resids.vs.time.plots <- c("state.resids.xtT", "std.state.resids.xtT")
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      # make plot of process residuals; set form='marxss' to get process resids
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1]
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
      df$.rownames <- factor(df$.rownames) # drop levels
      df$.rownames <- paste0("State ", df$.rownames)
      p1 <- ggplot2::ggplot(df[!is.na(df$.resids), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size
        ) +
        ggplot2::xlab("Time") +
        ggplot2::ylab(ifelse(stringr::str_detect(i, "std"), "Standardized state residuals, xtT - E[x]", "State residuals, xtT - E[x]")) +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle(ifelse(stringr::str_detect(i, "std"), "Standardized state residuals", "State residuals"))
      if (decorate) {
        p1 <- p1 + ggplot2::stat_smooth(formula = y ~ x, method = "loess", se = FALSE, na.rm = TRUE)
        if (stringr::str_detect(i, "std")) {
          df$sigma <- 1
          note <- "Cholesky standardized state smoothation residuals. Residuals outside the +/- 2 limits are potential outliers of x(t) to x(t+1)."
        } else {
          df$sigma <- df$.sigma
          df$sigma[is.na(df$value)] <- 0
          note <- "State smoothation residuals. Note, standardized state smoothation residuals are the more typical outlier diagnostic."
        }
        df$ymin.resid <- stats::qnorm(alpha / 2) * df$sigma
        df$ymax.resid <- -stats::qnorm(alpha / 2) * df$sigma
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin.resid, ymax = ~ymax.resid), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize) +
          ggplot2::labs(caption = note) + ggplot2::theme(plot.caption = element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        return(p1)
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
      df$.rownames <- factor(df$.rownames) # drop levels
      slope <- tapply(df$.resids, df$.rownames, slp)
      intercept <- tapply(df$.resids, df$.rownames, int)
      abline.dat <- data.frame(
        .rownames = names(slope),
        slope = slope,
        intercept = intercept,
        stringsAsFactors = FALSE
      )
      ylab.val <- paste0(
        ifelse(stringr::str_detect(i, "std"), "cholesky standardized ", ""),
        cname, " ",
        switch(ctype,
          xtT = "smoothation ",
          ytt1 = "innovation ",
          ytt = "ytt ",
          ytT = "smoothation "
        ),
        "residuals"
      )
      ylab.val <- stringr::str_to_sentence(ylab.val)
      p1 <- ggplot2::ggplot(df) +
        ggplot2::geom_qq(ggplot2::aes_(sample = ~.resids), na.rm = TRUE) +
        ggplot2::xlab("Theoretical quantiles") +
        ggplot2::ylab(ylab.val) +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle("Residuals normality test")
      if (decorate) p1 <- p1 + ggplot2::geom_abline(data = abline.dat, ggplot2::aes_(slope = ~slope, intercept = ~intercept), color = "blue")
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        return(p1)
      }
    }

    # ACF plots
    acffun <- function(x, y) {
      bacf <- acf(x, plot = FALSE, na.action = na.pass)
      bacfdf <- with(bacf, data.frame(lag, acf))
      return(bacfdf)
    }
    acfci <- function(x) {
      ciline <- qnorm((1 - conf.level) / 2) / sqrt(length(x))
      return(ciline)
    }
    acf.plots <- plot.type[stringr::str_detect(plot.type, "acf")]
    for (i in acf.plots[acf.plots %in% plot.type]) {
      ctype <- rev(stringr::str_split(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(stringr::str_detect(i, "model"), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (stringr::str_detect(i, "std")) df$.resids <- df$.std.resids
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

      title.val <- paste0(
        ifelse(stringr::str_detect(i, "std"), "Cholesky standardized ", ""),
        cname, " ",
        switch(ctype,
          xtT = "smoothation ",
          ytt1 = "innovation ",
          ytt = "ytt ",
          ytT = "smoothation "
        ),
        "residuals ACF"
      )
      title.val <- stringr::str_to_sentence(title.val)

      p1 <- ggplot2::ggplot(acf.dat, mapping = ggplot2::aes(x = lag, y = acf)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
        ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0), na.rm = TRUE) +
        ggplot2::xlab("Lag") +
        ggplot2::ylab("ACF") +
        ggplot2::facet_wrap(~.rownames, scales = "free_y") +
        ggplot2::ggtitle(title.val)
      p1 <- p1 +
        ggplot2::geom_hline(data = ci.dat, ggplot2::aes_(yintercept = ~ci), color = "blue", linetype = 2) +
        ggplot2::geom_hline(data = ci.dat, ggplot2::aes_(yintercept = ~ -ci), color = "blue", linetype = 2)
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        return(p1)
      }
    }

    for (i in plot.type) {
      print(plts[[i]])
      if (!silent) cat(paste("plot.type =", i, "\n"))
      if (i != plot.type[length(plot.type)] && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      } else {
        if (!silent) cat("Finished plots.\n")
      }
    }
  }
