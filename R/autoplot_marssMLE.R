autoplot.marssMLE <-
  function(x,
           plot.type = c(
             "fitted.ytT", "fitted.ytt", "fitted.ytt1",
             "ytT", "ytt", "ytt1",
             "fitted.xtT", "fitted.xtt1",
             "xtT", "xtt", "xtt1",
             "model.resids.ytt1", "qqplot.model.resids.ytt1", "acf.model.resids.ytt1",
             "std.model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
             "model.resids.ytT", "qqplot.model.resids.ytT", "acf.model.resids.ytT",
             "std.model.resids.ytT", "qqplot.std.model.resids.ytT", "acf.std.model.resids.ytT",
             "model.resids.ytt", "qqplot.model.resids.ytt", "acf.model.resids.ytt",
             "std.model.resids.ytt", "qqplot.std.model.resids.ytt", "acf.std.model.resids.ytt",
             "state.resids.xtT", "qqplot.state.resids.xtT", "acf.state.resids.xtT",
             "std.state.resids.xtT", "qqplot.std.state.resids.xtT", "acf.std.state.resids.xtT",
             "residuals", "all"
           ),
           form = c("marxss", "marss", "dfa"),
           standardization = c("Cholesky", "marginal", "Block.Cholesky"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE, pi.int = FALSE,
           fig.notes = TRUE, plot.par = list(), ..., silent = FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for autoplot.marssMLE. Please install it.", call. = FALSE)
    }

    # Argument checks
    standardization <- match.arg(standardization)

    if (!is.numeric(conf.level) || length(conf.level) > 1 || conf.level > 1 || conf.level < 0) stop("autoplot.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)

    # Test logical arguments
    largs <- list(conf.int = conf.int, pi.int = pi.int, decorate = decorate, fig.notes = fig.notes, silent = silent)
    for (i in names(largs)) {
      if (!(identical(largs[[i]], TRUE) || identical(largs[[i]], FALSE))) stop(paste("autoplot.marssMLE:", i, "must be TRUE/FALSE"), call. = FALSE)
    }

    # Argument checks: plot.type
    if (missing(plot.type)) {
      plot.type <- c(
        "fitted.ytT", "xtT",
        "model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
        "std.model.resids.ytT",
        "std.state.resids.xtT", "qqplot.std.state.resids.xtT"
      )
    }
    plot.type <- match.arg.exact(plot.type, several.ok = TRUE)
    if (identical(plot.type, "residuals")) {
      plot.type <- c(
        "model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
        "std.model.resids.ytT", "std.state.resids.xtT", "qqplot.std.state.resids.xtT"
      )
    }
    plot.all <- FALSE
    if (identical(plot.type, "all")) {
      plot.all <- TRUE
      plot.type <- eval(formals()$plot.type)
      plot.type <- plot.type[!(plot.type %in% c("residuals", "all"))]
    }

    # Check class and alter plot.type as needed
    if (!inherits(x, "marssMLE")) {
      if (inherits(x, "marssResiduals")) {
        plot.type <- plot.type[grepl("resids", plot.type)]
        # Make sure that plot types are possible for the object that the user passed in
        ctype <- unique(x$type)
        plot.type <- plot.type[sapply(plot.type, function(x) {
          any(sapply(ctype, function(s) grepl(s, x)))
        })]
        cname <- unique(x$name)
        plot.type <- plot.type[sapply(plot.type, function(x) {
          any(sapply(cname, function(s) grepl(s, x)))
        })]
        if (length(plot.type) == 0) {
          message("Nothing to plot. Either your MARSSresiduals object does not include model or state residuals or you have passed in the wrong plot.type, i.e. model residual plots when your MARRSSresiduals object only includes state residuals.")
          return()
        }
        # Set up the residuals object
        resids <- x
        cstan <- attr(resids, "standardization")
      } else {
        stop("autoplot.marssMLE: x must be class marssMLE or marssResiduals.", call. = FALSE)
      }
    }

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

    # If user requests any residuals plots, set up the residuals data frames unless x is marssResiduals object
    if (!inherits(x, "marssResiduals")) {
      resids <- c()
      cstan <- standardization
      if (any(grepl("tt1", plot.type))) {
        resids <- residuals.marssMLE(x, type = "tt1", standardization = cstan)
      }
      if (any(grepl("tt", plot.type) & !grepl("tt1", plot.type))) {
        resids <- rbind(resids, residuals.marssMLE(x, type = "tt", standardization = cstan))
      }
      if (any(grepl("tT", plot.type))) {
        resids <- rbind(resids, residuals.marssMLE(x, type = "tT", standardization = cstan))
      }
    }

    state.plots <- c("xtT", "xtt", "xtt1")
    for (i in state.plots[state.plots %in% plot.type]) {
      ctype <- i
      if (model_form == "dfa" && "rotate" %in% names(extras)) {
        rotate <- extras[["rotate"]]
        if (!(rotate %in% c(TRUE, FALSE))) stop("autoplot.marssMLE: rotate must be TRUE/FALSE. \n")
      } else {
        rotate <- FALSE
      }
      df <- tsSmooth.marssMLE(x, type = i, ifelse(conf.int, "confidence", "none"), level = conf.level, ...)
      df$.rownames <- factor(df$.rownames, levels = unique(df$.rownames))
      rottext <- ifelse(rotate, "the rotated ", "")
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~.estimate))
      if (conf.int) {
        p1 <- p1 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~.conf.low, ymax = ~.conf.up), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle(paste(str_to_sentence(paste0(rottext, "States")), i))
      if (fig.notes) {
        note <- paste0("This is the estimate of ", rottext, "X conditioned on the data from t=1 to ", switch(ctype,
          xtT = "T.",
          xtt = "t.",
          xtt1 = "t-1."
        ), ifelse(ctype != "xtT", " If you want the smoothed state estimates conditioned on all the data, use xtT instead.", ""), ifelse(conf.int, paste(" Confidence intervals are for the expected value of X (conditioned on the data up to ", switch(ctype,
          xtT = "T).",
          xtt = "t).",
          xtt1 = "t-1)."
        )), ""))
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
      }
    }

    fitted.plots <- paste0("fitted.", c("ytt1", "ytt", "ytT", "xtT", "xtt1"))
    for (i in fitted.plots[fitted.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1]
      cname <- ifelse(grepl("y", i), "model", "state")
      tit <- paste("Fitted", ctype)
      if (conf.int) tit <- paste(tit, "+ CI")
      if (pi.int) tit <- paste(tit, "+ PI (dashed)")
      df <- fitted.marssMLE(x, type = ctype, interval = "confidence", level = conf.level)
      df$ymin <- df$.conf.low
      df$ymax <- df$.conf.up
      # Drop levels and make sure ggplot doesn't rearrange levels
      df$.rownames <- factor(df$.rownames, levels = unique(df$.rownames))

      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~.fitted))
      if (conf.int) {
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin, ymax = ~ymax), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      if (pi.int) {
        df2 <- fitted.marssMLE(x, type = ctype, interval = "prediction", level = conf.level)
        df$ymin.pi <- df2$.lwr
        df$ymax.pi <- df2$.upr
        p1 <- p1 + ggplot2::geom_line(data = df, ggplot2::aes_(~t, ~ymin.pi), linetype = "dashed")
        p1 <- p1 + ggplot2::geom_line(data = df, ggplot2::aes_(~t, ~ymax.pi), linetype = "dashed")
      }
      # Add data points if y
      if (decorate & grepl("y", i)) {
        p1 <- p1 + ggplot2::geom_point(
          data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y),
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size, na.rm = TRUE
        )
      }
      if (decorate & grepl("x", i)) {
        p1 <- p1 + ggplot2::geom_point(
          data = df, ggplot2::aes_(~t, ~.x),
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size, na.rm = TRUE
        )
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle(tit)
      if (fig.notes) {
        if (cname == "model") {
          note <- paste("This is the model fitted value of Y conditioned on the data from t=1 to", switch(ctype,
            ytT = "T.",
            ytt = "t.",
            ytt1 = "t-1."
          ), ifelse(ctype != "ytt1", "Use fitted.ytt1 if you want the one-step-ahead predictions instead.", "These are known as the one-step-ahead predictions."))
          if (decorate && conf.int && !pi.int) note <- paste(note, "The CI is for the expected value of Y and the data points will not fall within the CI. Use prediction intervals to compare the data to intervals.")
          if (decorate && conf.int && pi.int) note <- paste(note, "The CI is for the expected value of Y and the data points will not fall within the CI. Data points should fall with the PI.")
        }
        if (cname == "state") {
          note <- paste("This is the model fitted value of X conditioned on the data from t=1 to", switch(ctype,
            xtT = "T.",
            xtt = "t.",
            xtt1 = "t-1."
          ), "This is not the model estimate of X (i.e. the states). It is the expected value of X(t) given the E[X(t-1)|y] where y is the data from t=1 to", switch(ctype,
            xtT = "T.",
            xtt = "t.",
            xtt1 = "t-1."
          ), "Use plot.type='xtT' if you want the traditional states estimates for a state-space model. These estimates are used in state residuals calculations.")
        }
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
      }
    }

    y.plots <- c("ytT", "ytt", "ytt1")
    for (i in y.plots[y.plots %in% plot.type]) {
      ctype <- i
      # make plot of expected value of Y condtioned on y(1)
      if (ctype %in% c("ytT", "ytt1") | !conf.int) {
        df <- tsSmooth.marssMLE(x, type = i, interval = ifelse(conf.int, "confidence", "none"), level = conf.level)
      } else {
        if (plot.all) next # If plot.type="all" then just skip the problematic plots
        if (conf.int) message(paste("Confidence intervals for", i, "are not implemented in MARSS.\nNo confidence intervals shown for the", i, "plot.\n"))
        df <- tsSmooth.marssMLE(x, type = i, interval = "none")
      }
      # make sure that ggplot doesn't re-order the levels
      df$.rownames <- factor(df$.rownames, levels = unique(df$.rownames))
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~.estimate)) +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size)
      if (conf.int & ctype != "ytt") {
        df$ymin <- df$.conf.low
        df$ymax <- df$.conf.up
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
      if (fig.notes) {
        note <- paste("This is the estimate of Y conditioned on the data from t=1 to", switch(ctype,
          ytT = "T.",
          ytt = "t.",
          ytt1 = "t-1."
        ), ifelse(ctype == "ytT", "Use this if you need an estimate of missing data points. For non-missing data, it will simply return the observed y. E(Y|y) = y. If you want the model estimate of Y, use fitted.ytT.", "If you need an estimate of missing data, use ytT."))
        if (conf.int) note <- paste(note, "Confidence intervals are for the expected value of Y and will be 0 if data are not missing.")
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
      }
    }

    resids.vs.time.plots <- paste0(rep(c("model.resids", "std.model.resids"), each = 3), c(".ytT", ".ytt", ".ytt1"))
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1] # ytT, ytt or ytt1
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
      # drop x levels and make sure ggplot doesn't rearrange the levels
      df$.rownames <- factor(df$.rownames, levels = unique(df$.rownames))
      title.val <- paste0(
        ifelse(grepl("std", i), paste(cstan, "standardized "), ""),
        cname, " ",
        switch(ctype,
          xtT = "smoothation ",
          ytt1 = "innovation ",
          ytt = "ytt ",
          ytT = "smoothation "
        ),
        "residuals"
      )
      substr(title.val, 1, 1) <- toupper(substr(title.val, 1, 1))
      p1 <- ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$value)), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size
        ) +
        ggplot2::xlab("Time") +
        ggplot2::ylab(ifelse(grepl("std", i), "Standardized observation residuals, y - E[y]", "Observation residuals, y - E[y]")) +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle(title.val)
      if (decorate) {
        p1 <- p1 + ggplot2::stat_smooth(formula = y ~ x, method = "loess", se = FALSE, na.rm = TRUE)
      }
      if (conf.int) {
        if (grepl("std", i)) {
          df$sigma <- 1
        } else {
          df$sigma <- df$.sigma
          df$sigma[is.na(df$value)] <- 0 # will only apply for y; x value is never NA
        }
        df$ymin.resid <- stats::qnorm(alpha / 2) * df$sigma
        df$ymax.resid <- -stats::qnorm(alpha / 2) * df$sigma
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin.resid, ymax = ~ymax.resid), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      if (fig.notes) {
        if (grepl("std", i)) {
          cstan2 <- ifelse(cstan == "marginal", "Marginal", cstan)
          if (i == "std.model.resids.ytt1") note <- paste(cstan2, "standardized innovations residuals. Use standardized model smoothation (ytT) residuals (std.model.resids.ytT) for outlier detection.")
          if (i == "std.model.resids.ytt") note <- paste(cstan2, "standardized ytt residuals. Use standardized model smoothation (ytT) residuals (std.model.resids.ytT) for outlier detection.")
          if (i == "std.model.resids.ytT" & cstan == "Cholesky") note <- paste(cstan2, "standardized model smoothation (ytT) residuals. Residuals outside the +/- 2 limits are potential outliers.")
          if (i == "std.model.resids.ytT" & cstan != "Cholesky") note <- paste(cstan2, "standardized model smoothation (ytT) residuals. Use Cholesky standardized residuals for outlier detection.")
        } else {
          if (i == "model.resids.ytt1") note <- paste("Innovations (one-step ahead) residuals.", ifelse(conf.int, "(1-alpha) fraction of residuals should fall within the CIs. A violation of this indicates problems with the normality assumption.", ""))
          if (i == "model.resids.ytt") note <- "Model residuals conditioned on data up to t. These are not typically used. See model.resids.ytt1 and std.model.resids.ytT for more standard residuals diagnostics."
          if (i == "model.resids.ytT") note <- paste("Model smoothation (ytT) residuals (y - E[y|all data]). Note, Cholesky standardized model smoothation residuals are the more typical outlier diagnostic.", ifelse(conf.int, "(1-alpha) fraction of residuals should fall within the CIs. A violation of this indicates potential outliers.", ""))
        }
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
      }
    }

    resids.vs.time.plots <- c("state.resids.xtT", "std.state.resids.xtT")
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      # make plot of process residuals; set form='marxss' to get process resids
      ctype <- rev(strsplit(i, "[.]")[[1]])[1]
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
      # Drop y levels and make sure ggplot doesn't rearrange levels
      df$.rownames <- factor(df$.rownames, levels = unique(df$.rownames))
      title.val <- paste0(
        ifelse(grepl("std", i), paste(cstan, "standardized "), ""),
        cname, " ",
        switch(ctype,
          xtT = "smoothation ",
          ytt1 = "innovation ",
          ytt = "ytt ",
          ytT = "smoothation "
        ),
        "residuals"
      )
      substr(title.val, 1, 1) <- toupper(substr(title.val, 1, 1))
      p1 <- ggplot2::ggplot(df[!is.na(df$.resids), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size
        ) +
        ggplot2::xlab("Time") +
        ggplot2::ylab(ifelse(grepl("std", i), "Standardized state residuals, xtT - E[x]", "State residuals, xtT - E[x]")) +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle(title.val)
      if (decorate) {
        p1 <- p1 + ggplot2::stat_smooth(formula = y ~ x, method = "loess", se = FALSE, na.rm = TRUE)
      }
      if (conf.int) {
        if (grepl("std", i)) {
          df$sigma <- 1
        } else {
          df$sigma <- df$.sigma
          df$sigma[is.na(df$value)] <- 0 # will never be the case for x so not really needed here
        }
        df$ymin.resid <- stats::qnorm(alpha / 2) * df$sigma
        df$ymax.resid <- -stats::qnorm(alpha / 2) * df$sigma
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin.resid, ymax = ~ymax.resid), alpha = plotpar$ci.alpha, fill = plotpar$ci.fill, color = plotpar$ci.col, linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
      }
      if (fig.notes) {
        if (grepl("std", i)) {
          cstan2 <- ifelse(cstan == "marginal", "Marginal", cstan)
          note <- paste(cstan2, "standardized state smoothation (xtT) residuals.")
          if (cstan == "Cholesky" & conf.int) note <- paste(note, "Residuals outside the +/- 2 limits are potential outliers of x(t) to x(t+1).")
          if (cstan != "Cholesky" & conf.int) note <- paste(note, "Confidence intervals are for the x(t) to x(t+1) step.")
        } else {
          note <- "State smoothation (xtT) residuals. Note, Cholesky standardized state smoothation residuals are the more typical outlier diagnostic."
          if (conf.int) note <- paste(note, "Confidence intervals are for the x(t) to x(t+1) step.")
        }
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
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

    qqplot.plots <- plot.type[grepl("qqplot", plot.type)]
    for (i in qqplot.plots[qqplot.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
      # Drop levels and make sure ggplot doesn't rearrange levels
      df$.rownames <- factor(df$.rownames, levels = unique(df$.rownames))
      slope <- tapply(df$.resids, df$.rownames, slp)
      intercept <- tapply(df$.resids, df$.rownames, int)
      abline.dat <- data.frame(
        .rownames = names(slope),
        slope = slope,
        intercept = intercept,
        stringsAsFactors = FALSE
      )
      # Drop levels and make sure ggplot doesn't rearrange levels
      abline.dat$.rownames <- factor(abline.dat$.rownames, levels = unique(abline.dat$.rownames))
      ylab.val <- paste0(
        ifelse(grepl("std", i), "standardized ", ""),
        cname, " ",
        switch(ctype,
          xtT = "smoothation ",
          ytt1 = "innovation ",
          ytt = "ytt ",
          ytT = "smoothation "
        ),
        "residuals"
      )
      ylab.val <- str_to_sentence(ylab.val)
      p1 <- ggplot2::ggplot(df) +
        ggplot2::geom_qq(ggplot2::aes_(sample = ~.resids), na.rm = TRUE) +
        ggplot2::xlab("Theoretical quantiles") +
        ggplot2::ylab(ylab.val) +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle("Residuals normality test") +
        ggplot2::geom_abline(data = abline.dat, ggplot2::aes_(slope = ~slope, intercept = ~intercept), color = "blue")
      if (fig.notes) {
        note <- paste0(
          ifelse(grepl("std", i), paste(cstan, "standardized "), ""),
          cname, " ",
          switch(ctype,
            xtT = "smoothation (xtT) ",
            ytt1 = "innovation (ytt1) ",
            ytt = "ytt ",
            ytT = "smoothation (ytT) "
          ),
          "residuals. The residuals should be Gaussian."
        )
        note <- str_to_sentence(note, ignore = c("Cholesky", "Gaussian", "Block.Cholesky", "(ytT)", "(xtT)"))
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
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
    acf.plots <- plot.type[grepl("acf", plot.type)]
    for (i in acf.plots[acf.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
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
      acf.dat$.rownames <- factor(acf.dat$.rownames, levels = unique(df$.rownames))

      cidf <- tapply(df$.resids, df$.rownames, acfci)
      ci.dat <- data.frame(.rownames = names(cidf), ci = cidf, stringsAsFactors = FALSE)
      ci.dat$.rownames <- factor(ci.dat$.rownames, levels = unique(df$.rownames))


      title.val <- paste0(
        ifelse(grepl("std", i), paste(cstan, "standardized "), ""),
        cname, " ",
        switch(ctype,
          xtT = "smoothation ",
          ytt1 = "innovation ",
          ytt = "ytt ",
          ytT = "smoothation "
        ),
        "residuals ACF"
      )
      title.val <- str_to_sentence(title.val, ignore = c("Cholesky", "Block.Cholesky", "(ytT)", "(xtT)"))

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
      if (fig.notes) {
        note <- paste0(
          ifelse(grepl("std", i), paste(cstan, "standardized "), ""),
          cname, " ",
          switch(ctype,
            xtT = "smoothation (xtT) ",
            ytt1 = "innovation (ytt1) ",
            ytt = "ytt ",
            ytT = "smoothation (ytT) "
          ),
          "residuals.",
          ifelse(grepl("ytt1", i), " These residuals should be temporally uncorrelated.", " These residuals are not expected to be temporally uncorrelated. Use innovation (ytt1) residuals to check for temporal correlation in the residuals.")
        )
        note <- str_to_sentence(note, ignore = c("Cholesky", "Block.Cholesky", "(ytT)", "(xtT)"))
        p1 <- p1 + ggplot2::labs(caption = paste0(strwrap(note), collapse = "\n")) + ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      plts[[i]] <- p1
      if (length(plot.type) == 1) {
        invisible(plts)
      }
    }

    for (i in names(plts)) {
      if(inherits(try(print(plts[[i]]), silent=TRUE), "try-error")){
        if (!silent) cat(paste("plot.type =", i, "returned an error.\n"))
      }else{
        if (!silent) cat(paste("plot.type =", i, "\n"))
        print(plts[[i]])
      }
      if (i != plot.type[length(plot.type)] && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      } else {
        if (!silent) cat("Finished plots.\n")
      }
    }

    invisible(plts)
  }
