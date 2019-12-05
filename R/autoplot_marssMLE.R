autoplot.marssMLE <-
  function(x,
           plot.type = c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqplot", "state.residuals.qqplot", "expected.value.observations"),
           form = c("marxss", "marss", "dfa"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for plot.marssMLE. Please install it.", call. = FALSE)
    }
    
    # Argument checks
    plot.type <- match.arg(plot.type, several.ok = TRUE)
    if (!is.numeric(conf.level) || length(conf.level)>1 || conf.level > 1 || conf.level < 0) stop("plot.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
    if (!(conf.int %in% c(TRUE, FALSE))) stop("plot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)
    
    if (missing(form)) {
      model_form <- attr(x[["model"]], "form")[1]
    } else {
      model_form <- match.arg(form)
    }
    
    extras <- list()
    if (!missing(...)) {
      extras <- list(...)
      allowednames <- c("rotate", "method", "hessian.fun", "nboot")
      bad.names <- names(extras)[!(names(extras) %in% allowednames)]
      if (!all(names(extras) %in% allowednames)) stop(paste("plot.marssMLE:", paste(bad.names, collapse = " "), "is/are unknown argument(s). See ?tidy.marssMLE for allowed arguments.\n"), call. = FALSE)
    }
    if (model_form != "dfa" & "rotate" %in% names(extras)) {
      cat("plot.marssMLE: 'rotate' argument is ignored if form!='dfa'\n Pass in form='dfa' if your model is a DFA model, but the form \n attribute is not set (because you set up your DFA model manually).\n\n")
      rotate <- FALSE
    }
    # End Argument checks
    
    alpha <- 1 - conf.level
    plts <- list()
    
    if ("states" %in% plot.type) {
      # make plot of states and CIs
      
      if ("rotate" %in% names(extras)) {
        rotate <- extras[["rotate"]]
        if (!(rotate %in% c(TRUE, FALSE))) stop("tidy.marssMLE: rotate must be TRUE/FALSE. \n")
      } else {
        rotate <- FALSE
      }
      
      states <- tidy.marssMLE(x, type = "states", conf.int = conf.int, conf.level = conf.level, ...)
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
      p1 <- ggplot2::ggplot(data = states, ggplot2::aes_(~t, ~estimate))
      if (conf.int) p1 <- p1 + ggplot2::geom_ribbon(data = states, ggplot2::aes_(ymin = ~conf.low, ymax = ~conf.high), alpha = 0.3, col = "grey")
      p1 <- p1 +
        ggplot2::geom_line() +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle("States")
      plts[["states"]] <- p1
      if (identical(plot.type, "states")) {
        return(p1)
      }
    }
    
    if ("observations" %in% plot.type) {
      # make plot of observations
      df <- augment.marssMLE(x, "observations", form = model_form)
      df$ymin <- df$.fitted - qnorm(alpha / 2) * df$.se.fit
      df$ymax <- df$.fitted + qnorm(alpha / 2) * df$.se.fit
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~.fitted))
      if (conf.int) {
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin, ymax = ~ymax), alpha = 0.3, col = "grey")
      }
      p1 <- p1 +
        ggplot2::geom_line() +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_point(data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y), col = "blue") +
        ggplot2::ggtitle("Model fitted values for Y")
      plts[["observations"]] <- p1
      if (identical(plot.type, "observations")) {
        return(p1)
      }
    }
    
    if ("expected.value.observations" %in% plot.type) {
      # make plot of observation residuals
      df <- tidy.marssMLE(MLEobj, type = "observations", form = "marxss")
      df$ymin <- df$conf.low
      df$ymax <- df$conf.high
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~estimate)) +
        ggplot2::geom_line()
      if (conf.int) {
        p1 <- p1 +
          ggplot2::geom_ribbon(data = df, ggplot2::aes_(ymin = ~ymin, ymax = ~ymax), alpha = 0.3, col = "grey")
      }
      p1 <- p1 +
        ggplot2::geom_line() +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_point(data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y), col = "blue") +
        ggplot2::ggtitle("Expected value of Y conditioned on data")
      plts[["expected.value.observations"]] <- p1
      if (identical(plot.type, "expected.value.observations")) {
        return(p1)
      }
    }
    
    if ("model.residuals" %in% plot.type) {
      # make plot of observation residuals
      df <- augment.marssMLE(x, "observations", form = "marxss")
      p1 <- ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$y)), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(col = "blue") +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Observation residuals, y - E[y]") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle("Model residual")
      if (decorate) p1 <- p1 + ggplot2::stat_smooth(method = "loess", se = conf.int, level = conf.level, na.rm = TRUE)
      plts[["model.residuals"]] <- p1
      if (identical(plot.type, "model.residuals")) {
        return(p1)
      }
    }
    
    if ("state.residuals" %in% plot.type) {
      # make plot of process residuals; set form='marxss' to get process resids
      df <- augment.marssMLE(x, "states", form = "marxss")
      df$.rownames <- paste0("State ", df$.rownames)
      p1 <- ggplot2::ggplot(df[!is.na(df$.resids), ], ggplot2::aes_(~t, ~.resids)) +
        ggplot2::geom_point(col = "blue") +
        ggplot2::xlab("Time") +
        ggplot2::ylab("State residuals, xtT - E[x]") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = 3) +
        ggplot2::ggtitle("State residuals")
      if (decorate) p1 <- p1 + ggplot2::stat_smooth(method = "loess", se = conf.int, level = conf.level, na.rm = TRUE)
      plts[["state.residuals"]] <- p1
      if (identical(plot.type, "state.residuals")) {
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
    
    if ("model.residuals.qqplot" %in% plot.type) {
      # make plot of observation residuals
      df <- augment.marssMLE(x, "observations", form = "marxss")
      slope <- tapply(df$.std.resid, df$.rownames, slp)
      intercept <- tapply(df$.std.resid, df$.rownames, int)
      abline.dat <- data.frame(.rownames = names(slope), slope = slope, intercept = intercept)
      p1 <- ggplot2::ggplot(df) +
        ggplot2::geom_qq(ggplot2::aes_(sample = ~.std.resid), na.rm = TRUE) +
        ggplot2::xlab("Theoretical Quantiles") +
        ggplot2::ylab("Standardized Model Residuals") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle("Model residuals")
      if (decorate) p1 <- p1 + ggplot2::geom_abline(data = abline.dat, ggplot2::aes_(slope = ~slope, intercept = ~intercept), color = "blue")
      plts[["model.residuals.qqplot"]] <- p1
      if (identical(plot.type, "model.residuals.qqplot")) {
        return(p1)
      }
    }
    
    if ("state.residuals.qqplot" %in% plot.type) {
      # make qqplot of state residuals
      df <- augment.marssMLE(x, "states", form = "marxss")
      df$.rownames <- paste0("State ", df$.rownames)
      slope <- tapply(df$.std.resid, df$.rownames, slp)
      intercept <- tapply(df$.std.resid, df$.rownames, int)
      abline.dat <- data.frame(.rownames = names(slope), slope = slope, intercept = intercept)
      p1 <- ggplot2::ggplot(df) +
        ggplot2::geom_qq(ggplot2::aes_(sample = ~.std.resid), na.rm = TRUE) +
        ggplot2::xlab("Theoretical Quantiles") +
        ggplot2::ylab("Standardized State Residuals") +
        ggplot2::facet_wrap(~.rownames, scales = "free_y") +
        ggplot2::ggtitle("State residuals")
      if (decorate) p1 <- p1 + ggplot2::geom_abline(data = abline.dat, ggplot2::aes_(slope = ~slope, intercept = ~intercept), color = "blue")
      plts[["state.residuals.qqplot"]] <- p1
      if (identical(plot.type, "state.residuals.qqplot")) {
        return(p1)
      }
    }
    for (i in plot.type) {
      print(plts[[i]])
      if (i != plot.type[length(plot.type)]) {
        cat(paste("plot.type =", i, "\n"))
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      } else {
        cat("Finished plots.\n")
      }
    }
  }