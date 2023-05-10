autoplot.marssPredict <-
  function(x, include, decorate = TRUE,
           plot.par = list(), ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")",
        call. = FALSE
      )
    }
    if (!inherits(x, "marssPredict")) {
      stop("autoplot.marssPredict requires a marssPredict object.", call. = FALSE)
    }
    if(!is.logical(decorate)){
      stop("autoplot.marssPredict: decorate must be TRUE/FALSE.", call. = FALSE)
    }
    pi.int <- ifelse(decorate, TRUE, FALSE)
    plotpar <- list(
      point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
      line.col = "black", line.size = 1, line.type = "solid",
      ci.fill = NULL, ci.col = NULL,
      ci.linetype = "blank",
      ci.linesize = 0, ci.alpha = 0.6,
      f.col = "#0000AA", f.linetype = "solid", f.linesize = 0.5,
      theme = ggplot2::theme_bw()
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

    h <- x$h
    df <- x$pred
    lev <- sort(x$level)
    nx <- length(x$t) - h # length of pre-forecast data
    idx <- rev(order(lev))
    nint <- length(lev)
    if (missing(include) || include > nx) {
      include <- nx
    }
    if (!is.numeric(include) || include < 1 || (include %% 1) != 0) {
      stop("auotplot.marssPredict: include must be positive integer greater than 0.\n", call. = FALSE)
    }
    if (is.null(lev) || nint == 0) {
      pi.int <- FALSE
    }
    else if (!is.finite(max(df[, -1 * c(1:4)]))) {
      pi.int <- FALSE
    }
    # Palette from colorspace::sequential_hcl(n) and kept in sysdata.rda
    if (pi.int && (is.null(plotpar[["ci.fill"]]) || length(plotpar[["ci.fill"]]) != nint)) {
      if (min(lev) < 50) {
        plotpar[["ci.fill"]] <- hcl_palette_100[lev]
      }
      else {
        plotpar[["ci.fill"]] <- hcl_palette_52[lev - 49]
      }
    }
    if (is.null(plotpar[["ci.col"]]) || length(plotpar[["ci.col"]]) != nint) {
      plotpar[["ci.col"]] <- plotpar[["ci.fill"]]
    }

    if (h == 0) {
      df <- subset(x$pred, x$pred$t >= nx - include + 1)
      df$.rownames <- factor(df$.rownames, levels=unique(df$.rownames))
      # Set up the plot
      plottitle <- paste(ifelse(substr(x$type, 1, 1)=="x", "State", "Data"), x$type, "Predictions")
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~estimate)) + plotpar[["theme"]]
      if (pi.int) {
        for (i in length(lev):1) {
          tmp <- df
          colnames(tmp)[colnames(tmp) == paste("Lo", lev[i])] <- "conf.low"
          colnames(tmp)[colnames(tmp) == paste("Hi", lev[i])] <- "conf.up"
          p1 <- p1 + ggplot2::geom_ribbon(
            data = tmp, ggplot2::aes_(ymin = ~conf.low, ymax = ~conf.up),
            alpha = plotpar$ci.alpha, fill = plotpar$ci.fill[i], color = plotpar$ci.col[i],
            linetype = plotpar$ci.linetype, size = plotpar$ci.linesize
          )
        }
        plottitle <- paste(plottitle, switch(x$interval.type, none="", confidence="+ CI", prediction = "+ PI"))
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.type, color = plotpar$line.col, linewidth = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~ df$.rownames, scale = "free_y") +
        ggplot2::ggtitle(plottitle)
      if (substr(x$type, 1, 1)=="y" && decorate) {
        p1 <- p1 + ggplot2::geom_point(ggplot2::aes_(~t, ~y),
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size,
          na.rm = TRUE
        )
      }
      note <- NULL
      if(!all(unlist(lapply(x$newdata, is.null)))){
        tmp <- names(x$newdata)[!unlist(lapply(x$newdata, function(x){ all(is.na(x))}))]
        tmp <- tmp[tmp != "t"]
        note <- paste0("Prediction used newdata: ", paste0(tmp, collapse=", "), ". ")
      }
      if(substr(x$type, 1, 1)=="x"){
        note <- paste(note, "State (x) prediction is the expected value of x(t) conditioned on x(t-1) where x(t-1) is computed conditioned on all the data (xtT) or data up to t-1 (xtt1). This type of state prediction is not commonly used. If you are looking for the estimated states (with CIs), then use autoplot(fit, plot.type='xtT').")
      }
      if(!is.null(note)){
        p1 <- p1 + 
          ggplot2::labs(caption = paste0(strwrap(note, width=grDevices::dev.size(units = "px")[1]/6), collapse = "\n")) + 
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      
      return(p1)
    }
    if (h != 0) {
      df <- subset(x$pred, x$t >= nx - include + 1)
      df$.rownames <- factor(df$.rownames, levels=unique(df$.rownames))
      # set up the plot
      plottitle <- paste(ifelse(substr(x$type, 1, 1)== "x", "State", "Data"), x$type, "Forecasts")
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~estimate)) + plotpar[["theme"]]
      if (pi.int) {
        for (i in length(lev):1) {
          tmp <- df
          colnames(tmp)[colnames(tmp) == paste("Lo", lev[i])] <- "conf.low"
          colnames(tmp)[colnames(tmp) == paste("Hi", lev[i])] <- "conf.up"
          tmp$conf.low[tmp$t <= nx] <- NA
          tmp$conf.up[tmp$t <= nx] <- NA
          p1 <- p1 +
            ggplot2::geom_ribbon(
              data = tmp, ggplot2::aes_(ymin = ~conf.low, ymax = ~conf.up),
              alpha = plotpar$ci.alpha, fill = plotpar$ci.fill[i], na.rm = TRUE
            )
        }
        plottitle <- paste(plottitle, switch(x$interval.type, none="", 
                                             confidence = paste0("+ ", paste0(lev, collapse=", "), "% CI"), 
                                             prediction = paste0("+ ", paste0(lev, collapse=", "), "% PI")))
      }
      # make the estimate line
      tmp <- df
      tmp$estimate[tmp$t > nx] <- NA # subset makes facet_wrap fail
      p1 <- p1 + ggplot2::geom_line(data = tmp, 
                                    linetype = plotpar$line.type, 
                                    color = plotpar$line.col, 
                                    linewidth = plotpar$line.size, na.rm = TRUE)
      # Add the forecast line
      tmp <- df
      tmp$estimate[tmp$t <= nx] <- NA
      p1 <- p1 + ggplot2::geom_line(data = tmp, linetype = plotpar$f.linetype, color = plotpar$f.col, linewidth = plotpar$f.linesize, na.rm = TRUE)
      # Add labels and titles
      p1 <- p1 +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~ df$.rownames, scale = "free_y") +
        ggplot2::ggtitle(plottitle)
      
      # Add the data points
      if (substr(x$type, 1, 1) == "y" && decorate) {
        p1 <- p1 + ggplot2::geom_point(
          data = df, ggplot2::aes_(~t, ~y),
          shape = plotpar$point.pch, fill = plotpar$point.fill,
          col = plotpar$point.col, size = plotpar$point.size, na.rm = TRUE
        )
      }
      
      if(!all(unlist(lapply(x$newdata, is.null)))){
        tmp <- names(x$newdata)[!unlist(lapply(x$newdata, function(x){ all(is.na(x))}))]
        tmp <- tmp[tmp != "t"]
        note <- paste("Forecast used newdata in the forecast period:", paste0(tmp, collapse=", "))
      p1 <- p1 + 
        ggplot2::labs(caption = paste0(strwrap(note, width=grDevices::dev.size(units = "px")[1]/6), collapse = "\n")) + 
        ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.5, hjust = 0))
      }
      
      
      return(p1)
    }
  }
