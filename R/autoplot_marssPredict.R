autoplot.marssPredict <-
  function (x, include, pi.int = TRUE, decorate = TRUE,
            plot.par = list(),  ...) 
  {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")", 
           call. = FALSE)
    }
    if (!inherits(x, "marssPredict")) {
        stop("autoplot.marssPredict requires a marssPredict object.", call.=FALSE)
    }
    plotpar <- list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
                    line.col = "black", line.size = 1, line.linetype = "solid",
                    ci.fill = NULL, ci.col = NULL, 
                    ci.linetype = "blank", 
                    ci.linesize = 0, ci.alpha = 0.6, 
                    f.col = "#0000AA", f.linetype = "solid", f.linesize=0.5,
                    theme = theme_bw()
                    )
    if (!is.list(plot.par)) stop("autoplot.marssMLE: plot.par must be a list.", call. = FALSE)
    if (!missing(plot.par)){
      if (!all(names(plot.par) %in% names(plotpar))){ 
        stop(paste0("autoplot.marssMLE: Allowed plot.par names are ", paste(names(plotpar), collapse=", "), ".\n"), call. = FALSE) } else {
          for( i in names(plot.par))
            plotpar[[i]] <- plot.par[[i]]
        }
    }

    h <- x$h
    df <- x$pred
    lev <- sort(x$level)
    nx <- length(x$t)-h # length of pre-forecast data
    idx <- rev(order(lev))
    nint <- length(lev)
    if (missing(include) || include > nx) {
      include <- nx
    }
    if( !is.numeric(include) || include < 1 || (include %% 1) != 0 )
      stop("auotplot.marssPredict: include must be positive integer greater than 0.\n", call.=FALSE)
    if ( is.null(lev) || nint == 0 ) {
      pi.int <- FALSE
    }
    else if (!is.finite(max(df[,-1*c(1:4)]))) {
      pi.int <- FALSE
    }
    if (pi.int && (is.null(plotpar[["ci.fill"]]) || length(plotpar[["ci.fill"]]) != nint)) {
      if (min(lev) < 50) {
        plotpar[["ci.fill"]] <- colorspace::sequential_hcl(100)[lev]
      }
      else {
        plotpar[["ci.fill"]] <- colorspace::sequential_hcl(52)[lev - 49]
      }
    }
    if(is.null(plotpar[["ci.col"]]) || length(plotpar[["ci.col"]]) != nint)
      plotpar[["ci.col"]] <- plotpar[["ci.fill"]]
    
    if(h == 0){
      df <- subset(x$pred, t > nx-include+1)
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~estimate)) + plotpar[["theme"]]
      loc <- which(colnames(df)=="se")
      if (pi.int){
        for(i in length(lev):1){
          tmp <- df
          colnames(tmp)[colnames(tmp)==paste("Lo",lev[i])] <- "conf.low"
          colnames(tmp)[colnames(tmp)==paste("Hi",lev[i])] <- "conf.high"
          p1 <- p1 + ggplot2::geom_ribbon(data = tmp, ggplot2::aes_(ymin = ~conf.low, ymax = ~conf.high), 
                                          alpha = plotpar$ci.alpha, fill = plotpar$ci.fill[i], color = plotpar$ci.col[i], 
                                          linetype = plotpar$ci.linetype, size = plotpar$ci.linesize)
        }
      }
      p1 <- p1 +
        ggplot2::geom_line(linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$line.size) +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle(paste(ifelse(x$type=="xtT", "State", "Data"), "Predictions"))
      if(x$type=="ytT" && decorate) 
          p1 <- p1 + ggplot2::geom_point(data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y), 
                                         shape = plotpar$point.pch, fill = plotpar$point.fill, 
                                         col = plotpar$point.col, size = plotpar$point.size)
      return(p1)
    }
    if(h != 0){
      df <- subset(x$pred, t > nx-include+1)
      p1 <- ggplot2::ggplot(data = df, ggplot2::aes_(~t, ~estimate)) + plotpar[["theme"]]
      tmp <- subset(df, t <= nx)
      p1 <- p1 + ggplot2::geom_line(data= tmp, linetype = plotpar$line.linetype, color = plotpar$line.col, size = plotpar$f.linesize)
      loc <- which(colnames(df)=="se")
      if (pi.int){
        for(i in length(lev):1){
          tmp <- subset(df, t > nx) 
          colnames(tmp)[colnames(tmp)==paste("Lo",lev[i])] <- "conf.low"
          colnames(tmp)[colnames(tmp)==paste("Hi",lev[i])] <- "conf.high"
          p1 <- p1 + ggplot2::geom_ribbon(data = tmp, ggplot2::aes_(ymin = ~conf.low, ymax = ~conf.high), 
                                          alpha = plotpar$ci.alpha, fill = plotpar$ci.fill[i])
          p1 <- p1 + ggplot2::geom_line(data= tmp, linetype = plotpar$f.linetype, color = plotpar$f.col, size = plotpar$line.size)
          
        }
      }
      p1 <- p1 +
        ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
        ggplot2::facet_wrap(~.rownames, scale = "free_y") +
        ggplot2::ggtitle(paste(ifelse(x$type=="xtT", "State", "Data"), "Predictions"))
      if(x$type=="ytT" && decorate) 
        p1 <- p1 + ggplot2::geom_point(data = df[!is.na(df$y), ], ggplot2::aes_(~t, ~y), 
                                       shape = plotpar$point.pch, fill = plotpar$point.fill, 
                                       col = plotpar$point.col, size = plotpar$point.size)
      return(p1)
    }
    
  }