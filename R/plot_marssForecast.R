plot.marssForecast <- function(x, include, PI = TRUE, showgap = TRUE, shaded = TRUE, shadebars = (x$h < 5), shadecols = NULL, col = 1, 
          fcol = 4, pi.col = 1, pi.lty = 2, ylim = NULL, main = NULL, 
          xlab = "", ylab = "", type = "l", flty = 1, flwd = 2, ...) 
{
  nx <- x$T
  h <- x$h
  xf <- x$forecast
  idx <- rev(order(x$level))
  nint <- length(x$level)
  n <- length(unique(xf$.rownames))
  if (missing(include) || include > nx) {
    include <- nx
  }
  if (length(x$level) == 0 || is.null(x$level)) {
    PI <- FALSE
  }
  else if (!is.finite(max(xf[,-1*c(1:4)]))) {
    PI <- FALSE
  }
  if (!shaded) {
    shadebars <- FALSE
  }
  if (is.null(shadecols)) {
    if (min(x$level) < 50) {
      shadecols <- rev(colorspace::sequential_hcl(100)[x$level])
    }
    else {
      shadecols <- rev(colorspace::sequential_hcl(52)[x$level - 49])
    }
  }
  if (length(shadecols) == 1) {
    if (shadecols == "oldstyle") {
      shadecols <- heat.colors(nint + 2)[switch(1 + (nint > 1), 2, nint:1) + 1]
    }
  }
  
  nX <- min(9, n)
  plot.nrow <- round(sqrt(nX))
  plot.ncol <- ceiling(nX / plot.nrow)
  par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
  for (plt in unique(xf$.rownames)) {
    if (is.null(main)) {
      if(x$type=="xtT") main <- paste("Forecasts for state ", plt, sep = "")
      if(x$type=="ytT") main <- paste("Forecasts for ", plt, sep = "")
    }
    tmp <- subset(xf, xf$.rownames == plt)
      if (is.null(ylim))
        ylim <- c(min(tmp[,-1*c(1:4)], na.rm = TRUE), max(tmp[,-1*c(1:4)], na.rm = TRUE))
      xxx <- 1:(nx+h)
      if(!showgap) xxx <- c(1:nx, nx:(nx+h-1))
      plot(xxx, c(tmp$estimate[(nx - include + 1):nx], rep(NA, h)), xlab = "", ylab = "Estimate", ylim = ylim, main=main, col=col, type=type, ...)
      if(x$type=="ytT") points(xxx[(nx - include + 1):nx], tmp$y[(nx - include + 1):nx], col = fcol, pch = 19, cex=0.75)
      if (PI){
        xxx <- (nx+1):(nx+h)
        if(!showgap) xxx <- xxx-1
        loc <- which(colnames(tmp)==".sd")
        for (i in 1:nint) {
          if (shadebars) {
            for (j in 1:h) {
              polygon(xxx[j] + c(-0.5, 0.5, 0.5, -0.5), 
                      c(rep(tmp[nx+j, loc+(idx[i]-1)*2+1], 2), 
                        rep(tmp[nx+j, loc+1+(idx[i]-1)*2+1], 2)), 
                      col = shadecols[i], border = FALSE)
            }
          }
          else if (shaded) {
            polygon(c(xxx, rev(xxx)), c(tmp[(nx+1):(nx+h), loc+(idx[i]-1)*2+1], 
                                        rev(tmp[(nx+1):(nx+h), loc+1+(idx[i]-1)*2+1])),
                    col = shadecols[i], border = FALSE)
          }
          else if (h == 1) {
            lines(c(xxx) + c(-0.5, 0.5), rep(tmp[(nx+1):(nx+h), loc+(idx[i]-1)*2+1], 2), col = pi.col, lty = pi.lty)
            lines(c(xxx) + c(-0.5, 0.5), rep(tmp[(nx+1):(nx+h), loc+1+(idx[i]-1)*2+1], 2), col = pi.col, lty = pi.lty)
          }
          else {
            lines(xxx, tmp[(nx+1):(nx+h), loc+(idx[i]-1)*2+1], col = pi.col, lty = pi.lty)
            lines(xxx, tmp[(nx+1):(nx+h), loc+1+(idx[i]-1)*2+1], col = pi.col, lty = pi.lty)
          }
        }
      } # end PI
      if (h > 1 && !shadebars) {
        lines(xxx, tmp$estimate[(nx+1):(nx+h)], lty = flty, lwd = flwd, col = fcol)
      }
      else {
        points(xxx, tmp$estimate[(nx+1):(nx+h)], col = fcol, pch = 19)
      }
    }
  
    invisible(x$forecast)
}