\name{plot.marssPredict}
\alias{plot.marssPredict}
\alias{autoplot.marssPredict}

\title{ Plot MARSS Forecast and Predict objects }
\description{
  Plots forecasts with prediction (default) or confidence intervals using base R graphics (\code{plot}) and ggplot2 (\code{autoplot}).  The plot function is built to mimic \code{plot.forecast} in the forecast package in terms of arguments and look.
}
\usage{
\method{plot}{marssPredict}(x, include, decorate = TRUE, main = NULL, showgap = TRUE, 
    shaded = TRUE, shadebars = (x$h < 5 & x$h != 0), shadecols = NULL, col = 1, 
    fcol = 4, pi.col = 1, pi.lty = 2, ylim = NULL, 
    xlab = "", ylab = "", type = "l", flty = 1, flwd = 2, ...)
\method{autoplot}{marssPredict}(x, include, decorate = TRUE, plot.par = list(),  ...)
  }
  \arguments{
\item{x}{ marssPredict produced by \code{\link{forecast.marssMLE}} or \code{\link{predict.marssMLE}}. }
\item{include}{ number of time step from the training data to include before the forecast. Default is all values.}
\item{ main }{ Text to add to plot titles. }
\item{ showgap }{ If showgap=FALSE, the gap between the training data and the forecasts is removed.}
\item{shaded}{ Whether prediction intervals should be shaded (TRUE) or lines (FALSE).}
\item{shadebars}{ Whether prediction intervals should be plotted as shaded bars (if TRUE) or a shaded polygon (if FALSE). Ignored if shaded=FALSE. Bars are plotted by default if there are fewer than five forecast horizons.}
\item{shadecols}{ Colors for shaded prediction intervals. }
\item{col}{ Color for the data line.}
\item{fcol}{ Color for the forecast line.}
\item{pi.col}{ If shaded=FALSE and PI=TRUE, the prediction intervals are plotted in this color.}
\item{pi.lty}{ If shaded=FALSE and PI=TRUE, the prediction intervals are plotted using this line type.}
\item{ylim}{ Limits on y-axis.}
\item{xlab}{ X-axis label. }
\item{ylab}{ Y-axis label. }
\item{type}{ Type of plot desired. As for plot.default.}
\item{flty}{ Line type for the forecast line. }
\item{flwd}{ Line width for the forecast line. }
\item{...}{ Other arguments, not used. }
\item{decorate}{ TRUE/FALSE. Add data points and CIs or PIs to the plots. }
\item{plot.par}{ A list of plot parameters to adjust the look of the plot.  The default is 
               list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
                    line.col = "black", line.size = 1, line.type = "solid",
                    ci.fill = NULL, ci.col = NULL, 
                    ci.linetype = "blank", 
                    ci.linesize = 0, ci.alpha = 0.6, 
                    f.col = "#0000AA", f.linetype = "solid", f.linesize=0.5,
                    theme = theme_bw()
                    ). }
  }
\value{
  None. Plots are plotted
}
\author{ 
  Eli Holmes and based off of \code{plot.forecast} in the forecast package written by Rob J Hyndman & Mitchell O'Hara-Wild.
}
\seealso{ \code{\link{predict.marssMLE}()} }
\examples{
data(harborSealWA)
dat <- t(harborSealWA[, -1])
fit <- MARSS(dat[1:2,])
fr <- predict(fit, n.ahead=10)
plot(fr, include=10)

# forecast.marssMLE does the same thing as predict with h
fr <- forecast(fit, n.ahead=10)
plot(fr)

# without h, predict will show the prediction intervals
fr <- predict(fit)
plot(fr)

# you can fit to a new set of data using the same model and same x0
fr <- predict(fit, newdata=list(y=dat[3:4,]), x0="use.model")
plot(fr)

# but you probably want to re-estimate x0
fr <- predict(fit, newdata=list(y=dat[3:4,]), x0="reestimate")
plot(fr)

# forecast; note h not n.ahead is used for forecast()
fr <- forecast(fit, h=10)
}
