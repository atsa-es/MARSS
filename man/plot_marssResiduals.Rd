\name{plot.marssResiduals}
\alias{plot.marssResiduals}
\alias{autoplot.marssResiduals}

\title{ Plot MARSS marssResiduals objects }
\description{
  Plots residuals using the output from a \code{\link{residuals.marssResiduals}()} call. By default all available residuals plots are plotted.  Individual plots can be plotted by passing in \code{plot.type}.  If an individual plot is made using \code{autoplot()}, the ggplot object is returned which can be further manipulated. Plots are only shown for those residual types in the \code{\link{marssResiduals}} object.
}
\usage{
\method{plot}{marssMLE}(x, plot.type = c("model.resids",
                 "state.resids", "qqplot.model.resids",
                 "qqplot.state.resids", "acf.model.resids"), 
                 conf.int = TRUE, conf.level = 0.95, decorate = TRUE, 
                 pi.int = FALSE, plot.par = list(), silent = FALSE)
\method{autoplot}{marssMLE}(x, plot.type = c("model.resids", 
    "state.resids", "qqplot.model.resids", "qqplot.state.resids", 
    "acf.model.resids"),
    conf.int=TRUE, conf.level=0.95, decorate=TRUE, pi.int = FALSE,
    plot.par = list(), silent = FALSE)
  }
\arguments{
  \item{x}{ A \code{\link{marssResiduals}} object.  }
  \item{plot.type}{ Type of plot.  If not passed in, all plots are drawn. See details for plot types.  }
  \item{conf.int}{ TRUE/FALSE. Whether to include a confidence interval. }
  \item{pi.int}{ TRUE/FALSE. Whether to include a prediction interval on the observations plot }
  \item{conf.level}{ Confidence level for CIs. }
  \item{decorate}{ TRUE/FALSE. Add smoothing lines to residuals plots or qqline to qqplots and add data points plus residuals confidence intervals to states and observations plots. }
  \item{plot.par}{ A list of plot parameters to adjust the look of the plots.  The default is list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
                    line.col = "black", line.size = 1, line.linetype = "solid",
                    ci.fill = "grey70", ci.col = "grey70", ci.linetype = "solid", 
                    ci.linesize = 0, ci.alpha = 0.6). }
  \item{silent}{ No console interaction or output. }
  }
\value{
  If an individual plot is selected using \code{plot.type} and \code{autoplot()} is called, then the ggplot object is returned invisibly.
}
\details{
If \code{resids <- residuals(x)} is used (default) where \code{x} is a \code{marssMLE} object from a \code{MARSS()} call, then \code{resids} has the innovations residuals, or one-step-ahead residuals. These are what are commonly used for residuals diagnostics in state-space modeling. However, other types of residuals are possible for state-space models; see \code{\link{MARSSresiduals}()} for details. The plot function for \code{marssResiduals} objects will handle all types of residuals that might be in the \code{marssResiduals} object. However if you simply use the default behavior, \code{resids <- residuals(x)} and \code{plot(resids)}, you will get the standard model residuals diagnostics plots for state-space models, i.e. only model residuals plots and only plots for innovations model residuals (no smoothations model residuals).

The plot types are as follows:
\describe{
\item{\code{"model.resids"}}{ Model residuals (data minus fitted values) are plotted versus time. The intervals are the CIs for the smoothation model residuals if \code{x$type is "ytT"} and the data are not expected to fall inside these CIS. If \code{x$type is "ytt1"}, these are one-step-ahead (aka innovations) residuals and the data would be expected to fall in the CIs. }
\item{\code{"state.resids"}}{ State smoothation residuals (E(x(t) | xtT(t-1)) minus xtT(t)). The intervals are the CIs for the smoothation residuals not one-step-ahead residuals. }
\item{\code{"qqplot.model.resids"}}{ Visual normality test for the model residuals.} 
\item{\code{"qqplot.state.resids"}}{ Visual normality test for the state residuals.} 
\item{\code{"acf.model.resids"}}{ ACF of the model residuals. If \code{x$type is "ytT"}, then these are smoothation model residuals and they are not temporally independent. If \code{x$type is "ytt1"}, these are the one-step-ahead (aka innovations) residuals and they should be temporally independent. }
}
}
\author{ 
  Eli Holmes
}
\examples{
data(harborSealWA)
model.list <- list( Z = as.factor(c(1, 1, 1, 1, 2)), R = "diagonal and equal")
fit <- MARSS(t(harborSealWA[, -1]), model = model.list)
resids <- residuals(fit)

require(ggplot2)
# plot versus time, QQ-norm plot, and ACF
autoplot(resids)

# only the ACF plot
autoplot(resids, plot.type = "acf.model.resids")

}
