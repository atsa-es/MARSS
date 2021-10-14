\name{plot.marssMLE}
\alias{plot.marssMLE}
\alias{autoplot.marssMLE}

\title{ Plot MARSS MLE objects }
\description{
  Plots fitted observations and estimated states with confidence intervals using base R graphics (\code{plot}) and ggplot2 (\code{autoplot}).  Diagnostic plots also shown.  By default a subset of standard diagnostic plots are plotted.  Individual plots can be plotted by passing in \code{plot.type}.  If an individual plot is made using \code{autoplot()}, the ggplot object is returned which can be further manipulated.
}
\usage{
\method{plot}{marssMLE}(x, plot.type = c(
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
    "residuals", "all"),
    form=c("marxss", "marss", "dfa"),
    standardization = c("Cholesky", "marginal", "Block.Cholesky"),
    conf.int=TRUE, conf.level=0.95, decorate=TRUE, pi.int = FALSE,
    plot.par = list(), ..., silent = FALSE)
\method{autoplot}{marssMLE}(x, plot.type = c(
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
    "residuals", "all"),
    form=c("marxss", "marss", "dfa"),
    standardization = c("Cholesky", "marginal", "Block.Cholesky"),
    conf.int=TRUE, conf.level=0.95, decorate=TRUE, pi.int = FALSE,
    fig.notes = TRUE, plot.par = list(), ..., silent = FALSE)
  }
\arguments{
  \item{x}{ A \code{\link{marssMLE}} object.  }
  \item{plot.type}{ Type of plot.  If not passed in, a subset of the standard plots are drawn. See details for plot types.  }
  \item{standardization}{ The type of standardization to be used plots, if the user wants to specify a specific standardization. Otherwise Cholesky standardization is used. }
  \item{form}{ Optional. Form of the model.  This is normally taken from the form attribute of the MLE object (x), but the user can specify a different form. }
  \item{conf.int}{ TRUE/FALSE. Whether to include a confidence interval. }
  \item{pi.int}{ TRUE/FALSE. Whether to include a prediction interval on the observations plot }
  \item{conf.level}{ Confidence level for CIs. }
  \item{decorate}{ TRUE/FALSE. Add smoothing lines to residuals plots or qqline to qqplots and add data points plus residuals confidence intervals to states and observations plots. }
  \item{plot.par}{ A list of plot parameters to adjust the look of the plots.  See details. }
  \item{fig.notes}{ Add notes to the bottom of the plots (only for \code{autoplot()}. }
  \item{silent}{ No console interaction or output. }
  \item{...}{ Other arguments, not used. }
  }
\value{
  \code{autoplot()} will invisibly return the list of ggplot2 plot objects. Use \code{plts <- autoplot()} to obtain that list.
}
\details{
The plot types are as follows:
\describe{
\item{\code{"fitted.y"}}{ This plots the fitted y, which is the expected value of Y conditioned on the data from \eqn{t=1} to \eqn{t-1}, \eqn{t} or \eqn{T}. It is \eqn{Z*xtT + a}. The data are plotted for reference but note that the lines and intervals are for new data not the observed data.}
\item{\code{"fitted.x"}}{ This plots the fitted x, which is the expected value of X conditioned on the data from \eqn{t=1} to \eqn{t-1} or \eqn{T}. It is \eqn{B*E[x(t-1)|y] + u}. The \eqn{E[x(t)|y]} are plotted for reference but note that the lines and intervals are for new X. This is not the estimated states; these are used for residuals calculations. If you want the state estimates use \code{xtT} (or \code{xtt}).}
\item{\code{"xtT"}}{  The estimated states from the Kalman smoother (conditioned on all the data). }
\item{\code{"xtt1"}}{  The estimated states conditioned on the data up to \eqn{t-1}. Kalman filter output. }
\item{\code{"model.resids.ytT"}, \code{"model.resids.ytt1"}, \code{"model.resids.ytt"}}{Model residuals (data minus fitted y). \code{ytT} indicates smoothation residuals, \code{ytt1} indicates innovation residuals (the standard state-space residuals), and \code{ytt} are the residuals conditioned on data up to \eqn{t}. }
\item{\code{"state.resids.xtT"}}{ State smoothation residuals (E(x(t) | xtT(t-1)) minus xtT(t)). The intervals are the CIs for the smoothation residuals not one-step-ahead residuals. }
\item{\code{"std"}}{ \code{std} in front of any of the above plot names indicates that the plots are for the standardized residuals. }
\item{\code{"qqplot"}}{ Visual normality test for the residuals, model or state. } 
\item{\code{"acf"}}{ ACF of the residuals. The only residuals that should be temporally independent are the innovation residuals: \code{acf.model.residuals.ytt1} and \code{acf.std.model.residuals.ytt1}. This ACF is a standard residuals diagnostic for state-space models. The other ACF plots will show temporal dependence and are not used for diagnostics. }
\item{\code{"ytT"}}{ The expected value of \eqn{Y} conditioned on all the data. Use this for estimates of the missing data points. Note for non-missing \eqn{y} values, the expected value of \eqn{Y} is \eqn{y}. } 
\item{\code{"ytt"}, \code{ytt1}}{ The expected value of \eqn{Y} conditioned on the data from 1 to \eqn{t} or \eqn{t-1}. } 
}

The plot parameters can be passed in as a list to change the look of the plots. For \code{plot.marssMLE()}, the default is \code{plot.par = list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1, line.col = "black", line.size = 1, line.linetype = "solid", ci.col = "grey70", ci.border = NA, ci.lwd = 1, ci.lty = 1)}. For \code{autoplot.marssMLE}, the default is \code{plot.par = list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1, line.col = "black", line.size = 1, line.linetype = "solid", ci.fill = "grey70", ci.col = "grey70", ci.linetype = "solid", ci.linesize = 0, ci.alpha = 0.6)}.
}
\author{ 
  Eric Ward and Eli Holmes
}
\examples{
data(harborSealWA)
model.list <- list( Z = as.factor(c(1, 1, 1, 1, 2)), R = "diagonal and equal")
fit <- MARSS(t(harborSealWA[, -1]), model = model.list)
plot(fit, plot.type = "fitted.ytT")

require(ggplot2)
autoplot(fit, plot.type = "fitted.ytT")

\dontrun{
# DFA example
dfa <- MARSS(t(harborSealWA[, -1]), model = list(m = 2), form = "dfa")
plot(dfa, plot.type = "xtT")
}
}
