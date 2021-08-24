\name{plot.marssMLE}
\alias{plot.marssMLE}
\alias{autoplot.marssMLE}

\title{ Plot MARSS MLE objects }
\description{
  Plots fitted observations and estimated states with confidence intervals using base R graphics (\code{plot}) and ggplot2 (\code{autoplot}).  Diagnostic plots also shown.  By default a subse of standard diagnostic plots are plotted.  Individual plots can be plotted by passing in \code{plot.type}.  If an individual plot is made using \code{autoplot()}, the ggplot object is returned which can be further manipulated.
}
\usage{
\method{plot}{marssMLE}(x, plot.type = c(
             "model.ytT", "xtT",
             "model.resids.ytt1", "qqplot.model.resids.ytt1", "acf.model.resids.ytt1",
             "std.model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
             "model.resids.ytT", "qqplot.model.resids.ytT",
             "std.model.resids.ytT", "qqplot.std.model.resids.ytT",
             "model.resids.ytt", "qqplot.model.resids.ytt",
             "std.model.resids.ytt", "qqplot.std.model.resids.ytt",
             "state.resids.xtT", "qqplot.state.resids.xtT",
             "std.state.resids.xtT", "qqplot.std.state.resids.xtT",
             "ytT", "residuals", "all"),
                 form = c("marxss", "marss", "dfa"), 
                 conf.int = TRUE, conf.level = 0.95, decorate = TRUE, 
                 pi.int = FALSE, plot.par = list(), ..., silent = FALSE)
\method{autoplot}{marssMLE}(x, plot.type = c(
             "model.ytT", "xtT",
             "model.resids.ytt1", "qqplot.model.resids.ytt1", "acf.model.resids.ytt1",
             "std.model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
             "model.resids.ytT", "qqplot.model.resids.ytT",
             "std.model.resids.ytT", "qqplot.std.model.resids.ytT",
             "model.resids.ytt", "qqplot.model.resids.ytt",
             "std.model.resids.ytt", "qqplot.std.model.resids.ytt",
             "state.resids.xtT", "qqplot.state.resids.xtT",
             "std.state.resids.xtT", "qqplot.std.state.resids.xtT",
             "ytT", "residuals", "all"),
    form=c("marxss", "marss", "dfa"),
    conf.int=TRUE, conf.level=0.95, decorate=TRUE, pi.int = FALSE,
    plot.par = list(), ..., silent = FALSE)
  }
\arguments{
  \item{x}{ A \code{\link{marssMLE}} object.  }
  \item{plot.type}{ Type of plot.  If not passed in, a subset of the standard plots are drawn. See details for plot types.  }
  \item{form}{ Optional. Form of the model.  This is normally taken from the form attribute of the MLE object (x), but the user can specify a different form. }
  \item{conf.int}{ TRUE/FALSE. Whether to include a confidence interval. }
  \item{pi.int}{ TRUE/FALSE. Whether to include a prediction interval on the observations plot }
  \item{conf.level}{ Confidence level for CIs. }
  \item{decorate}{ TRUE/FALSE. Add smoothing lines to residuals plots or qqline to qqplots and add data points plus residuals confidence intervals to states and observations plots. }
  \item{plot.par}{ A list of plot parameters to adjust the look of the plots.  The default is list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
                    line.col = "black", line.size = 1, line.linetype = "solid",
                    ci.fill = "grey70", ci.col = "grey70", ci.linetype = "solid", 
                    ci.linesize = 0, ci.alpha = 0.6). }
  \item{silent}{ No console interaction or output. }
  \item{...}{ Other arguments, not used. }
  }
\value{
  If an individual plot is selected using \code{plot.type} and \code{autoplot()} is called, then the ggplot object is returned invisibly.
}
\details{
The plot types are as follows:
\describe{
\item{\code{"model.ytT"}}{ This plots the fitted ytT, which is the expected value of new data. It is \eqn{Z*xtT + a}. The data are plotted for reference but note that the lines and intervals are for new data not the observed data.}
\item{\code{"xtT"}}{  The estimated states from the Kalman smoother (conditioned on all the data). }
\item{\code{"model.resids.ytT"}, \code{"model.resids.ytt1"}, \code{"model.resids.ytt"}}{Model residuals (data minus fitted y). \code{ytT} indicates smoothation residuals, \code{ytt1} indicates innovation residuals (the standard state-space residuals), and \code{ytt} are the residuals conditioned on data up to \eqn{t}. }
\item{\code{"state.resids.xtT"}}{ State smoothation residuals (E(x(t) | xtT(t-1)) minus xtT(t)). The intervals are the CIs for the smoothation residuals not one-step-ahead residuals. }
\item{\code{"std"}}{ \code{std} in front of any of the above plot names indicates that the plots are for the Cholesky standardized residuals. }
\item{\code{"qqplot.model.resids"}}{ Visual normality test for the model residuals. These can be shown for the standardized residuals and for any of the conditioning options (\code{ytT}, \code{ytt1}, \code{ytt}). } 
\item{\code{"qqplot.state.resids"}}{ Visual normality test for the state residuals. These can be shown for the standardized residuals also. } 
\item{\code{"ytT"}}{ The expected value of \eqn{y} conditioned on the data. Use this for estimates of the missing data points. Note for non-missing \eqn{y} values, the expected value of \eqn{y} is \eqn{y}. } 
\item{\code{"acf.model.resids.ytt1"}}{ ACF of the one-step-ahead model residuals. These can be shown for the standardized residuals also.}
}
}
\author{ 
  Eric Ward and Eli Holmes
}
\examples{
data(harborSealWA)
model.list <- list( Z = as.factor(c(1, 1, 1, 1, 2)), R = "diagonal and equal")
fit <- MARSS(t(harborSealWA[, -1]), model = model.list)
plot(fit, plot.type = "model.ytT")

require(ggplot2)
autoplot(fit, plot.type = "model.ytT")

\dontrun{
# DFA example
dfa <- MARSS(t(harborSealWA[, -1]), model = list(m = 2), form = "dfa")
plot(dfa, plot.type = "xtT")
}
}
