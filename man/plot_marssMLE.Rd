\name{plot.marssMLE}
\alias{plot.marssMLE}
\alias{autoplot.marssMLE}

\title{ Plot MARSS MLE objects }
\description{
  Plots fitted observations and estimated states with confidence intervals using base R graphics (\code{plot}) and ggplot2 (\code{autoplot}).  Diagnostic plots also shown.  By default all plots are plotted.  Individual plots can be plotted by passing in \code{type.plot}.  If an individual plot is made using \code{autoplot}, the ggplot object is returned which can be further manipulated.
}
\usage{
\method{plot}{marssMLE}(x, plot.type = c("model.ytT", "xtT", "model.resids",
                 "state.resids", "qqplot.model.resids",
                 "qqplot.state.resids", "ytT", "acf.model.resids",
                 "acf.state.resids"), form = c("marxss", "marss",
                 "dfa"), conf.int = TRUE, conf.level = 0.95, 
                 decorate = TRUE, pi.int = FALSE, plot.par = list(), ...)
\method{autoplot}{marssMLE}(x, plot.type = c("model.ytT", "xtT", "model.resids", 
    "state.resids", "qqplot.model.resids", "qqplot.state.resids", "ytT", 
    "acf.model.resids", "acf.state.resids"),
    form=c("marxss", "marss", "dfa"),
    conf.int=TRUE, conf.level=0.95, decorate=TRUE, pi.int = FALSE,
    plot.par = list(), ...)
  }
\arguments{
  \item{x}{ A \code{\link{marssMLE}} object.  }
  \item{plot.type}{ Type of plot.  If not passed in, all plots are drawn. See details for plot types.  }
  \item{form}{ Optional. Form of the model.  This is normally taken from the form attribute of the MLE object (x), but the user can specify a different form. }
  \item{conf.int}{ TRUE/FALSE. Whether to include a confidence interval. }
  \item{pi.int}{ TRUE/FALSE. Whether to include a prediction interval on the observations plot }
  \item{conf.level}{ Confidence level for CIs. }
  \item{decorate}{ TRUE/FALSE. Add smoothing lines to residuals plots or qqline to qqplots and add data points plus residuals confidence intervals to states and observations plots. }
  \item{plot.par}{ A list of plot parameters to adjust the look of the plots.  The default is list(point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
                    line.col = "black", line.size = 1, line.linetype = "solid",
                    ci.fill = "grey70", ci.col = "grey70", ci.linetype = "solid", 
                    ci.linesize = 0, ci.alpha = 0.6). }
  \item{...}{ Other arguments, not used. }
  }
\value{
  If an individual plot is selected using \code{plot.type} and \code{autoplot} is called, then the ggplot object is returned invisibly.
}
\details{
The plot types are as follows:
\describe{
\item{\code{"model.ytT"}}{ \code{"observations"} can also be used. This plots the fitted ytT, which is the expected value of new data. It is \eqn{Z*xtT + a}. The data are plotted for reference but note that the lines and intervals are for new data not the observed data.}
\item{\code{"xtT"}}{  \code{"states"} can also be used. The estimates states from the Kalman smoother (conditioned on all the data). }
\item{\code{"model.resids"}}{ \code{"model.residuals"} can also be used. Model smoothation residuals (data minus fitted ytT). The intervals are the CIs for the smoothation residuals not one-step-ahead (aka innovations) residuals. }
\item{\code{"state.resids"}}{ \code{"state.residuals"} can also be used. State smoothation residuals (E(x(t) | xtT(t-1)) minus xtT(t)). The intervals are the CIs for the smoothation residuals not one-step-ahead residuals. }
\item{\code{qqplot.model.resids}}{\code{"model.residuals.qqplot"} can also be used. Visual normality test for the model residuals.} 
\item{\code{qqplot.state.resids}}{\code{"state.residuals.qqplot"} can also be used. Visual normality test for the state residuals.} 
\item{\code{"ytT"}}{ The expected value of \eqn{y} conditioned on the data. Use this for estimates of the missing data points. Note for non-missing \eqn{y} values, the expected value of \eqn{y} is \eqn{y}. } 
\item{\code{"acf.model.resids"}}{ ACF of the model residuals.}
\item{\code{"acf.state.resids"}}{ ACF of the state residuals.}
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
