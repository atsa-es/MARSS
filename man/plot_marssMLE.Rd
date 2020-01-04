\name{plot.marssMLE}
\alias{plot}
\alias{plot.marssMLE}
\alias{autoplot}
\alias{autoplot.marssMLE}
\title{ plotting functions for MARSS MLE objects }
\description{
  Plots fitted observations and estimated states with confidence intervals using base R graphics (\code{plot}) and ggplot2 (\code{autoplot}).  Diagnostic plots also shown.  By default all plots are plotted.  Individual plots can be plotted by passing in \code{type.plot}.  If an individual plot is made using \code{autoplot}, the ggplot object is returned which can be further manipulated.
}
  \usage{
  \method{plot}{marssMLE}(x, 
           plot.type=c("observations", "states", "model.residuals", 
                       "state.residuals", "model.residuals.qqplot", 
                       "state.residuals.qqplot"), 
           form=c("marxss", "marss", "dfa"),
           conf.int=TRUE, conf.level=0.95, decorate=TRUE,
           plot.par = list(), ...)
  \method{autoplot}{marssMLE}(x, 
           plot.type=c("observations", "states", "model.residuals", 
                       "state.residuals", "model.residuals.qqplot", 
                       "state.residuals.qqplot", "expected.value.observations"),
           form=c("marxss", "marss", "dfa"),
           conf.int=TRUE, conf.level=0.95, decorate=TRUE,
           plot.par = list(), ...)
  }
  \arguments{
  \item{x}{ A marssMLE object.  }
  \item{plot.type}{ Type of plot.  If not passed in, all plots are drawn. Options for arguments include 'observations' (fits to the raw data), 'states' (estimates of the hidden or latent trends), 'model.residuals' (residuals for the observation error), 'state.residuals' (residuals associated with the process model), 'model.residuals.qqplot' (qq plot for the observation residuals), 'state.residuals.qqplot' (qq plot for the state residuals)}
  \item{form}{ Optional. Form of the model.  This is normally taken from the form attribute of the MLE object (x), but the user can specify a different form. }
  \item{conf.int}{ TRUE/FALSE. Whether to include a confidence interval. }
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
\author{ 
  Eric Ward and Eli Holmes
}
\examples{
data(harborSealWA)
fit <- MARSS(t(harborSealWA[, -1]), 
             model = list(
                      Z = as.factor(c(1, 1, 1, 1, 2)), 
                      R = "diagonal and equal"))
plot(fit, plot.type = "observations")

require(ggplot2)
autoplot(fit, plot.type = "observations")

\dontrun{
# DFA example
dfa <- MARSS(t(harborSealWA[, -1]), model = list(m = 2), form = "dfa")
plot(dfa, plot.type = "states")
}
}
