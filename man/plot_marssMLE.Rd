\name{plot.marssMLE}
\alias{plot}
\alias{plot.marssMLE}
\title{ plot function for MARSS MLE objects }
\description{
  Plots fitted observations and estimated states with confidence intervals.  Diagnostic plots also shown.  By default all plots are plotted.  Individual plots can be plotted by passing in \code{type.plot}.  If an individual plot is made, the ggplot object is returned which can be further manipulated.
}
  \usage{
  \method{plot}{marssMLE}(x, 
           plot.type=c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqnorm", "state.residuals.qqnorm"), 
           form=c("marxss", "marss", "dfa"), ...)
  }
  \arguments{
  \item{x}{ A marssMLE object.  }
  \item{plot.type}{ Type of plot.  If not passed in, all plots are drawn. Options for arguments include 'observations' (fits to the raw data), 'states' (estimates of the hidden or latent trends), 'model.residuals' (residuals for the observation error), 'state.residuals' (residuals associated with the process model), 'model.residuals.qqnorm' (qq plot for the observation residuals), 'state.residuals.qqnorm' (qq plot for the state residuals)}
  \item{form}{ Optional. Form of the model.  This is normally taken from the form attribute of the MLE object (x), but the user can specify a different form.}
  \item{...}{ Other arguments, not used. }
  }
\value{
  If an individual plot is selected using \code{plot.type} then the ggplot object is returned invisibly.
}
\author{ 
  Eric Ward and Eli Holmes
}
\examples{
data(harborSealWA)
kemfit <- MARSS(t(harborSealWA[,-1]), model=list(Z=as.factor(c(1,1,1,1,2)), R="diagonal and equal"))
plot(kemfit)   
 
# DFA example   
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
plt = plot(dfa, type="states")
}
