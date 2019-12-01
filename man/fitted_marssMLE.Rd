\name{fitted.marssMLE}
\alias{fitted}
\alias{fitted.marssMLE}
\title{ fitted function for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs marssMLE objects.  \code{fitted(marssMLE)}, where marssMLE is the output from a \code{\link{MARSS}} call,   will return the modeled value of \eqn{y(t)} or \eqn{x(t)}.  See details.
}
  \usage{
  \method{fitted}{marssMLE}(object, ..., one.step.ahead=FALSE, type=c("observations", "states"))
  }
  \arguments{
  \item{object}{ A marssMLE object.  }
  \item{...}{ Other arguments, not used. }
  \item{one.step.ahead}{ If TRUE, expected value of the states is computed conditioned on the data up to t-1}
  \item{type}{ Fitted values for the observations or the states.}
  }
\value{
  A T column matrix of the fitted values with one row for each observation (or state) time series.
}
\details{
\strong{observation fitted values}

The model predicted (fitted) y at time t is \eqn{Z(t)\hat{x}(t)+a(t)}{Z(t)\hat{x}(t)+a(t)}, where the model is written in marss form. See \code{\link{MARSS.marxss}} for a discussion of the conversion of MARSS models with covariates (c and d) into marss form (which is how models are written in the internal MARSS algorithms).
  
\eqn{\hat{x}(t)} is the expected value of the states at time t. If \code{one.step.ahead=FALSE}, \eqn{\hat{x}(t)} is the expected value conditioned on all the data, i.e. xtT returned by \code{\link{MARSSkf}}.  If \code{one.step.ahead=TRUE}, then the predictions of the observations or states at time t use only the data up to time t-1.  These are commonly known as the one step ahead predictions for a state-space model.

\strong{state fitted values}

The model predicted \eqn{\hat{x}(t)} given \eqn{\hat{x}(t-1)} is \eqn{B(t)\hat{x}(t-1)+u(t)}, where the model is written in marss form. This type of state fitted value is used for process outlier detection and shock detection.  See \code{\link{residuals.marssMLE}} and read the references cited.
  
\eqn{\hat{x}(t-1)} used in the prediction is the expected value of the states at time t-1. If \code{one.step.ahead=FALSE}, \eqn{\hat{x}(t-1)} is the expected value conditioned on all the data, i.e. xtT returned by \code{\link{MARSSkf}}.  If \code{one.step.ahead=TRUE}, \eqn{\hat{x}(t-1)} is the expected value conditioned on the data up to time t-1, i.e. xtt (returned by \code{\link{MARSSkf}}) evaluated at time t-1.  

If you want estimates of the states, then you'll want either the states estimate conditioned on all the data (or conditioned on the data up to time t-1 or up to time t).  These are returned by \code{\link{MARSSkf}} in \code{xtT} (and \code{xtt1} and \code{xtt}) respectively.  Which one you want depends on your objective and application.  You can also use the \code{\link{tidy.marssMLE}} function to return a data.frame with the states, standard errors and intervals. \code{\link{tidy.marssMLE}} has much more information on interval calculation for MARSS models.
}
\seealso{ \code{\link{MARSSkf}}, \code{\link{tidy.marssMLE}}}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11,12),]
  MLEobj <- MARSS(dat, model=list(Z=factor(c("WA","OR","OR"))))
  fitted(MLEobj)
  
  # Example of fitting a stochastic level model to the Nile River flow data
  # red line is smooothed level estimate
  # grey line is one-step-ahead prediction using xtt1
  nile <- as.vector(datasets::Nile)
  mod.list = list(
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  x0=matrix("pi") )
  fit <- MARSS(nile, model=mod.list, silent=TRUE)
  plot(nile, type="p", pch=16, col="blue")
  lines(fitted(fit)[1,], col="red", lwd=2)
  lines(fitted(fit, one.step.ahead=TRUE)[1,], col="grey", lwd=2)
  }