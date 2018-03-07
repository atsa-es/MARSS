\name{fitted.marssMLE}
\alias{fitted}
\alias{fitted.marssMLE}
\title{ fitted function for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs marssMLE objects.  \code{fitted(marssMLE)}, where marssMLE is the output from a \code{\link{MARSS}} call, 
  will return the modeled value of y(t) or x(t).  See details.
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

The model predicted (fitted) y at time t is \eqn{Z(t)\hat{x}(t)+a(t)}, where the model is written in marss form. See \code{\link{MARSS.marxss}} for a discussion of the conversion of MARSS models with covariates (c and d) into marss form (which is how models are written in the internal MARSS algorithms).
  
\eqn{\hat{x}(t)} is the expected value of the states at time t. If \code{one.step.ahead=FALSE}, \eqn{\hat{x}(t)} is the expected value conditioned on all the data, i.e. xtT returned by \code{\link{MARSSkf}}.  If \code{one.step.ahead=TRUE}, \eqn{\hat{x}(t)} is the expected value conditioned on the data up to time t-1, i.e. xtt1 returned by \code{\link{MARSSkf}}.  These are commonly known as the one step ahead predictions for a state-space model.

\strong{state fitted values}

The model predicted \eqn{\hat{x}(t)} given \eqn{\hat{x}(t-1)} is \eqn{B(t)\hat{x}(t-1)+u(t)}, where the model is written in marss form. 
  
\eqn{\hat{x}(t)} is the expected value of the states at time t. If \code{one.step.ahead=FALSE}, \eqn{\hat{x}(t)} is the expected value conditioned on all the data, i.e. xtT returned by \code{\link{MARSSkf}}.  If \code{one.step.ahead=TRUE}, \eqn{\hat{x}(t)} is the expected value conditioned on the data up to time t-1, i.e. xtt1 returned by \code{\link{MARSSkf}}.  This type of state fitted value is used for process outlier detection and shock detection.  See \code{\link{residuals.marssMLE}} and read the references cited.  

If you want estimates of the states, which is want one typically wants not the fitted values, then you'll want either the states estimate conditioned on all the data or conditioned on the data up to time t-1.  These are returned by \code{\link{MARSSkf}} in \code{xtT} and \code{xtt1} respectively.  Which one you want depends on your objective and application.  You can also use the \code{\link{tidy.marssMLE}} function to return a data.frame with the xtT, standard errors and confidence intervals.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
