\name{marssPredict-class}
\Rdversion{1.1}
\docType{class}
\alias{marssPredict}

\title{Class "marssPredict"}
\description{
\code{marssPredict} objects are returned by \code{\link{predict.marssMLE}} and \code{\link{forecast.marssMLE}}.

  A marssPredict object in the  \code{\link{MARSS-package}} has the output with intervals, the original model and values needed for plotting. The object is mainly used for \code{\link{plot.marssPredict}} and \code{\link{print.marssPredict}}.
  }

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "marssPredict")}: ... }
    \item{plot}{\code{signature(object = "marssPredict")}: ... }
	 }
}
\seealso{
\code{\link{plot.marssPredict}}, \code{\link{predict.marssMLE}}, \code{\link{forecast.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, WA.
}
\keyword{classes}
