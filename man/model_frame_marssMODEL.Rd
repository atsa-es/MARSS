\name{model.frame.marssMODEL}
\alias{model.frame}
\alias{model.frame.marssMODEL}
\alias{model.frame.marssMLE}
\keyword{internal}

\title{ model.frame method for marssMLE and marssMODEL objects }
\description{
  \code{model.frame(marssMLE)} or \code{model.frame(marssMODEL)}, where marssMLE is the output from a \code{\link{MARSS}} call marssMODEL is the model element of a marssMLE object, 
  will return a data.frame with the data (y) and inputs/covariates (c and d elements) for a MARXSS model. See \code{\link{MARSS.marxss}}.  This is mainly a utility function to help with the \code{broom} functions (\code{tidy}, \code{augment} and \code{glance}).
}
\usage{
  \method{model.frame}{marssMODEL}(formula, ...)
}
\arguments{
  \item{formula}{ A marssMODEL object.  }
  \item{...}{ Other arguments not used. }
}

\value{
  A data.frame witht the data and inputs (c and d) in a MARXSS model.  See \code{\link{MARSS.marxss}}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
  
  eli(dot)holmes(at)noaa(dot)gov
}
