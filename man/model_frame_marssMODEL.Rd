\name{model.frame.marssMODEL}
\alias{model.frame}
\alias{model.frame.marssMODEL}
\alias{model.frame.marssMLE}
\keyword{internal}

\title{ model.frame method for marssMLE and marssMODEL objects }
\description{
  \code{model.frame(M]LEobj)} or \code{model.frame(MODELobj)}, where MLEobj is a \code{\link{marssMLE}} object output by a \code{\link{MARSS}()} call and MODELobj is a \code{\link{marssMODEL}} object in the model element of a \code{\link{marssMLE}} object, 
  will return a data frame with the data (y) and inputs/covariates (c and d elements) for a MARSS model in "marxss" form. See \code{\link{MARSS.marxss}}.  This is mainly a utility function to help with the functions (\code{\link{tidy.marssMLE}}, and \code{\link{glance.marssMLE}}).
}
\usage{
\method{model.frame}{marssMODEL}(formula, ...)
}
\arguments{
  \item{formula}{ A \code{\link{marssMODEL}} object.  }
  \item{...}{ Other arguments not used. }
}

\value{
  A data frame with the data and inputs (c and d) in a MARSS model in "marxss" form.  See \code{\link{MARSS.marxss}}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
  
  eli(dot)holmes(at)noaa(dot)gov
}
