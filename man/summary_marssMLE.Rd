\name{summary.marssMLE}
\alias{summary.marssMLE}

\title{ Summary methods for marssMLE objects }
\description{
A brief summary of the fit: number of state and observation time series and the estimates. See also \code{\link{glance.marssMLE}()} and \code{\link{tidy.marssMLE}()} for other summary like output.
}
\usage{
\method{summary}{marssMLE}(object, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{digits}{ Number of digits for printing.  }
  \item{...}{ Not used. }
}
\value{
  Returns `object` invisibly.
}
\seealso{ \code{\link{glance.marssMLE}()}, \code{\link{tidy.marssMLE}()} } 
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{ 
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  fit <- MARSS(dat)
  
  summary(fit)
  glance(fit)
  tidy(fit)
  }