\name{logLik.marssMLE}
\alias{logLik.marssMLE}
\alias{logLik}

\title{ logLik method for MARSS MLE objects }
\description{
  Returns a logLik class object with attributes nobs and df.
}
  \usage{
  \method{logLik}{marssMLE}(object, ...)
  }
  \arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{...}{ Other arguments. Not used. }
  }
\value{
 An object of class logLik.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\seealso{ \code{\link{MARSSkf}} }
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
MLEobj <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))
logLik(MLEobj)

  }