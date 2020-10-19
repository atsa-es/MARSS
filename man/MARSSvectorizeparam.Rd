\name{MARSSvectorizeparam}
\alias{MARSSvectorizeparam}
\keyword{internal}

\title{ Vectorize or Replace the par List }
\description{
  Converts \code{MLEobj[["what"]]} to a vector or assigns a vector to \code{MLEobj[["what"]]}.  This is a utility function in the \code{\link{MARSS-package}} for \code{\link{marssMODEL}} objects of form="marss" and is not exported.  Users achieve this functionality with \code{\link[=coef.marssMLE]{coef}}.
}
\usage{
MARSSvectorizeparam(MLEobj, parvec = NA, what="par")
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
  \item{parvec}{ NA or a vector. See Value. }
  \item{what}{ What part of the MLEobj is being replaced or vectorized.  Need to be a par list. }
}
\details{
  Utility function to generate parameter vectors for optimization functions, and to set \code{MLEobj[[what]]} using a vector of values.  The function bases the unlisting and naming order on \code{names(MLEobj$marss$fixed)}.  Appends matrix name to the row names in the par list.
}
\value{
  If parvec=NA, a vector of the elements of the \code{what} element. Otherwise, a \code{\link{marssMLE}} object with \code{MLEobj[["what"]]} set by parvec.
}
\author{ 
  Eli Holmes and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ \code{\link{marssMLE}} }
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4, ]
kem <- MARSS(dat)
paramvec <- MARSS:::MARSSvectorizeparam(kem)
paramvec
}

