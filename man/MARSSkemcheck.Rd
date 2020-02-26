\name{MARSSkemcheck}
\alias{MARSSkemcheck}
\keyword{internal}

\title{ Model Checking for MLE objects Passed to MARSSkem }
\description{
  This is a helper function in the \code{\link{MARSS-package}} that checks that the model can be handled by the \code{\link{MARSSkem}} algorithm. It also returns the structure of the model as a list of text strings.
}
\usage{
MARSSkemcheck(MLEobj)
}
\arguments{
  \item{ MLEobj }{ An object of class \code{\link{marssMLE}}. }
}
\value{
  A list with of the model elements A, B, Q, R, U, x0, Z, V0 specifying the structure of the model using text strings).
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{marssMODEL}}, \code{\link{MARSSkem}} }

