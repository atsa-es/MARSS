\name{MARSSapplynames}
\alias{MARSSapplynames}
\keyword{internal}

\title{ Names for marssMLE Object Components }
\description{
Puts names on the par, start, par.se, init components of \code{\link{marssMLE}} objects.  This is a utility function in the \code{\link{MARSS-package}} and is not exported.
}
\usage{
MARSSapplynames(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
}
\details{
The X.names and Y.names are attributes of \code{\link{marssMODEL}} objects (which would be in \code{$marss} and \code{$model} in the \code{\link{marssMLE}} object).  These names are applied to the par elements in the \code{\link{marssMLE}} object.
}
\value{
  The object passed in, with row and column names on matrices as specified.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\seealso{ \code{\link{marssMLE}}, \code{\link{marssMODEL}} }


