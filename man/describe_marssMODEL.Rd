\name{describe.marssMODEL}
\alias{describe.marssMODEL}
\keyword{internal}

\title{ Describe a marssMODEL Objects }
\description{
  \code{describe.marssMODEL} will print out information on the model in short form (e.g. 'diagonal and equal'). It is used by \code{\link{print.marssMODEL}}. It calls form specific functions: \code{describe_dfa}, \code{describe_marss}, and \code{describe_marxss}.
}
\usage{
describe.marssMODEL(x)
}
\arguments{
  \item{x}{ A marssMODEL object.  }
}

\value{
  \code{describe.marssMODEL(marssMODEL)} returns a list with the structure of each parameter matrix in 'English' (e.g. 'diagonal and unequal'). 
  }
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{ 
dat <- t(harborSeal)
dat <- dat[c(2, 11), ]
MLEobj <- MARSS(dat)

MARSS:::describe.marssMODEL(MLEobj$model)
}