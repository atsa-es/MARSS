\name{match.arg.exact}
\alias{match.arg.exact}
\title{ match.arg with exact matching }
\description{
  The base R \code{\link{match.arg}()} uses \code{\link{pmatch}()} and does partial matching. This is a problem for many functions where \code{"xtt1"} is different than \code{"xtt"}, say. This function implements exact matching.
}
\usage{
match.arg.exact(arg, choices, several.ok = FALSE, exact = TRUE) 
}
\arguments{
  \item{arg}{ a character vector (of length one unless \code{several.ok} is \code{TRUE}). }
  \item{choices}{ a character vector of candidate values }
  \item{several.ok}{ logical specifying if \code{arg} should be allowed to have more than one element }
  \item{exact}{ require exact matching }
}
\keyword{internal}
\examples{
# this fails
# MARSS:::match.arg.exact(c("a"), c("aa", "bb"))
# this does not
match.arg(c("a"), c("aa", "bb"))
}

