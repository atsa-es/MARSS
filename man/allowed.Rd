\name{allowed}
\alias{allowed.methods}
\alias{kem.methods}
\alias{optim.methods}
\alias{alldefaults}
\keyword{internal}
\title{ MARSS Function Defaults and Allowed Methods }
\description{
  Defaults and allowed fitting methods for the \code{\link{MARSS}()} function are specified in the file \code{onLoad.R}. These are hidden package globals that are assigned to the package environment when the library is loaded either via `library(MARSS)`, `require(MARSS)` or a call to a MARSS function using `MARSS::`. 
  
\details{
\code{allowed.methods} is a vector with the allowed \code{method} arguments for the \code{\link{MARSS}()} function.  \code{kem.methods} and \code{optim.methods} are vectors of \code{method} arguments that fall in each of these two categories; used by \code{\link{MARSS}()}.  \code{alldefaults} is a list that specifies the defaults for \code{\link{MARSS}()} arguments if they are not passed in.
}


