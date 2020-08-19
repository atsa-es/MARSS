\name{marssResiduals-class}
\Rdversion{1.1}
\docType{class}
\alias{marssResiduals}

\title{Class "marssResiduals"}
\description{
\code{\link{marssResiduals}} are the objects returned by \code{\link{residuals.marssMLE}} in the package \code{\link{MARSS-package}}. It is a data frame in tibble format (but not tibble class).

\code{standardization} 
\itemize{
\item "Cholesky" means it is standardized by the Cholesky transformation of the full variance-covariance matrix of the model and state residuals. 
\item "marginal" means that the residual is standardized by its standard deviation, i.e. the square root of the value on the diagonal of the variance-covariance matrix of the model and state residuals. 
}

\code{type} 
\itemize{
\item "tT" means the fitted values are computed conditioned on all the data. See \code{\link[=fitted.marssMLE]{fitted}} with \code{type="ytT"} or \code{type="xtT"}.
\item "tt1" means the fitted values are computed conditioned on the data up to t-1. See \code{\link[=fitted.marssMLE]{fitted}} with \code{type="ytt1"} or \code{type="xtt1"}.

}

}

\seealso{
\code{\link{residuals.marssMLE}}, \code{\link{MARSSresiduals}} 
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA
  
  eli(dot)holmes(at)noaa(dot)gov
}
\keyword{classes}
