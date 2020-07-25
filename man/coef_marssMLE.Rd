\name{coef.marssMLE}
\alias{coef.marssMLE}

\title{ Coefficient function for MARSS MLE objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{coef(MLEobj)}, where MLEobj is one's output from a \code{\link{MARSS}()} call, will print out the estimated parameters.  The default output is a list with values for each parameter, however the output can be altered using the \code{type} argument to output a vector of all the estimated values (\code{type="vector"}) or a list with the full parameter matrix with the estimated and fixed elements (\code{type="matrix"}).
}
\usage{
\method{coef}{marssMLE}(object, ..., type="list", form=NULL, what="par")
}
\arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{...}{ Other arguments for coef. }
  \item{type}{ What to print.  Default is "list".  If you input type as a vector, coef returns a list of output. See examples.
  \itemize{
    \item{ "list" }{ A list of only the estimated values in each matrix. Each model matrix has it's own list element.}
    \item{ "vector" }{ A vector of all the estimated values in each matrix. }
    \item{ "matrix" }{ A list of the parameter matrices each parameter with fixed values at their fixed values and the estimated values at their estimated values. Time-varying parameters, including d and c in a marxss form model, are returned as an array with time in the 3rd dimension. }
    \item{ parameter name }{ Returns the parameter matrix for that parameter with fixed values at their fixed values and the estimated values at their estimated values. Note, time-varying parameters, including d and c in a marxss form model, are returned as an array with time in the 3rd dimension.}
  } }
  \item{form}{ By default, the model form specified in the call to \code{\link{MARSS}()} is used to determine how to display the coefficients.  This information is in \code{ attr(object$model,"form") }.  The default form is \code{"marxss"}; see \code{\link{MARSS.marxss}}.  However, the internal functions convert this to form \code{"marss"}; see \code{\link{MARSS.marss}}. The marss form of the model is stored (in \code{object$marss}). You can look at the coefficients in marss form by passing in \code{form="marss"}. This is mainly useful is for debugging numerical problems since the error reports will be for the "marss" form.}
  \item{what}{ By default, coef shows the parameter estimates.  Other options are "par.se", "par.lowCI", "par.upCI", "par.bias", and "start".}
}

\value{
  A list of the estimated parameters for each model matrix.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{
\code{\link{tidy.marssMLE}}, \code{\link{print.marssMLE}} 
}
\examples{ 
dat <- t(harborSeal)
dat <- dat[c(2, 11), ]
fit <- MARSS(dat)

coef(fit)
coef(fit, type = "vector")
coef(fit, type = "matrix")
# to retrieve just the Q matrix
coef(fit, type = "matrix")$Q
}