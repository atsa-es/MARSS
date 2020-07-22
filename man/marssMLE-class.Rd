\name{marssMLE-class}
\Rdversion{1.1}
\docType{class}
\alias{marssMLE}
%\alias{print,marssMLE-method}
%\alias{summary,marssMLE-method}
%\alias{coef,marssMLE-method}
%\alias{residuals,marssMLE-method}
%\alias{predict,marssMLE-method}
%\alias{simulate,marssMLE-method}

\title{Class "marssMLE"}
\description{
\code{\link{marssMLE}} objects specify fitted multivariate autoregressive state-space models (maximum-likelihood) in the package \code{\link{MARSS-package}}.

  A marssMLE object in the  \code{\link{MARSS-package}} that has all the elements needed for maximum-likelihood estimation of multivariate autoregressive state-space model: the data, model, estimation methods, and any control options needed for the method.  If the model has been fit and parameters estimated, the object will also have the MLE parameters.  Other functions add other elements to the marssMLE object, such as CIs, s.e.'s, AICs, and the observed Fisher Information matrix. There are print, summary, coef, residuals, predict and simulate methods for \code{\link{marssMLE}} objects and a bootstrap function.  Rather than working directly with the elements of a \code{\link{marssMLE}} object, use \code{\link{print.marssMLE}}, \code{\link{tidy.marssMLE}}, \code{\link{fitted.marssMLE}}, or \code{\link{residuals.marssMLE}} to extract output.
  }

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "marssMLE")}: ... }
    \item{summary}{\code{signature(object = "marssMLE")}: ... }
    \item{coef}{\code{signature(object = "marssMLE")}: ... }
    \item{predict}{\code{signature(object = "marssMLE")}: ... }
    \item{simulate}{\code{signature(object = "marssMLE")}: ... }
	 }
}
\seealso{
\code{\link{is.marssMLE}}
}
\author{ 
  Eli Holmes and Kellie Wills, NOAA, Seattle, USA
  
  eli(dot)holmes(at)noaa(dot)gov
}
\keyword{classes}
