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

  A marssMLE object in the  \code{\link{MARSS-package}} that has all the elements needed for maximum-likelihood estimation of multivariate autoregressive state-space model: the data, model, estimation methods, and any control options needed for the method.  If the model has been fit and parameters estimated, the object will also have the MLE parameters.  Other functions add other elements to the marssMLE object, such as CIs, s.e.'s, AICs, and the observed Fisher Information matrix. There are \code{print}, \code{summary}, \code{coef}, \code{fitted}, \code{residuals}, \code{predict} and \code{simulate} methods for \code{\link{marssMLE}} objects and a bootstrap function.  Rather than working directly with the elements of a \code{\link{marssMLE}} object, use \code{\link[=print.marssMLE]{print}()}, \code{\link[=tidy.marssMLE]{tidy}()}, \code{\link[=fitted.marssMLE]{fitted}()}, \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}, \code{\link[=predict.marssMLE]{predict}()}, or \code{\link[=residuals.marssMLE]{residuals}()} to extract output.
  }

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "marssMLE")}: ... }
    \item{summary}{\code{signature(object = "marssMLE")}: ... }
    \item{coef}{\code{signature(object = "marssMLE")}: ... }
    \item{residuals}{\code{signature(object = "marssMLE")}: ... }
    \item{predict}{\code{signature(object = "marssMLE")}: ... }
    \item{fitted}{\code{signature(object = "marssMLE")}: ... }
    \item{logLik}{\code{signature(object = "marssMLE")}: ... }
    \item{simulate}{\code{signature(object = "marssMLE")}: ... }
    \item{forecast}{\code{signature(object = "marssMLE")}: ... }
    \item{accuracy}{\code{signature(object = "marssMLE")}: ... }
    \item{toLatex}{\code{signature(object = "marssMLE")}: ... }
	 }
}
\seealso{
\code{\link{is.marssMLE}()}, \code{\link{print.marssMLE}()}, \code{\link{summary.marssMLE}()}, \code{\link{coef.marssMLE}()}, \code{\link{residuals.marssMLE}()}, \code{\link{fitted.marssMLE}()}, \code{\link{tsSmooth.marssMLE}()}, \code{\link{logLik.marssMLE}()}, \code{\link{simulate.marssMLE}()}, \code{\link{predict.marssMLE}()}, \code{\link{forecast.marssMLE}()}, \code{\link{accuracy.marssMLE}()}, \code{\link{toLatex.marssMLE}()}
}
\author{ 
  Eli Holmes and Kellie Wills, NOAA, Seattle, USA
}
\keyword{classes}
