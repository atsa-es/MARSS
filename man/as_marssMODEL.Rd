\name{marss.conversion}
\alias{marxss_to_marss}
\alias{marss_to_marxss}
\keyword{internal}

\title{ Convert Model Objects between Forms }
\description{
These are utility functions for model objects in the package \code{\link{MARSS-package}}.  
Users would not normally work directly with these functions.
}
\usage{
marss_to_marxss(x, C.and.D.are.zero = FALSE)
marxss_to_marss(x, only.par = FALSE)
}

\arguments{
  \item{x}{ An object of class \code{\link{marssMLE}}. }
  \item{C.and.D.are.zero}{ If the C and D matrices are all 0, then a marss model can be converted to marxss without further information besides the marss model. }
  \item{only.par}{ If only.par=TRUE then only the par element is changed and marss is used for the marss object. }
}
\details{
As the name of the functions imply, these convert \code{\link{marssMODEL}} objects of different forms into other forms.  form=marss is the base form needed for the internal algorithms, thus other (more user friendly forms) must have a \code{form_to_marss} function to convert to the base form.  The printing functions are customized to show output in the user-friendly form, thus a \code{marss_to_form} function is needed for \code{print} and \code{coef} methods for \code{\link{marssMLE}} objects.

}
\value{
  A \code{\link{marssMODEL}} object of the appropriate form.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}

\seealso{ \code{\link{marssMODEL}} } 


