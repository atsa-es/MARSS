\name{fitted.marssMLE}
\alias{fitted}
\alias{fitted.marssMLE}
\title{ fitted function for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs marssMLE objects.  \code{fitted(marssMLE)}, where marssMLE is the output from a \code{\link{MARSS}} call, 
  will print out the fitted values for y conditioned on all the data or the data up to t-1.  If the former and y(t) is not missing,
  then the fitted value is simply y(t).  If y(t) is missing, the value is Z(t)x(t|y_1^T)+a(t) conditioned on any observed y at time t.
  The latter is important if R is not diagonal. If the latter, fitted is always  Z(t)x(t|y_1^{t-1})+a(t)).
}
  \usage{
  \method{fitted}{marssMLE}(object, ..., one.step.ahead=FALSE)
  }
  \arguments{
  \item{object}{ A marssMLE object.  }
  \item{...}{ Other arguments for fitted. }
  \item{one.step.ahead}{ If true, the fitted value is computed conditioned on the data up to  }
  }

\value{
  A n x T matrix of the fitted values. If one.step.ahead=FALSE, fitted = Z(t) kf$xtT + a(t).  
  If one.step.ahead=TRUE, fitted = Z(t) kf$xtt1 + a(t).  \code{kf} is the output from \code{\link{MARSSkf}}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
