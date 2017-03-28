\name{print.marssMODEL}
\alias{print.marssMODEL}
\alias{describe.marssMODEL}
\title{ Printing marssMODEL Objects }
\description{
  \code{print(marssMODEL)}., where \code{\link{marssMODEL}} is a marssMODEL object, will print out summary information on the model.
}
\usage{
\method{print}{marssMODEL}(x, ...)
}
\arguments{
  \item{x}{ A marssMODEL object.  }
  \item{...}{ Other arguments for print. }
}

\value{
  A print out of information and invisibly a list with the structure of each parameter matrix.  If you assign the print call to a value, then you can reference the output.  }
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{ 
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)

  print(MLEobj$model)
  #this is identical to 
  print(MLEobj, what="model")
}