\name{print.marssMODEL}
\alias{print.marssMODEL}
\alias{summary.marssMODEL}
\alias{describe.marssMODEL}
\title{ Printing marssMODEL Objects }
\description{
  \code{print(marssMODEL)}, where \code{\link{marssMODEL}} is a marssMODEL object, will print out information on the model in short form (e.g. 'diagonal and equal').
  
  \code{summary(marssMODEL)}, where \code{\link{marssMODEL}} is a marssMODEL object, will print out detailed information on each parameter matrix showing where the estimated values (and their names) occur.
}
\usage{
\method{print}{marssMODEL}(x, ...)
\method{summary}{marssMODEL}(object, ..., silent = FALSE)
}
\arguments{
  \item{x}{ A marssMODEL object.  }
  \item{object}{ A marssMODEL object.  }
  \item{...}{ Other arguments . }
  \item{silent}{ TRUE/FALSE Whether to print output. }
}

\value{
  \code{print(marssMODEL)} prints out of the structure of each parameter matrix in 'English' (e.g. 'diagonal and unequal') and returns invisibly the list.  If you assign the print call to a value, then you can reference the output. 

  \code{summary(marssMODEL)} prints out of the structure of each parameter matrix in as list matrices showing where each estimated value occurs in each matrix and returns invisibly the list.  The output can be verbose, especially if parameter matrices are time-varying.  Pass in \code{silent=TRUE} and assign output (a list with each parameter matrix) to a variable.  Then specific parameters can be looked at.

  }
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{ 
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  MLEobj <- MARSS(dat)

  print(MLEobj$model)
  #this is identical to 
  print(MLEobj, what="model")
}