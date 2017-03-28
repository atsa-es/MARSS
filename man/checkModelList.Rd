\name{checkModelList}
\alias{checkModelList}
\title{ Check model List Passed into MARSS Call }
\description{
  This is a helper function to check the model list passed in to a MARSS() call for any errors. Not exported.
}
\usage{
checkModelList( model, defaults, this.form.allows)
}
\arguments{
  \item{model}{ A list from which a marssMODEL model will be constructed. }
  \item{defaults}{ A list with the defaults for the elements in the model list in case the user leaves out any. }
  \item{this.form.allows}{ A list of what inputs are allowed for each element in the model list. }
}
\details{
This is a helper function to check that all the model list that will be passed to a \code{MARSS.form} function to make the marssMODEL object.  If elements in the list are left off, they will be filled in by defaults.
}
\value{
  If the function does not stop due to errors, it returns an updated model list with missing elements filled in by the defaults.
}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov   
}
\seealso{ 
  \code{\link{MARSS}} \code{\link{marssMODEL}} \code{\link{checkModelList}} 
}


