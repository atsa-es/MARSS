\name{MARSSsimulate}
\alias{MARSSsimulate}
\alias{simulate.marssMLE}
\title{ Simulate Data from a MARSS Model }
\description{
  Generates simulated data from a MARSS model with specified parameter estimates. This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSsimulate(object, tSteps = NULL, nsim = 1, silent = TRUE, 
  miss.loc = NULL)
}
\arguments{
  \item{object}{ A fitted \code{\link{marssMLE}} object, as output by \code{\link{MARSS}()}. }
  \item{tSteps}{ Number of time steps in each simulation. If left off, it is taken to be consistent with \code{MLEobj}.}
  \item{nsim}{ Number of simulated data sets to generate. }
  \item{silent}{ Suppresses progress bar. }
  \item{miss.loc}{ Optional matrix specifying where to put missing values. See Details. }
}
\details{
  Optional argument \code{miss.loc} is an array of dimensions n x tSteps x nsim, specifying where to put missing values
  in the simulated data. If missing, this would be constructed using \code{MLEobj$marss$data}.  If the locations of the missing values are the same for all simulations, \code{miss.loc} can be a matrix of dim=c(n, tSteps) (the original data for example).  The default, if \code{miss.loc} is left off, is that there are no missing values even if \code{MLEobj$marss$data} has missing values. 
}
\value{
  \item{sim.states}{ Array (dim m x tSteps x nsim) of state processes simulated from parameter estimates. m is the number of states (rows in X).}
  \item{sim.data}{ Array (dim n x tSteps x nsim) of data simulated from parameter estimates. n is the number of rows of data (Y).}
  \item{MLEobj}{ The \code{\link{marssMLE}} object from which the data were simulated. }
  \item{miss.loc}{ Matrix identifying where missing values were placed. It should be exactly the same dimensions as the data matrix. The location of NAs in the miss.loc matrix indicate where the missing values are. }
  \item{tSteps}{ Number of time steps in each simulation. }
  \item{nsim}{ Number of simulated data sets generated. }
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{marssMODEL}}, \code{\link{marssMLE}}, \code{\link{MARSSboot}} 
}
\examples{
d <- harborSeal[, c(2, 11)]
dat <- t(d)
MLEobj <- MARSS(dat)

# simulate data that are the
# same length as original data and no missing data
sim.obj <- MARSSsimulate(MLEobj, tSteps = dim(d)[1], nsim = 5)

# simulate data that are the
# same length as original data and have missing data in the same location
sim.obj <- MARSSsimulate(MLEobj, tSteps = dim(d)[1], nsim = 5, miss.loc = dat)
}

