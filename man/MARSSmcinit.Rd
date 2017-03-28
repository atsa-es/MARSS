\name{MARSSmcinit}
\alias{MARSSmcinit}
\title{ Monte Carlo Initialization }
\description{
  Performs a Monte Carlo search for optimal initial conditions iterative maximization algorithms (\code{\link{MARSSkem}} and \code{\link{MARSSoptim}}).   This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSmcinit(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
}
\details{
  It is recommended that initial parameter values be set using \code{MARSSmcinit()}, particularly if the model is not a good fit to the data.  This requires more compuation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.

  Options for \code{MARSSmcinit()} may be set using \code{MLEobj$control}, as follows:
  \describe{ 
    \item{\code{MLEobj$control$numInits}}{ Number of random initial value draws. } 
    \item{\code{MLEobj$control$numInitSteps}}{ Maximum number of EM iterations for each random initial value draw. }
    \item{\code{MLEobj$control$boundsInits}}{ Length 6 list. Each component is a length 2 vector of bounds on the uniform distributions from which initial values will be drawn (for A, B, U, and Z). For R and Q, variance-covariance matrices are generated from a wishart distribution with df=bound[1] and S=diag(bound[2],m).  Note, random initial conditions are only used for parameters that are not fixed.} 
  }
  The default values for these are given in \code{MARSSsettings.R} and listed in \code{\link{MARSS}}.
}
\value{
  A list with 8 matrices Z, A, R, B, U, Q, x0, V0, specifying initial values for parameters for iteration 1 of the EM algorithm.  Note the output is the initial values for a marssMODEL in marss form.
}
\references{ 
  The user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.    

 }
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}

\seealso{ \code{\link{MARSSkem}} \code{\link{marssMLE}} \code{\link{MARSS}} }
\examples{
\dontrun{ 
  #Note doing a Monte-Carlo search takes a long, long time
  dat = t(harborSeal)
  dat = dat[c(2,nrow(dat)),]
  fit1=MARSS(dat, control=list(MCInit=TRUE))
  fit1
  #Show the inits that were used
  fit1$start
  #Try fewer initial start locations 
  #and different mean variance (0.1 instead of 1) for R and Q
  cntl.list = list(MCInit=TRUE, numInits=10, 
                   numInitSteps = 10, 
                   boundsInits=list(Q=c(1,0.1),R=c(1,0.1)))
  fit2=MARSS(dat, control=cntl.list)
  fit2
  #Show the inits that were used
  fit2$start
  #ignore the values for Z,B, and V0; those parameters are fixed
}
  }


