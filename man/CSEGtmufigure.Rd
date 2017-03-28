\name{CSEGtmufigure}
\alias{CSEGtmufigure}
\title{ Plot Forecast Uncertainty }
\description{
  Plot the uncertainty in the probability of hitting a percent threshold (quasi-extinction) for a single random walk trajectory.  This is the quasi-extinction probability used in Population Viability Analysis.  The uncertainty is shown as a function of the forecast, where the forecast is defined in terms of the forecast length (number of time steps) and forecasted decline (percentage).  This is a function used by one of the vignettes in the \code{\link{MARSS-package}}. 
}
\usage{
CSEGtmufigure(N = 20, u = -0.1, s2p = 0.01, make.legend = TRUE)
}
\arguments{
  \item{N}{ Time steps between the first and last population data point (positive integer) }
  \item{u}{ Per time-step decline (-0.1 means a 10\% decline per time step; 1 means a doubling per time step.) }
  \item{s2p}{ Process variance (Q). (a positive number)}
  \item{make.legend}{ Add a legend to the plot? (T/F) }
}
\details{
  This figure shows the region of high uncertainty in dark grey. In this region, the minimum 95 percent confidence intervals on the probability of quasi-extinction span 80 percent of the 0 to 1 probability. Green hashing indicates where the 95 percent upper bound does not exceed 5\% probability of quasi-extinction.  The red hashing indicates, where the 95 percent lower bound is above 95\% probability of quasi-extinction. The light grey lies between these two certain/uncertain extremes.   The extinction calculation is based on Dennis et al. (1991). The minimum theoretical confidence interval is based on Fieberg and Ellner (2000).  This figure was developed in Ellner and Holmes (2008).
  
 Examples using this figure are shown in the user guide (\code{RShowDoc("UserGuide",package="MARSS")}) in the PVA case study.
}

\references{ 
Dennis, B., P. L. Munholland, and J. M. Scott. (1991) Estimation of growth and extinction parameters for endangered species. Ecological Monographs 61:115-143.

Fieberg, J. and Ellner, S.P. (2000) When is it meaningful to estimate an extinction probability? Ecology, 81, 2040-2047.

Ellner, S. P. and E. E. Holmes. (2008) Resolving the debate on when extinction risk is predictable. Ecology Letters 11:E1-E5.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA, and Steve Ellner, Cornell Univ.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{CSEGriskfigure}} }
\examples{
CSEGtmufigure(N = 20, u = -0.1, s2p = 0.01)
}
\keyword{ hplot }

