\name{CSEGriskfigure}
\alias{CSEGriskfigure}
\keyword{experimental}
\title{ Plot Extinction Risk Metrics }
\description{
  Generates a six-panel plot of extinction risk metrics used in Population Viability Analysis (PVA). This is a function used by one of the vignettes in the \code{\link{MARSS-package}}. 
}
\usage{
CSEGriskfigure(data, te = 100, absolutethresh = FALSE, threshold = 0.1, 
  datalogged = FALSE, silent = FALSE, return.model = FALSE, 
  CI.method = "hessian", CI.sim = 1000)
}
\arguments{
  \item{data}{ A data matrix with 2 columns; time in first column and counts in second column. Note time is down rows, which is different than the base \code{\link{MARSS-package}} functions.}
  \item{te}{ Length of forecast period (positive integer) }
  \item{absolutethresh}{ Is extinction threshold an absolute number? (T/F)}
  \item{threshold}{ Extinction threshold either as an absolute number, if \code{absolutethresh=TRUE}, or as a fraction of current population count, if \code{absolutethresh=FALSE}. }
  \item{datalogged}{ Are the data already logged? (T/F)}
  \item{silent}{ Suppress printed output? (T/F) }
  \item{return.model}{ Return state-space model as \code{\link{marssMLE}} object? (T/F)}
  \item{CI.method}{ Confidence interval method: "hessian", "parametrc", "innovations", or "none". See \code{\link{MARSSparamCIs}}. }
  \item{CI.sim}{ Number of simulations for bootstrap confidence intervals (positive integer). }
}
\details{
  Panel 1: Time-series plot of the data.
  Panel 2: CDF of extinction risk.
  Panel 3: PDF of time to reach threshold.
  Panel 4: Probability of reaching different thresholds during forecast period.
  Panel 5: Sample projections.
  Panel 6: TMU plot (uncertainty as a function of the forecast).
}
\value{
  If \code{return.model=TRUE}, an object of class \code{\link{marssMLE}}.
}
\references{ 
Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.

(theory behind the figure) Holmes, E. E., J. L. Sabo, S. V. Viscido, and W. F. Fagan. (2007) A statistical approach to quasi-extinction forecasting. Ecology Letters 10:1182-1198.

(CDF and PDF calculations) Dennis, B., P. L. Munholland, and J. M. Scott. (1991) Estimation of growth and extinction parameters for endangered species. Ecological Monographs 61:115-143.

(TMU figure) Ellner, S. P. and E. E. Holmes. (2008) Resolving the debate on when extinction risk is predictable. Ecology Letters 11:E1-E5.
 }
\author{ 
  Eli Holmes, NOAA, Seattle, USA, and Steve Ellner, Cornell Univ. 
}
\seealso{ \code{\link{MARSSboot}}, \code{\link{marssMLE}}, \code{\link{CSEGtmufigure}} }
\examples{
d <- harborSeal[, 1:2]
kem <- CSEGriskfigure(d, datalogged = TRUE)
}
\keyword{ hplot }
