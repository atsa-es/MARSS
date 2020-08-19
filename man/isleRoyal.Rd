\name{isleRoyal}
\alias{isleRoyal}
\keyword{datasets}

\docType{data}
\title{ Isle Royale Wolf and Moose Data }
\description{
Example data set for estimation of species interaction strengths. These are data on the number of wolves and 
  moose on Isle Royal, Michigan.  The data are unlogged.  The covariate data are the average Jan-Feb, average Apr-May and average July-Sept
  temperature (Fahrenheit) and precipitation (inches).  Also included are 3-year running means of these covariates.  The choice of 
  covariates is based on those presented in the Isle Royale 2012 annual report.
}
\usage{
data(isleRoyal)
}
\format{
The data are supplied as a matrix with years in the first column.
}
\source{
Peterson R. O., Thomas N. J., Thurber J. M., Vucetich J. A. and Waite T. A. (1998) Population limitation and the wolves of Isle Royale. In: Biology and Conservation of Wild Canids (eds. D. Macdonald and C. Sillero-Zubiri). Oxford University Press, Oxford, pp. 281-292.

Vucetich, J. A. and R. O. Peterson. (2010) Ecological studies of wolves on Isle Royale. Annual Report 2009-10. School of Forest Resources and Environmental Science, Michigan Technological University, Houghton, Michigan USA 49931-1295

The source for the covariate data is the Western Regional Climate Center (http://www.wrcc.dri.edu) using their data for the NE Minnesota division.  The website used was http://www.wrcc.dri.edu/cgi-bin/divplot1_form.pl?2103 and www.wrcc.dri.edu/spi/divplot1map.html.
}
\examples{
str(isleRoyal)
}
