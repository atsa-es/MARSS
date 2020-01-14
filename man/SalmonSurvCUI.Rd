\name{SalmonSurvCUI}
\alias{SalmonSurvCUI}
\keyword{internal}
\keyword{datasets}

\docType{data}
\title{ Salmon Survial Indices }
\description{
  Example data set for use in MARSS vignettes for the DLM chapter in the \code{\link{MARSS-package}} User Guide. This is a 42-year time-series of the logit of juvenile salmon survival along with an index of April coastal upwelling.  See the source for details.
}
\usage{
data(SalmonSurvCUI)
}
\format{
The data are provided as a matrix with time running down the rows.  Column 1 is year, column 2 is the logit of the proportion of juveniles that survive to adulthood, column 3 is an index of the April coastal upwelling index.
}
\source{
Scheuerell, Mark D., and John G. Williams. "Forecasting climate-induced changes in the survival of Snake River spring/summer Chinook salmon (Oncorhynchus tshawytscha)." Fisheries Oceanography 14.6 (2005): 448-457.
}
\examples{
str(SalmonSurvCUI)
}

