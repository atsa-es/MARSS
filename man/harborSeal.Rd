\name{harborSeal}
\alias{harborSeal}
\alias{harborSealWA}
\keyword{datasets}

\docType{data}
\title{ Harbor Seal Population Count Data (Log counts) }
\description{
Data sets used in MARSS vignettes in the \code{\link{MARSS-package}}.  These are data sets based on LOGGED count data from Oregon, Washington and California sites where harbor seals were censused while hauled out on land. "harborSealnomiss" is an extrapolated data set where missing values in the original data set have been extrapolated so that the data set can be used to demonstrate fitting population models with different underlying structures.
}
\usage{
data(harborSeal)
data(harborSealWA)
}
\format{
  Matrix "harborSeal" contains columns "Years", "StraitJuanDeFuca", "SanJuanIslands", "EasternBays", "PugetSound", "HoodCanal", "CoastalEstuaries", "OlympicPeninsula", "CA.Mainland",  "OR.NorthCoast", "CA.ChannelIslands", and "OR.SouthCoast".
    Matrix "harborSealnomiss" contains columns "Years", "StraitJuanDeFuca", "SanJuanIslands", "EasternBays", "PugetSound", "HoodCanal", "CoastalEstuaries", "OlympicPeninsula", "OR.NorthCoast", and "OR.SouthCoast".
 Matrix "harborSealWA" contains columns "Years", "SJF", "SJI", "EBays", "PSnd", and "HC", representing the same five sites as the first five columns of "harborSeal".
}
\details{
Matrix "harborSealWA" contains the original 1978-1999 LOGGED count data for five inland WA sites. Matrix "harborSealnomiss" contains 1975-2003 data for the same sites as well as four coastal sites, where missing values have been replaced with extrapolated values.  Matrix "harborSeal" contains the original 1975-2003 LOGGED data (with missing values) for the WA and OR sites as well as a CA Mainland and CA ChannelIslands time series.
}
\source{
  Jeffries et al. 2003.  Trends and status of harbor seals in Washington State: 1978-1999. Journal of Wildlife Management 67(1):208-219.
  
  Brown, R. F., Wright, B. E., Riemer, S. D. and Laake, J. 2005. Trends in abundance and current status of harbor seals in Oregon: 1977-2003. Marine Mammal Science, 21: 657-670.
  
  Lowry, M. S., Carretta, J. V., and Forney, K. A. 2008.  Pacific harbor seal census in California during May-July 2002 and 2004. California Fish and Game 94:180-193.
  
  Hanan, D. A. 1996. Dynamics of Abundance and Distribution for Pacific Harbor Seal, Phoca vitulina richardsi, on the Coast of California. Ph.D. Dissertation, University of California, Los Angeles. 158pp. DFO. 2010. Population Assessment Pacific Harbour Seal (Phoca vitulina richardsi). DFO Can. Sci. Advis. Sec. Sci. Advis. Rep. 2009/011.
  }
\examples{
str(harborSealWA)
str(harborSeal)
}
\keyword{datasets}
