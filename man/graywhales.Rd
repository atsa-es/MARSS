\name{population-count-data}
\alias{population-count-data}
\alias{graywhales}
\alias{grouse}
\alias{prairiechicken}
\alias{wilddogs}
\alias{rockfish}
\alias{okanaganRedds}
\alias{kestrel}
\alias{redstart}
\keyword{datasets}

\docType{data}
\title{ Population Data Sets }
\description{
Example data sets for use in the \code{\link{MARSS-package}} User Guide. Some are logged and some unlogged population counts. See the details below on each data set.

The data sets are matrices with year in the first column and counts in other columns.  Since MARSS functions require time to be across columns, these data matrices must be transposed before passing into MARSS functions.
}
\usage{
data(graywhales)
data(grouse)
data(prairiechicken)
data(wilddogs)
data(kestrel)
data(okanaganRedds)
data(rockfish)
data(redstart)
}
\format{
The data are supplied as a matrix with years in the first column and counts in the second (and higher) columns.
}
\source{
\itemize{
\item{graywhales }{ Gerber L. R., Master D. P. D. and Kareiva P. M. (1999) Gray whales and the value of monitoring data in implementing the U.S. Endangered Species Act.  Conservation Biology, 13, 1215-1219.}
\item{grouse }{ Hays D. W., Tirhi M. J. and Stinson D. W. (1998) Washington state status report for the sharptailed grouse. Washington Department Fish and Wildlife, Olympia, WA. 57 pp.}
\item{prairiechicken }{ Peterson M. J. and Silvy N. J. (1996) Reproductive stages limiting productivity of the endangered Attwater's prairie chicken. Conservation Biology, 10, 1264-1276.}
\item{wilddogs }{ Ginsberg, J. R., Mace, G. M. and Albon, S. (1995). Local extinction in a small and declining population: Wild Dogs in the Serengeti. Proc. R. Soc. Lond. B, 262, 221-228.}
\item{okanaganRedds }{ A data set of Chinook salmon redd (egg nest) surveys.  This data comes from the Okanagan River in Washington state, a major tributary of the Columbia River (headwaters in British Columbia). Unlogged.}
\item{rockfish }{ LOGGED catch per unit effort data for Puget Sound total total rockfish (mix of species) from a series of different types of surveys.}
\item{kestrel }{ Three time series of American kestrel logged abundance from adjacent Canadian provinces along a longitudinal gradient (British Columbia, Alberta, Saskatchewan).  Data have been collected annually, and corrected for changes in observer coverage and detectability. LOGGED.}
\item{redstart}{ 1966 to 1995 counts for American Redstart from the North American Breeding Bird Survey (BBS record number 0214332808636; Peterjohn 1994) used in Dennis et al. (2006).  Peterjohn, B.G. 1994. The North American Breeding Bird Survey. Birding 26, 386--398. and Dennis et al. 2006. Estimating density dependence, process noise, and observation error. Ecological Monographs 76:323-341.}
}
}
\examples{
str(graywhales)
str(grouse)
str(prairiechicken)
str(wilddogs)
str(kestrel)
str(okanaganRedds)
str(rockfish)
}

