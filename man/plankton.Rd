\name{xplankton}
\alias{plankton}
\alias{ivesDataLP}
\alias{ivesDataByWeek}
\alias{lakeWAplankton}
\alias{lakeWAplanktonRaw}
\alias{lakeWAplanktonTrans}
%\keyword{internal}
\keyword{datasets}

\docType{data}
\title{ Plankton Data Sets }
\description{
  Example data sets for use in MARSS vignettes for the \code{\link{MARSS-package}}. 
  
The lakeWAplankton dataset consists for two datasets: lakeWAplanktonRaw and a dataset derived from the raw dataset, lakeWAplanktonTrans. lakeWAplanktonRaw is a 32-year time series (1962-1994) of monthly plankton counts from Lake Washington, Washington, USA. Columns 1 and 2 are year and month.  Column 3 is temperature (C), column 4 is total phosphorous, and column 5 is pH.  The next columns are the plankton counts in units of cells per mL for the phytoplankton and organisms per L for the zooplankton.   Since MARSS functions require time to be across columns, these data matrices must be transposed before passing into MARSS functions.

lakeWAplanktonTrans is a transformed version of lakeWAplanktonRaw.  Zeros have been replaced with NAs (missing).  The logged (natural log) raw plankton counts have been standardized to a mean of zero and variance of 1 (so logged and then z-scored). Temperature, TP & pH were also z-scored but not logged (so z-score of the untransformed values for these covariates). The single missing temperature value was replaced with -1 and the single missing TP value was replaced with -0.3.
  
  The Ives data are from Ives et al. (2003) for West Long Lake (the low planktivory case).  The Ives data are unlogged.  \code{ivesDataLP} and \code{ivesDataByWeek} are the same data with LP having the missing weeks in winter removed while in ByWeek, the missing values are left in.  The phosporous column is the experimental input rate + the natural input rate for phosphorous, and Ives et al. used 0.1 for the natural input rate when no extra phosporous was added.  The phosporous input rates for weeks with no sampling (and no experimental phosporous input) have been filled with 0.1 in the "by week" data.
}
\usage{
data(ivesDataLP)
data(ivesDataByWeek)
data(lakeWAplankton)
}
\format{
The data are provided as a matrix with time running down the rows.
}
\source{
\itemize{
\item{ivesDataLP and ivesDataByWeek}{ Ives, A. R. Dennis, B. Cottingham, K. L. Carpenter, S. R. (2003) Estimating community stability and ecological interactions from time-series data. Ecological Monographs, 73, 301-330. }
\item{lakeWAplanktonTrans}{ Hampton, S. E. Scheuerell, M. D. Schindler, D. E. (2006)  Coalescence in the Lake Washington story: Interaction strengths in a planktonic food web. Limnology and Oceanography, 51, 2042-2051. }
\item{lakeWAplanktonRaw}{ Adapted from the Lake Washington database of Dr. W. T. Edmondson, as funded by the Andrew Mellon Foundation; data courtesy of Dr. Daniel Schindler, University of Washington, Seattle, WA.}
}
}
\examples{
str(ivesDataLP)
str(ivesDataByWeek)
}

