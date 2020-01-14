\name{loggerhead}
\alias{loggerhead}
\alias{loggerheadNoisy}
\keyword{internal}
\keyword{datasets}

\docType{data}
\title{ Loggerhead Turtle Tracking Data }
\description{
Data used in MARSS vignettes in the \code{\link{MARSS-package}}.  Tracking data from ARGOS tags on eight individual loggerhead turtles, 1997-2006.
}
\usage{
data(loggerhead)
data(loggerheadNoisy)
}
\format{
  Data frames "loggerhead" and "loggerheadNoisy" contain the following columns:
  \describe{
    \item{turtle}{ Turtle name. }
    \item{day}{ Day of the month (character). }
    \item{month}{ Month number (character). }
    \item{year}{ Year (character). }
    \item{lon}{ Longitude of observation. }
    \item{lat}{ Latitude of observation. } 
  }
}
\details{
  Data frame "loggerhead" contains the original latitude and longitude data. Data frame "loggerheadNoisy" has noise added to the lat and lon data to represent data corrupted by errors.
}
\source{
  Gray's Reef National Marine Sanctuary (Georgia) and WhaleNet:
  http://whale.wheelock.edu/whalenet-stuff/stop_cover_archive.html
}
\examples{
str(loggerhead)
str(loggerheadNoisy)
}
