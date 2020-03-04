\name{accuracy.marssMLE}
\alias{accuracy.marssMLE}
\alias{accuracy.marssPredict}

\title{Return accuracy metrics}
\usage{
accuracy.marssPredict(f, x, test=NULL, type="ytt1", verbose=FALSE, ...)
accuracy.marssMLE(f, x, test=NULL, type="ytt1", verbose=FALSE, ...)
}
\arguments{
  \item{f}{A \code{\link{marssMLE}} or \code{\link{marssPredict}} object}
  \item{x}{A matrix or tibble with data to test against the h steps of a forecast.}
  \item{test}{Which time steps in training data (data model fit to) to compute accuracy for.}
  \item{type}{type="ytt1" is the one-step-ahead predictions. type="ytT" is the fitted ytT predictions. The former are standardly used for training data prediction metrics.}  
  \item{verbose}{Show metrics for each time series of data.}
  \item{...}{Not used.}
}
\description{
 This is a method for the generic \code{accuracy} function in the forecast package. The forecast package is not part of the MARSS package and you will need to load it separately (See examples). Alternatively you can use \code{accuracy.marssMLE()} or \code{accuracy.marssPredict()}.
 
The measures calculated are:
\itemize{
  \item ME: Mean Error
  \item RMSE: Root Mean Squared Error
  \item MAE: Mean Absolute Error
  \item MPE: Mean Percentage Error
  \item MAPE: Mean Absolute Percentage Error
  \item MASE: Mean Absolute Scaled Error
  \item ACF1: Autocorrelation of errors at lag 1.
}
The MASE calculation is scaled using MAE of training set naive
forecasts which is simply y(i) at t-1.

For the training data, the metrics are shown for the one-step-ahead predictions by default (\code{type="ytt1"}). This is the prediction of y(t) conditioned on the data up to t-1 (and the model estimated from all the data). With type="ytT", you can compute the metrics for the fitted ytT, which is the expected value of new data at t conditioned on all the data. \code{type} does not affect test data (forecasts past the end of the training data).
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12),]
fit<- MARSS(dat[,1:12], model = list(Z = factor(c("WA", "OR", "OR"))))

library(forecast)
accuracy(fit)

fr <- predict(fit, h=10)
accuracy(fit, x=dat[,13:22])
}
\references{
Hyndman, R.J. and Koehler, A.B. (2006) "Another look at measures of forecast accuracy". International Journal of Forecasting, 22(4), 679-688. Hyndman, R.J. and Athanasopoulos, G. (2018) "Forecasting: principles and practice", 2nd ed., OTexts, Melbourne, Australia. Section 3.4 "Evaluating forecast accuracy". https://otexts.org/fpp2/accuracy.html.
}