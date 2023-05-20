\name{accuracy}
\alias{accuracy}
\alias{accuracy.marssPredict}
\alias{accuracy.marssMLE}

\title{Return accuracy metrics}
\usage{
\method{accuracy}{marssPredict}(object, x, test = NULL, type = "ytt1", verbose = FALSE, ...)
\method{accuracy}{marssMLE}(object, x, test = NULL, type = "ytt1", verbose = FALSE, ...)
}
\arguments{
  \item{object}{A \code{\link{marssMLE}} or \code{\link{marssPredict}} object}
  \item{x}{A matrix or data frame with data to test against the h steps of a forecast.}
  \item{test}{Which time steps in training data (data model fit to) to compute accuracy for.}
  \item{type}{type="ytt1" is the one-step-ahead predictions. type="ytT" is the fitted ytT predictions. The former are standardly used for training data prediction metrics.}  
  \item{verbose}{Show metrics for each time series of data.}
  \item{...}{Not used.}
}
\description{
 This is a method for the generic \code{accuracy} function in the generics package. It is written to mimic the output from the accuracy function in the forecast package. See that package for details.
 
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
The MASE calculation is scaled using MAE of the training set naive
forecasts which are simply \eqn{\mathbf{y}_{t-1}}{y(t-1)}.

For the training data, the metrics are shown for the one-step-ahead predictions by default (\code{type="ytt1"}). This is the prediction of \eqn{\mathbf{y}_t}{y(t)} conditioned on the data up to \eqn{t-1} (and the model estimated from all the data). With \code{type="ytT"}, you can compute the metrics for the fitted \code{ytT}, which is the expected value of new data at \eqn{t} conditioned on all the data. \code{type} does not affect test data (forecasts are past the end of the training data).
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12),]
train.dat <- dat[,1:12]
fit <- MARSS(train.dat, model = list(Z = factor(c("WA", "OR", "OR"))))

accuracy(fit)

# Compare to test data set
fr <- predict(fit, n.ahead=10)
test.dat <- dat[,13:22]
accuracy(fr, x=test.dat)
}
\references{
Hyndman, R.J. and Koehler, A.B. (2006) "Another look at measures of forecast accuracy". International Journal of Forecasting, 22(4), 679-688. 

Hyndman, R.J. and Athanasopoulos, G. (2018) "Forecasting: principles and practice", 2nd ed., OTexts, Melbourne, Australia. Section 3.4 "Evaluating forecast accuracy". \href{https://otexts.org/fpp2/accuracy.html}{https://otexts.org/fpp2/accuracy.html}.
}