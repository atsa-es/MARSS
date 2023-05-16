\name{MARSS.out}
\alias{MARSS.out}
\keyword{experimental}

\title{ Output and Intervals for MARSS MLE objects }
\description{
  There are a variety of types of outputs and intervals for MARSS models. This help file summarizes them and tells you which function to use to get the output you want.
}
\usage{
MARSS.out() 
}
\details{
When looking at estimates and intervals for x and y, keep in mind that the MARSS parameters are always fixed. The uncertainty arise because x and y are stochastic not because of uncertainty in the parameter estimates. The latter type of uncertainty arises because we could have collected different data and those data would give us different parameter estimates. If you looked at the intervals on model fits (a regression line say) from a classic linear regression those intervals reflect parameter uncertainty. The intervals on x (states) and y (observations) in a MARSS model are fundamentally different. X is stochastic and unobserved. Y is stochastic and observed (at least partially).
\itemize{
\item  If you want the estimate of x conditioned on all the data 1 to T, then use \code{\link{tidy.marssMLE}} with \code{type="states", conditioning="T"}. This gives \eqn{E[X|all data]}.  This is the states estimates from the Kalman smoother. You will can also output these \code{\link{MARSSkf}}. The confidence interval is based on  \eqn{\var[X|all data]}, the conditional variance of X. That variance is \eqn{VtT[,,t]}.

\item  If you want the estimate of x conditioned on the data up to t-1, then use \code{\link{tidy.marssMLE}} with \code{type="states", conditioning="t-1"}. This is the one step ahead prediction, \eqn{B x(t-1) + u}, and its variance. \code{\link{tidy.marssMLE}} returns a prediction interval based on the variance  \eqn{Vtt1[,,t]=B Vtt[,,t-1] B^T + Q}. If you wanted the confidence interval on the expected value, that would use the variance \eqn{B Vtt[,,t-1] B^T}. \code{\link[=fitted.marssMLE]{fitted}} will return that.

\item  If you want the estimate of y conditioned on all the data, then use \code{\link{tidy.marssMLE}} with \code{type="observations", conditioning="T"}. This gives \eqn{E[Y|all data]} and if you have no missing data, \eqn{E[Y|data]} = \eqn{data} so you will just be returning your data. The confidence intervals will also be 0. For missing y values, you will see prediction intervals since \code{\link{tidy.marssMLE}} returns intervals based on \eqn{var[Y|all data]}.
  
\item  If you want the estimate of y conditioned on the data up to t-1, then use \code{\link{tidy.marssMLE}} with \code{type="observations", conditioning="t-1"}. This is the one step ahead prediction, \eqn{Z x(t|t-1) + A}, and its variance. \code{\link{tidy.marssMLE}} returns a prediction interval  based on the variance \eqn{Z Vtt1[,,t] Z^T + R}. If you wanted the confidence interval on the expected value of the prediction, that would use \eqn{Z Vtt1[,,t] Z^T}. \code{\link[=fitted.marssMLE]{fitted}} will return that.
}

\strong{One step ahead predictions} Use \code{conditioning="t-1"} for the one-step-ahead predictions of \eqn{\mathbf{y}_t}{y_t}, using only data up to t-1. Because the one-step-ahead predictions do not include \eqn{y_t} for calculation, you can use your data to evaluate the intervals. Thus in this case (only) it would make sense to plot your data on top of the prediction intervals. But in general use \code{\link{residuals.marssMLE}} when you want statistics (such as intervals) associated with your observed data.

\strong{Output that \code{fitted.marssMLE} provides and \code{tidy.marssMLE} does not}


\strong{fitted Y(t)} The fitted values of \eqn{\mathbf{Y}_t}{Y_t} is the expected value of the right side of the \eqn{\mathbf{Y}_t}{Y_t} equation WITHOUT the \eqn{v_t}. Use \code{\link{fitted.marssMLE}} and \code{\link[=fitted.marssMLE]{fitted}()} to get the fitted values, their residuals, and their residual variances (for outlier analysis).

\strong{Residual intervals} For residuals analysis (outlier detection), you want to use the "residual" intervals provided with \code{\link[=fitted.marssMLE]{fitted}()}. These are are the quantiles for the model residuals for data that were used to fit the model (and estimate xtT).  Read more about residuals analysis for MARSS models at \code{\link{MARSSresiduals}}. The residuals intervals are also what you want to evaluate data that you leave out, for example for a leave-one-out-cross-validation.

\strong{Leave-one-out analysis} The one-step-ahead prediction intervals are not what you want if you left some of your data out and want to evaluate the model residuals for the left-out-data. The residuals intervals are also what you want for this case (with the left-out-data set to NA).

\strong{Prediction intervals for new data} The one-step-ahead prediction intervals are also not what you want if you want to compute intervals for data between \eqn{t=1} to \eqn{T} that you didn't collect (so not in the future but contemporaneous with your data).  To compute prediction intervals for new data within \eqn{t=1} to \eqn{T}, you need to decide what you mean by 'new' and how these new data at time \eqn{t} relate to existing data at time \eqn{t}. Are the new data are i.i.d relative to the old data? So if you made an \eqn{\mathbf{R}}{R} matrix, it would be diagonal with the same value on the diagonal. Or are they missing data which might have a different variance and be correlated with the old data. In this case the \eqn{\mathbf{R}}{R} matrix might not be diagonal with one value on the diagonal.

\strong{Deciding on \code{fitted.marssMLE}, \code{predict.marssMLE} or \code{residuals.marssMLE} for estimates and intervals calculations.}

The function you want to use depends on what you are trying to do.
\itemize{
\item{Get the smoothed state estimates from the Kalman smoother}{ \code{fitted.marssMLE(x, type="xtT", fit.type="estimate")}. These are also the values are in \code{xtT} returned by \code{\link{MARSSkf}()}. }
\item{Get the contemporaneous state estimates from the Kalman filter}{ \code{fitted.marssMLE(x, type="xtt", fit.type="estimate")}.  These are also the values are in \code{xtt} returned by \code{\link{MARSSkf}()}. }
\item{Get the one-step ahead state estimates from the Kalman filter}{ \code{fitted.marssMLE(x, type="xtt1", fit.type="estimate")}.  These are also the values are in \code{xtt1} returned by \code{\link{MARSSkf}()}. }
\item{Get the model predictions for new data that might be generated from your model}{Use \code{fitted.marssMLE} with \code{type="ytT"} and \code{fit.type="prediction"}) or use \code{\link{predict.marssMLE}} with \code{type="ytT"}. This returns the model predictions (Z x(t)+A) for \eqn{Y_t} conditioned on all the data. Confidence intervals and prediction intervals can be returned. The former is the interval for the mean of new data and the latter is the interval for new data (not the mean but data themselves).}
\item{Compare your observed data to model predictions}{In this case, you want the distribution of the model residuals for the data. Use \code{\link{residuals.marssMLE}}.  }
\item{Get estimates and variance of missing data in your data set}{Use \code{fitted.marssMLE} with \code{type="ytT"} and \code{fit.type="estimate"}). The observed data will have an expected value equal to the observed data and variance of 0, while the missing data will have an expected value and variance conditioned on all the observed data. }
\item{Do a leave-one-out cross-validation}{ In this case, you want the distribution of the model residuals for those left-out values. Use \code{\link{residuals.marssMLE}}. You want the standard errors for the left-out data minus the predictions which is what  \code{\link{residuals.marssMLE}} gives in the \code{.type="model"} column.  }
\item{One-step-ahead predictions of y}{Use \code{fitted.marssMLE} with \code{type="ytt1"} and \code{fit.type="prediction"} or \code{predict.marssMLE} with \code{type="ytt1"}. Confidence (mean prediction) and prediction intervals (new data) are returned.  }
\item{One-step-ahead forecast errors}{Same as one-step-ahead above or use \code{\link{residuals.marssMLE}} with \code{type="innovations"}. }
}


}
\seealso{
\code{\link{MARSSkf}}, \code{\link{fitted.marssMLE}}, \code{\link{tidy.marssMLE}}, \code{\link{residuals.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
