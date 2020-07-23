\name{fitted.marssMLE}
\alias{fited.marssMLE}
\title{Return fitted values for X(t) and Y(t) in a MARSS model}
\description{

\code{fitted.marssMLE} is a general function to return the different types of fitted values for \eqn{\mathbf{x}_t}{x(t)} and \eqn{\mathbf{y}_t}{y(t)} in a MARSS model. 

\deqn{\mathbf{x}_{t+1} = \mathbf{B} \mathbf{x}_t + \mathbf{u} + \mathbf{C} \mathbf{c}_t + \mathbf{G} \mathbf{w}_t, \textrm{ where } \mathbf{w}_t \sim \textrm{MVN}(0,\mathbf{Q})}{x(t+1) = B x(t) + u + C c(t) + G w(t)}
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{D} \mathbf{d}_t + \mathbf{H} \mathbf{v}_t, \textrm{ where } \mathbf{v}_t \sim \textrm{MVN}(0,\mathbf{R})}{y(t) = Z x(t) + a + D d(t) + H v(t)}

The fitted values are the expectation of the right side of the MARSS equations. There are two considerations when computing the fitted values. 

First, one can include the error term  or exclude it. Within the state-space literature, "fitted" might mean either. If the error term is included (\code{fit.type="estimate"}), the Kalman filter and smoother estimates are returned for \eqn{x} and the mathematically analogous estimates are returned for \eqn{y}.  Read about these fitted values in \code{\link{MARSSest}}. If the error term is excluded (\code{fit.type="prediction"}), then the model predictions are computed for \eqn{y} (Z x(t) + A + covariates) and the one-step ahead predictions are computed for \eqn{x} (B x(t-1) + u + covariates). Read about these fitted values in \code{\link{MARSSpredict}}.


Second, one can do the calculations by conditioning on the data from 1 to \eqn{t-1} (\code{type="tt1"}), \eqn{t} (\code{type="tt"}), or \eqn{T} (\code{type="tT"}). Again within the state-space literature all three time conditionings are used.  

}

\usage{
\method{fitted}{marssMLE}(object, 
    type = c("ytT", "ytt", "ytt1", "xtT", "xtt", "xtt1"),
    interval = c("none", "confidence", "prediction"),
    level = 0.95, 
    fit.type = c("prediction", "estimate"),
    ...)
}
\arguments{
  \item{object}{A \code{\link{marssMLE}} object.}
  
  \item{type}{ Type of time conditioning to use.  Data 1 to \eqn{t-1} (\code{tt1}), data 1 to \eqn{t} (\code{tt}) or data 1 to \eqn{T} (\code{tT}). }
  
\item{interval}{ If \code{interval="confidence"}, then the standard error and confidence intervals are returned. If \code{interval="prediction"}, then the standard deviation and prediction intervals are returned.  }

  \item{level}{ Confidence level. alpha=1-level}

  \item{...}{Optional arguments.  If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not the Z matrix).}

}

\value{
A data frame with estimates, standard errors (or deviations), and intervals.
}

\details{

See \code{\link{MARSSest}} and \code{\link{MARSSpredict}} for background on the two types of fitted values and how they are calculated.

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

\references{ 
R. H. Shumway and D. S. Stoffer (2000).  Time series analysis and its applications. Edition 4. Springer-Verlag, New York.

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
MLEobj <- MARSS(dat)

# Make a plot of the estimated states
library(ggplot2)
d <- fitted(MLEobj, type = "xtT")
ggplot(data = d) +
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), linetype = 2, alpha = 0.3) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count")

# Make a plot of the estimates for the missing values
library(ggplot2)
d <- fitted(MLEobj, type = "ytT")
ggplot(data = d) +
  geom_point(aes(t, estimate)) +
  geom_line(aes(t, estimate)) +
  geom_point(aes(t, y), color = "blue") +
  geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(aes(t, pred.low), linetype = 2) +
  geom_line(aes(t, pred.high), linetype = 2) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval")

# Make a plot of the model estimate of y(t), i.e., put a line through the points
# Intervals are for new data not the blue dots 
# (which were used to fit the model so are not new)
library(ggplot2)
d <- predict(MLEobj, type = "ytT", interval="confidence", level=0.95)
d2 <- predict(MLEobj, type = "ytT", interval="prediction", level=0.95)
d$lwr <- d2$`Lo 95`
d$upr <- d2$`Hi 95`
ggplot(data = d) +
  geom_line(aes(t, estimate), size=1) +
  geom_point(aes(t, y), color = "blue") +
  geom_ribbon(aes(x = t, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.3) +
  geom_line(aes(t, pred.low), linetype = 2) +
  geom_line(aes(t, pred.high), linetype = 2) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval") +
  geom_text(x=15, y=7, label="The intervals are for \n new data not the blue dots")
}