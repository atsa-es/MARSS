\name{residuals.marssMLE}
\alias{residuals.marssMLE}

\title{Model and state fitted values, residuals, and residual sigma}
\usage{
\method{residuals}{marssMLE}(object, ..., 
           type=c("tt1", "tT", "tt"),
           standardization=c("Cholesky", "marginal"),
           form=attr(object[["model"]], "form")[1])
}
\arguments{
  \item{object}{a \code{\link{marssMLE}} object}
  
  \item{type}{"tt1" means innovations residuals. The fitted values are computed conditioned on the data up to t-1. See \code{\link[=fitted.marssMLE]{fitted}} with \code{type="ytt1"} or \code{type="xtt1"}. "tT" means smoothation residuals. the fitted values are computed conditioned on all the data. See \code{\link[=fitted.marssMLE]{fitted}} with \code{type="ytT"} or \code{type="xtT"}. "tt" means contemporaneous residuals. The fitted values are computed conditioned on the data up to t. In MARSS functions, estimates at time t conditioned on data 1 to T are denoted "tT", conditioned on the data from 1 to t-1 are denoted "tt1" and conditioned on data 1 to t are "tt".}

  \item{standardization}{"Cholesky" means it is standardized by the lower triangle of the Cholesky transformation of the full variance-covariance matrix of the model and state residuals. 
"marginal" means that the residual is standardized by its standard deviation, i.e. the square root of the value on the diagonal of the variance-covariance matrix of the model and state residuals. }


  \item{form}{If you want the function to use a different function than \code{residuals_form}.  This might be useful if you manually specified a DFA model and want to use \code{residuals_dfa} for rotating. This can be ignored in almost all cases.}
  
  \item{...}{Not used.}

}
\description{
\code{residuals.marssMLE} returns a data frame with fitted values, residuals, residual standard deviation (sigma), and standardized residuals. A residual is the difference between the "value" of the model (\eqn{y}) or  state (\eqn{x}) and the "fitted" value.  At time t (in the returned data frame), the model residuals are for time t. For the the state residuals, the residual is for the transition from t to t+1 following the convention in Harvey, Koopman and Penzer (1998). For the the state innovation residuals, this means that \code{state.residual[,t]} is for the transition from \eqn{t} to \eqn{t+1} and is conditioned on data 1 to \eqn{t} while \code{model.residual[,t]} is is conditioned on data 1 to \eqn{t-1}. State innovation residuals are not normally used while state smoothation residuals are used in trend outlier analysis.
}
\value{
A data frame with the following columns:

\item{.rownames}{Fitted values of observations or states. See details.}
\item{type}{model or state}
\item{t}{time step}
\item{value}{The data value if type="model" or the x estimate if type="state". }
\item{.fitted}{Model predicted values of observations or states. See details.}
\item{.resids}{Model or states residuals. See details.}
\item{.sigma}{The standard error of the model or state residuals. Intervals for the residuals can be constructed from \code{.sigma} using qnorm(alpha/2)*.sigma + .fitted.}
\item{.std.resid}{Standardized residuals. See \code{\link{MARSSresiduals}} for a discussion of residual standardization.}

}
\details{

See \code{\link{MARSSresiduals}} for a discussion of the residuals calculations for MARSS models.

\strong{model residuals}

The model residuals are in the data frame with \code{.type="model"}.

The model residuals are the familiar type of residuals, they are the difference between the data at time t and the predicted value at time t. For the model residuals, the "value"" is the data (or NA if data are missing). If \code{type="tT"}, the predicted value is the expected value of \eqn{y} conditioned on all the data, i.e. is computed using the smoothed estimated of \eqn{x} at time t (\code{xtT}).  If \code{type="tt1"}, the predicted value is the expected value of \eqn{y} conditioned on the data up to time t-1, i.e. is computed using the estimate of \eqn{x} at time t conditioned on the data up to time t-1 (\code{xtt1}). These are known as the one-step-ahead predictions and the residuals are known as the innovations. 

The standard errors help visualize how well the model fits to the data.  See \code{\link[=fitted.marssMLE]{fitted}} for a discussion of the calculation of the model predictions for the observations.  The standardized smoothation residuals can be used for outlier detection.  See the references in \code{\link{MARSSresiduals}} and the chapter on shock detection in the MARSS User Guide.

\strong{state residuals}

The state residuals are in the data frame with \code{.type="state"}.

If you want the expected value of the states and an estimate of their standard errors (for confidence intervals), then \code{residuals()} is not what you want to use.  You want to use \code{\link[=tsSmooth.marssMLE]{tsSmooth}(..., type="xtT")} to return the smoothed estimate of the state or you can find the states in the \code{states} element of the \code{\link{marssMLE}} object returned by a MARSS() call. For the one-step-ahead state estimates, use \code{\link[=tsSmooth.marssMLE]{tsSmooth}(..., type="xtt1")}.

The state residuals are only for state-space models. At time \eqn{t}, the state residuals are the difference between the state estimate at time t+1 and the predicted value of the state at time \eqn{t+1} given the estimate of the state at time \eqn{t}. For smoothation state residuals, this is 
\deqn{ \hat{\mathbf{w}}_{t+1} = \tilde{\mathbf{x}}_{t+1}^T - \mathbf{B}\tilde{\mathbf{x}}_{t}^T - \mathbf{u} }{ hatw_t+1 = xt+1T - B xtT - u}
For innovations state residuals, this is 
\deqn{ \hat{\mathbf{w}}_{t+1} = \tilde{\mathbf{x}}_{t+1}^{t+1} - \mathbf{B}\tilde{\mathbf{x}}_{t}^t - \mathbf{u} }{ hatw_t+1 = xt+1t+1 - B xtt - u}

Smoothation state sresiduals are used for outlier detection or shock detection in the state process.  See \code{\link{MARSSresiduals}} and read the references cited. Note that the state residual at time T (the last time step) is NA since this would be the transition from T to T+1 (past the end of the data).

}


\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
fit <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))

library(ggplot2)
theme_set(theme_bw())

# Make a plot of the model residuals (innovations) with intervals
d <- residuals(fit, type="tt1")
d$.conf.low <- d$.fitted+qnorm(0.05/2)*d$.sigma
d$.conf.up <- d$.fitted-qnorm(0.05/2)*d$.sigma
ggplot(data = d) +
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, value), na.rm=TRUE) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.1) +
  ggtitle("Model residuals (innovations)") +
  xlab("Time Step") + ylab("Count")

# Make a plot of the smoothed residuals with intervals
d <- residuals(fit, type="tT")
d$.conf.low <- d$.fitted+qnorm(0.05/2)*d$.sigma
d$.conf.up <- d$.fitted-qnorm(0.05/2)*d$.sigma
ggplot(data = d) +
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, value), na.rm=TRUE) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.1) +
  facet_grid(~.rownames) +
  ggtitle("Smoothations") +
  xlab("Time Step") + ylab("Count")

# Make a plot of xtT versus prediction of xt from xtT[t-1]
# This is NOT the estimate of the smoothed states with CIs. Use tsSmooth() for that.
d <- residuals(fit, type = "tT")
ggplot(data = d) +
  geom_point(aes(t, value), na.rm=TRUE) +
  geom_line(aes(x = t, .fitted), color="blue") +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("xtT (points) and prediction (line)")

}

\references{
Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 

See also the discussion and references in \code{\link{MARSSresiduals.tT}}, \code{\link{MARSSresiduals.tt1}} and \code{\link{MARSSresiduals.tt}}.
}
