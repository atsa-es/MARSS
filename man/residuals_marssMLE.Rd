\name{residuals.marssMLE}
\alias{residuals.marssMLE}

\title{Return the model and state fitted values, residuals, and residual sigma}
\usage{
\method{residuals}{marssMLE}(x, type=c("smoothations", "innovations"),
                               standardization=c("Cholesky", "marginal"),
                               form=attr(x[["model"]], "form")[1])
}
\arguments{
  \item{x}{a \code{\link{marssMLE}} object}
  
  \item{type}{"Cholesky" means it is standardized by the Cholesky transformation of the full variance-covariance matrix of the model and state residuals. 
"marginal" means that the residual is standardized by its standard deviation, i.e. the square root of the value on the diagonal of the variance-covariance matrix of the model and state residuals.}

  \item{standardization}{"smoothations" means the fitted values are computed conditioned on all the data. See \code{\link{fitted}} with \code{type="ytT"} or \code{type="xtT"}. "innovations" means the fitted values are computed conditioned on the data up to t-1. See \code{\link{fitted}} with \code{type="ytt1"} or \code{type="xtt1"}.}

  \item{form}{If you want the function to use a different function than \code{residuals_form}.  This might be useful if you manually specified a DFA model and want to use \code{residuals_dfa} for rotating.}

}
\description{
\code{residuals.marssMLE} returns a data.frame (in tibble form) with fitted values, residuals, residual standard deviation (sigma), and standardized residuals. A residual is the difference between the "value" of the model (\eqn{y}) or  state (\eqn{x}) and the "fitted" value. 
}
\value{
A data frame with the following columns:

\item{.rownames}{Fitted values of observations or states. See details.}
\item{type}{model or state}
\item{t}{time step}
\item{value}{The data value if type="model" or the x estimate if type="state". }
\item{.fitted}{Fitted values of observations or states. See details.}
\item{.resids}{Model or states residuals. See details.}
\item{.sigma}{The standard error of the model or state residuals. Intervals for the residuals can be constructed from \code{.sigma} using qnorm(alpha/2)*.sigma + .fitted.}
\item{.std.resid}{Standardized residuals. Used for outlier detection. See \code{\link{MARSSresiduals}}.}

}
\details{

See \code{\link{MARSSresiduals}} for a discussion of the residuals calculations for MARSS models.

\strong{model residuals}

The model residuals are in the data frame with \code{.type="model"}.

The model residuals are the familar type of residuals, they are the difference between the data and the predictied value. So for the model residuals, the "value"" is the data (or is NA if data are missing). If \code{type="smoothations"}, the predicted value is the expected value of \eqn{y} conditioned on all the data, i.e. is computed using the smoothed estimated of \eqn{x} at time t (\code{xtT}).  If \code{type="innovations"}, the predicted value is the expected value of \eqn{y} conditioned on the data up to time t-1, i.e. is computed using the estimate of \eqn{x} at time t conditioned on the data up to time t-1 (\code{xtt1}). These are known as the one-step-ahead predictions and the residuals are known as the innovations. 

The standard errors help visualize how well the model fits to the data.  See \code{\link{fitted.marssMLE}} for a discussion of the calculation of the fitted values for the observations (the modeled values).  The standardized residuals can be used for outlier detection.  See \code{\link{MARSSresiduals}} and the chapter on shock detection in the MARSS User Guide.

\strong{state residuals}

The state residuals are in the data frame with \code{.type="state"}.

If you want the expected value of the states and an estimate of their standard errors (for confidence intervals), then \code{residuals()} is not what you want to use.  You want to use \code{\link{tidy.marssMLE}} to return the smoothed estimate of the state or you can find the states in the \code{states} element of the \code{\link{marssMLE}} object. 

The state residuals are only for state-space models. They are the difference between the state estimate at time t and the predicted value of the state given the estimate of the state at time t-1. For the state residuals, the "value"" depends on whether you are compute "smoothations" or "innovations" state residuals. If \code{type="smoothations"}, the value is the smoothed estimate (\eqn{xtT}) at time t and the predicted value is the expected value of \eqn{x} computed using the smoothed estimated of \eqn{x} at time t-1, so is \eqn{B xtT[,t-1] + u}).  If \code{type="innovations"}, the value is expected value of \eqn{x} conditioned on the data up to time t and the predicted value is computed using the expected value of \eqn{x} conditioned on the data up to time t-1, i.e. \eqn{B xtt[,t-1]} + u or \eqn{xtt1[,t]}.

These types of residuals are used for outlier detection or shock detection in the state process.  They are also used for model diagnostics.  See \code{\link{MARSSresiduals}} and read the references cited.

}


\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
MLEobj <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))

library(broom)
library(ggplot2)
theme_set(theme_bw())

# Make a plot of the observations and model fits
d <- residuals(MLEobj, interval = "confidence")
ggplot(data = d) +
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, y)) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.1) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count")

# Make a plot of xtT versus prediction of xt from xtT[t-1]
# This is NOT the estimate of the states with CIs. Use tidy() for that.
d <- residuals(MLEobj, type = "xtT")
ggplot(data = d) +
  geom_point(aes(t, xtT)) +
  geom_line(aes(x = t, .fitted)) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("xtT (points) and predition (line)")

}

\note{
Within the base code, a form-specific internal \code{residuals} function is called to allow the output to vary based on form: residuals_dfa, residuals_marss, residuals_marxss.
}