\name{fitted.marssMLE}
\alias{fited.marssMLE}
\title{Return expected value of X(t) and Y(t) with summary information}
\description{

\code{fitted.marssMLE} returns the estimated state and observations. This is the expected value of the left-hand side of the \eqn{\mathbf{x}}{x} and \eqn{\mathbf{y}}{y} equations.

For the states, these are the outputs from the Kalman filter and smoother and are the expected values of \eqn{\mathbf{x}_t}{x(t)} conditioned on the data from 1 to \eqn{t-1}, \eqn{t}, or \eqn{T}. The former are the one-step ahead filter estimates and the latter are the smoothed state estimates. In the state-space literature, the state estimate is normally the expected values of the left-side of the \eqn{\mathbf{x}}{x}. However, if you need the expected value of the left-side, that is provided in \code{\link{predict.marssMLE}}.

In the state-space literature, the \eqn{\mathbf{y}}{y} "estimates" would normally refer to the expected value of the right-side of the \eqn{\mathbf{y}}{y} equation (i.e. the expected value of \eqn{\mathbf{Z} \mathbf{X}_t + \mathbf{a}}{Z X(t) + a}). That is provided in \code{\link{predict.marssMLE}}. \code{fitted.marssMLE} provides the expected value of the left-side of the \eqn{\mathbf{y}}{y} equation conditioned on the data from 1 to \eqn{t-1}, \eqn{t}, or \eqn{T}. These estimates are used to estimate missing values in the data. If \eqn{\mathbf{y}}{y} is multivariate, some \eqn{y} are missing at time \eqn{t} and some not, and \eqn{\mathbf{R}}{R} is non-diagonal, then the expected value of \eqn{\mathbf{y}_t}{y(t)} from the right-side of the \eqn{\mathbf{y}}{y} would be incorrect because it would take into account the information in the observed data at time \eqn{t} on the missing data at time \eqn{t} (except as it influences \eqn{\mathbf{x}_t}{x(t)}. 

Note, if there are no missing values, the expected value of \eqn{\mathbf{y}_t}{y(t)} (left-side) conditioned on the data from 1 to \eqn{t} or \eqn{T} is simply \eqn{\mathbf{y}_t}{y(t)}. The expectation is only useful when there are missing values for which an estimate is needed. The expectation of the left-side of the \eqn{\mathbf{y}}{y} is used in the general EM algorithm for the missing values case and the base function is \code{\link{MARSShatyt}}.
}
\usage{
fitted.marssMLE(object, 
    type = c("xtT", "xtt", "xtt1", "ytT", "ytt", "ytt1"),
    interval = c("none", "confidence"),
    level = 0.95, ...)
}
\arguments{
  \item{object}{A \code{\link{marssMLE}} object.}
  
  \item{type}{ Type of estimates to return. Smoothed states (xtT), one-step-ahead states (xtt1), contemporaneous states (xtt), the model ytT (Z xtT + A), the model ytt (Z xtt + A)the model ytt1 (Z xtt1 + A), the expected value of y conditioned on data 1 to \eqn{t-1} (ytt1), the expected value of y conditioned on data 1 to \eqn{t} (ytt), or the expected value of y conditioned on data 1 to \eqn{T} (ytT). See details.}
  
  \item{conf.int}{ Whether to compute confidence and prediction intervals on the estimates. }
  
\item{interval}{ If \code{interval="confidence"}, then the standard error and confidence intervals are returned. There are no prediction intervals for the expected value of the left-side of the MARSS equations only the right-side. }

  \item{level}{ Confidence level. alpha=1-level}

  \item{...}{Optional arguments.  If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not the Z matrix).}

}

\value{
A data frame with estimates, standard errors, and confidence intervals.
}

\details{

Below, X and Y refers to the random variable and x and y refer to a specific realization from this random variable.
  
\strong{state estimates (x)}

For \code{type="xtT}, \code{fitted.marssMLE} returns the confidence intervals of the state at time \eqn{t} conditioned on the data from 1 to \eqn{T} using the estimated model parameters as true values. These are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the \eqn{\mathbf{x}}{x}). However, these intervals do not include parameter uncertainty. If you want state residuals (for residuals analysis), use \code{\link{MARSSresiduals}} or \code{\link{residuals.marssMLE}}.

\strong{Quantiles} The state \eqn{\mathbf{X}_t}{X(t)} in a MARSS model has a conditional multivariate normal distribution, that can be computed from the model parameters and data. In Holmes (2012, Eqn. 11) notation, its expected value conditioned on the observed data and the model parameters \eqn{\Theta}{Theta} is \eqn{\tilde{\mathbf{x}}_t}{tildex(t)}. In \code{\link{MARSSkf}}, this is \code{xtT[,t]}. The variance of \eqn{\mathbf{X}_t}{X(t)} conditioned on the observed data and \eqn{\Theta} is \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} (\code{VtT[,,t]}). Note that \code{VtT[,,t] != B VtT[,,t-1] t(B) + Q}, which you might think by looking at the MARSS equations. That is because the variance of \eqn{\mathbf{W}_t}{W(t)} conditioned on the data (past, current and FUTURE) is not equal to  \eqn{\mathbf{Q}}{Q} (\eqn{\mathbf{Q}}{Q} is the unconditional variance).

\eqn{\tilde{\mathbf{x}}_t}{tildex(t)} (\code{xtT}) is an estimate of \eqn{\mathbf{x}_t}{x(t)} (the true value), and the standard error of that estimate is given by \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} (\code{VtT[,,t]}). Let \code{se.xt} denote the sqrt of the diagonal of \code{VtT}. The equation for the \eqn{\alpha/2}{alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.xt + xtT}). \eqn{\mathbf{x}_t}{x(t)} is multivariate and this interval is for one of the \eqn{x}'s in isolation. You could compute the m-dimensional confidence region for the multivariate \eqn{\mathbf{x}_t}{x(t)}, also, but \code{fitted.marssMLE} returns the univariate confidence intervals. 

The variance \code{VtT} gives information on the uncertainty of the true location of \eqn{\mathbf{x}_t}{x(t)} conditioned on the observed data. As more data are collected (or added to the analysis), this variance will shrink since the data, especially data at time \eqn{t}, increases the information about the locations of \eqn{\mathbf{x}_t}{x(t)}. This does not affect the estimation of the model parameters, those are fixed (we are assuming), but rather our information about the states at time \eqn{t}.

If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}

For \code{type="xtt} and \code{type=="xtt1"}, the calculations and interpretations of the intervals are the same but the conditioning is for data 1 to \eqn{t} or 1 to \eqn{t-1}. 


\strong{observation estimates (y)}

For \code{type="ytT"}, this returns the expected value and standard error of \eqn{\mathbf{Y}_t}{Y_t} (left-hand side of the \eqn{\mathbf{y}}{y} equation) conditioned on \eqn{\mathbf{Y}_t=y_t}{Y_t=y_t}. If you have no missing data, this just returns your data set. But you have missing data, this what you want in order to estimate the values of missing data in your data set. The expected value of \eqn{\mathbf{Y}_t|\mathbf{Y}=\mathbf{y}(1:T)}{Y|Y=y(1:T)} is in \code{ytT} in \code{\link{MARSShatyt}} output and the variance is \code{OtT-tcrossprod(ytT)} from the \code{\link{MARSShatyt}} output. 

The intervals reported by \code{fitted.marssMLE} for the missing values take into account all the information in the data, specifically the correlation with other data at time \eqn{t} if \eqn{\mathbf{R}}{R} is not diagonal. This is what you want to use for interpolating missing data. You do not want to use \code{predict.marssMLE()} as those predictions are for entirely new data sets and thus will ignore relevant information if \eqn{\mathbf{y}_t}{y_t} is multivariate, not all \eqn{\mathbf{y}_t}{y_t} are missing, and the \eqn{\mathbf{R}}{R} matrix is not diagonal.

The standard error and confidence interval for the expected value of the missing data along with the standard deviation and prediction interval for the missing data are reported. The former uses the variance of \eqn{E[\mathbf{Y}_t]}{E[Y(t)]} conditioned on the data while the latter uses variance of \eqn{\mathbf{Y}_t}{Y(t)} conditioned on the data. \code{\link{MARSShatyt}} returns these variances and expected values. See Holmes (2012) for a discussion of the derivation of expectation and variance of \eqn{\mathbf{Y}_t}{Y(t)} conditioned on the observed data (in the section 'Computing the expectations in the update equations').

For \code{type="ytt"} and \code{type="ytt1"}, only the estimates are provided. \code{\link{MARSShatyt}} does not return the necessary variances matrices for the standard errors for these cases.

\strong{Deciding on \code{fitted.marssMLE}, \code{predict.marssMLE} or \code{residuals.marssMLE} for the observation process}

The types of intervals you want for the \eqn{\mathbf{y}}{y} part of the MARSS equation depends on what you are trying to do.
\itemize{
\item{Get the model predictions for new data that might be generated from your model}{Use \code{\link{predict.marssMLE}} with \code{type="ytT"}. This returns the model predictions (Z x(t)+A) for \eqn{Y_t} conditioned on all the data. Confidence intervals and prediction intervals are returned. The former is the interval for the mean of new data and the latter is the interval for new data (not the mean but data themselves).}
\item{Get the distribution of new data at time \eqn{t} that would be generated by the process}{Same as above.}
\item{Compare your observed data to model predictions}{In this case, you want the distribution of the model residuals for the data. Use \code{\link{residuals.marssMLE}}.  }
\item{Get estimates and variance of missing data in your data set}{Use \code{fitted.marssMLE} with \code{type="ytT"}. The observed data will have an expected value equal to the observed data and variance of 0, while the missing data will have an expected value and variance conditioned on all the observed data. }
\item{Do a leave-one-out cross-validation}{ In this case, you want the distribution of the model residuals for those left-out values. Use \code{\link{residuals.marssMLE}}. You want the standard errors for the left-out data minus the predictions which is what  \code{\link{residuals.marssMLE}} gives in the \code{.type="model"} column.  }
\item{One-step-ahead predictions}{Use \code{predict.marssMLE} with \code{type="ytt1"}. Confidence (mean prediction) and prediction intervals (new data) are returned.  }
\item{One-step-ahead forecast errors}{Same as one-step-ahead or use \code{\link{residuals.marssMLE}} with \code{type="innovations"}. }
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