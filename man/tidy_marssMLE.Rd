\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters, states, and observations with summary information}
\usage{
  tidy.marssMLE(x, type = c("parameters", "states", "observations", "x", "y"),
                interval = TRUE,
                interval.type = c("confidence","quantile","residual"),
                smoothing = c("T", "t1", "t"),
                alpha = 0.05,
                interval.type = 
                form=attr(x[["model"]], "form")[1], ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type}{What you want estimates and intervals for. Parameters, states (x) or observations(y).}
  
  \item{interval}{ Whether to compute intervals on the estimates. }
  
  \item{interval.type}{ Type of interval if interval==TRUE. If type=="parmeters", the only allowed interval type is "confidence". If type is states, the only option is "quantile". If type is observations, the options are "residual", "prediction" or "quantile". Default is "residual". See details for discussion of intervals for states and model fits.}
  
  \item{smoothing}{ If type is states or observations, how should the estimates be computed: using all data 1:T (T), using only data up to t-1 (t1), or using only data up to t (t). If smoothing="t1", then you get the one-step-ahead predictions (and their quantiles). If smoothing="T", then you get the smoothed state and model fits (and their quantiles). }

  \item{alpha}{ Confidence level if confidence interval is returned otherwise the quantile level for quantiles and prediction intervals.}

  \item{form}{If you want the tidy function to use a different form than that specified in attr(x$model, "form").  Useful if you have a DFA model that you manually set up, which does not have the form attribute set. Normally just ignore and let the function use the "form" set in the attributes.}
  
  \item{...}{Optional arguments.  If \code{interval=TRUE}, then arguments to specify how CIs are computed can be passed in. See details and \code{\link{MARSSparamCIs}}. If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not Z matrix).}

}
\description{

The tidy function is designed to work with the \code{broom} package and you will need to load that package if you want to call \code{tidy(fit)} instead of \code{tidy.marssMLE(fit)}. The \code{tidy.marssMLE} functions returns the parameter estimates. It will also return the estimated state (X) at time t conditioned on all the data, the data up to t-1 (one-step-ahead), or the data up t. Lastly, it will return the fitted Y estimate at time t, conditioned on all the data, the data up to t-1 (one-step-ahead), or the data up t.
}
\return{
\strong{type="parameters"} A data.frame with estimates, sample standard errors, and confidence intervals (unless interval=FALSE).

\strong{type="states"} A data.frame with estimates of x with the standard deviation and quantiles (unless interval=FALSE) computed using the estimated model parameters. The state returned depends on the smoothing specification. If smoothing="T", the expected value of x(t) conditioned on all the data is returned. If smoothing="t1", it is conditioned on the data up to t-1. These are called one-step-ahead predictions. If smoothing="t", it is conditioned on data up to t. 

\strong{type in "observations", "y"} A data.frame of the fitted y. This is the  the right side of the y equation (minus the \eqn{v_t}). Intervals are included if interval=TRUE. For residuals analysis (outlier detection), you want to use the "residual" intervals. If you want the estimate and distribution of your missing data (if any), you also want the "residual" intervals. If you want the distribution of new data sets, you want the "quantile" intervals---this would be the sort of distribution you get if you simulated new data sets from your estimated model. Are these prediction intervals? No. Prediction intervals include the uncertainty in the estimated model parameters.
Confidence and prediction intervals for fitted y are not provided at this time. See details below. 
}

\details{

\strong{type=="parameters"}

  If \code{type=="parameters"}, this returns a data.frame with the estimated parameters  of a MARSS model with, optionally, standard errors and confidence intervals. This assembles information available via the print and coef functions into a data.frame that summarizes the estimates.
  
  If interval=TRUE, \code{\link{MARSSparamCIs}} will be run to add confidence intervals to the fitted model object if these are not already added to the model object.  The default CIs are calculated using a analytically computed Hessian matrix.  This can be changed by passing in optional arguments for \code{\link{MARSSparamCIs}}. 
  
\strong{type="states"}

Currently \code{tidy.marssMLE} does not return the confidence intervals for the estimated states. Instead only the quantiles of the states at time t using the estimated model parameters are returned. These are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the x). However, see the comments below on confidence intervals for the expected value of the states.

\strong{Quantiles} The state \eqn{X_t} in a MARSS model is a random variable. In Holmes (2012, Eqn. 11) notation, its expected value conditioned on all the observed data (1:T) and the model parameters \eqn{\Theta} is \eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t}. In \code{\link{MARSSkf}} notation, this is \code{xtT[,t]}. The variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the observed data and \eqn{\Theta} is \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). Analogously, we can compute the expected value and variance of \eqn{\mathbf{X}_t}{X_t} conditioned on only the data up to t-1. Those values are \code{xtt1[,t]} and \code{Vtt1[,,t]}. We can also compute the expected value and variance of \eqn{\mathbf{X}_t}{X_t} conditioned on only the data up to t. Those are \code{xtt[,t]} and \code{Vtt[,,t]}. All these statistics are available from the Kalman filter and smoother function \code{\link{MARSSkf}}.

For states, \code{tidy.marssMLE} returns the empirical, i.e. estimated, quantiles of \eqn{\mathbf{X}_t}{X_t} based on one of the estimated conditional variance-covariance matrices and corresponding expected value of \eqn{\mathbf{X}_t}{X_t}. These are the estimated values since we do not know the true \eqn{\Theta} but rather its estimate \eqn{\hat{\Theta}}. Let's denote the former as \code{Vx} (which is either \code{VtT}, \code{Vtt1} or \code{Vtt}) and the latter as \code{Ex} (which is either \code{xtT}, \code{xtt1}, or \code{xtt}). Let \code{sd.xt} denote the sqrt of the diagonal of \code{Vx}.  The equation for the \eqn{\alpha/2} quantile is (\code{qnorm(alpha/2)*sd.xt + Ex}). To emphasize that these are quantile intervals on \eqn{X_t} and not confidence intervals on \eqn{E[X_t]}, the column headings are \code{.quant.low} and \code{.quant.up} and the column on its standard deviation is labelled \code{.std.dev} not \code{.std.error}. A standard error is a measure of the variability of some estimated value from a sample and this variability shrinks as the sample size increases. The standard deviation does not shrink with more data; it is simply a function of the model parameters.
  
If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}

What about confidence and prediction intervals? These are intervals that incorporate the uncertainty in the estimated model parameters. If the interval is for \eqn{E[X_t]} then it is called a confidence interval. If it is for \eqn{X_t} (not its expected value) then it is called a prediction interval. Currently these are not included in the MARSS R package. But see the notes below for how one would compute these.

\strong{type="observations"}

The fitted Y means the right side of the y equation in the MARSS model without the \eqn{\mathbf{v_t}}{v_t}. The fitted y can be computed using either xtT, xtt1, or xtt depending on whether you want the fitted y using all the data, using only data up to t-1, or using data up to t. The model residuals are difference between the y (data) and the fitted y.  For residuals analysis (outlier detection), you want to use the "residual" intervals. These are are the quantiles for the model residuals for data that were used to fit the model (and estimate xtT).  Read more about residuals analysis for MARSS models at \code{\link{residuals.marssMLE}}. If you want the distribution of new data sets, you want the "quantile" intervals. "quantile" intervals are the unconditional quantiles of data that can be generated from the model. These are for de novo data like you get if you simulated from your estimated model. 

What about getting the estimate and distribution for missing data? If your y is missing data, you might want an estimate of those values. In this case, you want to include all the other data that you have, that is you want the estimate to be conditioned on all the data. So in this case you also want the "residual" intervals. If you use the "quantile" intervals, then you are not using all the information in your observed data. For example, say you have 3 y observations for one x and you are missing one of the y observations at time t. The residual interval for the missing y will include the information from the other 2 y but the quantile interval will not.

Why are residual and prediction intervals different? Here is a simple example. Say you observe 100 data points (y) from a Normal distribution with mean u and variance q. The estimate of u is the mean of y. Let's call that \eqn{\hat{u}}. The estimate of q is the variance of the sample residuals: \eqn{var(y-\hat{u})}, but q is not the variance of \eqn{y-\hat{u}}! q is the variance of \eqn{y-u}. The variance of \eqn{y-\hat{u}} is smaller than q because y were used to estimate \eqn{\hat{u}}. The y are closer (on average) to the sample mean (mean of y) than y is to u (population mean). That is why residual intervals are narrower than prediction intervals. If you are doing residual analysis (looking for outliers), then you uses the residual intervals.

Currently confidence and prediction intervals on fitted y and y are not included in the MARSS R package. Read the notes about how to conceptually compute these for the states in the notes below. The process is similar for the fitted y and y.
}
\references{ 
R. H. Shumway and D. S. Stoffer (2000).  Time series analysis and its applications. Edition 4. Springer-Verlag, New York.

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]
}
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11,12),]
  MLEobj <- MARSS(dat, model=list(Z=factor(c("WA","OR","OR"))))

library(broom)
library(ggplot2)

# A data frame of the estimated parameters
 tidy(MLEobj)
 
# Make a plot of the estimated states
# Don't use augment.  States are not data.
d <- tidy(MLEobj, type="states")
ggplot(data = d) + 
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1) +
  facet_grid(~term) +
  xlab("Time Step") + ylab("Count")
}
\note{
\section{Confidence intervals for the estimate of \eqn{E[X_t]}}

Although \code{tidy.marssMLE} does not return these CIs currently, computation of them is conceptually straight forward. Let's start with a review of computing the standard error of a sample statistic.

\strong{Standard error} The standard error of the estimate of \eqn{E[X_t]}---whether conditioned on all the data (\code{xtT}), data up to t-1 (\code{xtt1}) or data up to t (\code{xtt})---is a measure of the how far the estimate (using the estimate \eqn{\hat{\Theta}}) is from the true \eqn{E[X_t]} using the true \eqn{\Theta}. The standard error of the estimate of \eqn{E[X_t]}} uses the information about our uncertainty in \eqn{\hat{\Theta}}. With more data, i.e. a longer time series or multiple \eqn{y_t} at each t, the uncertainty in our estimated \eqn{\Theta} shrinks and the standard error of the estimate of \eqn{E[X_t]}} shrinks. This is the standard notion of the standard error of a sample statistic.

So how do you get the standard error of the estimate of \eqn{E[X_t]}}? This can be computed analytically for MARSS models, but let's imagine doing this with a parametric bootstrap.  If \eqn{\Theta} is known (say you are doing a computer simulation and specify \eqn{\Theta}), the standard error can be computed from the variance-covariance matrix of \eqn{\hat{\Theta}}. The latter is computed from the Fisher Information matrix (see \code{\link{MARSShessian}}). Then we use a parametric bootstrap to get the distribution of estimates of \eqn{E[X_t]} from bootstrapped time series of a specific length, number of missing values, and structure. 

The process is: draw \eqn{\hat{\Theta}_b}'s from the variance-covariance matrix, simulate a time series from that \eqn{\hat{\Theta}_b}, compute \code{xtT} (or \code{xtt} or \code{xtt1}), and repeat 1000s of times. The function \code{\link{MARSSsimulate}} will do most of this process for you, everything but running the bootstrap data through \code{\link{MARSSkf}} to get \code{xtT} (or \code{xtt} or \code{xtt1}).  That is the general concept of how one gets the distribution of sample statistics and computes their standard error when working in a frequentist framework with maximum-likelihood parameter estimates (the statistical framework that the MARSS package operates within). From that the distribution of the sample statistics (in our case the xtT's), you can get their standard error at time t, which is the square root of the variance of the estimates at time t.

So that was the case where \eqn{\Theta} is known. If we only have an estimate of \eqn{\Theta} from our one data set, then the above procedure yields the \strong{sample} standard error of the estimate of \eqn{E[X_t]}, i.e. \eqn{E[X_t]|\hat{\Theta}}. 

\strong{Confidence interval} A confidence interval is associated with a specific sample. It is an interval on a statistic (like mean or in our case \eqn{E[X_t]}), computed from a sample. This interval is constructed such that it covers the true value of the statistic for \eqn{1-\alpha} percent of all dat samples. So for our bootstrapped data sets from the true \eqn{\Theta}, we could construct an interval for each bootstrap data set and if constructed correctly, the interval would include \eqn{E[X_t]|\Theta} in \eqn{1-\alpha} percent of the bootstrapped data sets (\eqn{|\Theta} was added to be explicit about which \eqn{E[X_t]} we are talking about). There are an infinite number of ways that you could contruct this interval. The one used overwhelmingly in frequentist statistics is one that has nice compact properties.

How do we construct the confidence interval? MARSS models are in the class of multivariate Gaussian models and we can compute the confidence interval from the sample standard error of the estimate of \eqn{E[X_t]|\hat{\Theta}}. \eqn{E[X_t]} is multivariate. We can compute univariate confidence intervals that apply to only one of the X, we can compute multivariate confidence regions, or we can compute the simultaneous univariate confidence regions. The latter means an interval such that if the interval for one X in \eqn{E[X_t]} covers then all cover (all or none).  

The univariate confidence intervals are the easy ones to compute. Let's use \code{se.xt} for the sample standard error of the estimate of \eqn{E[X_t]|\hat{\Theta}}---the square root of the diagonal of the variance-covariance matrix for the estimates of \eqn{E[X_t]}. The approximate CIs can be computed by \code{qnorm(alpha/2)*se.xt + Ex}, where \code{Ex} is \code{xtT}, \code{xtt} or \code{xtt1} depending on what \eqn{E[X_t]} we are interested in.

Using \code{qnorm()} is an approximation and is treating \eqn{\Theta} as known instead of as an estimate. The Normal approximation does not have heavy enough tails and leads to intervals that are too narrow. In a univariate Normal case, you'd use a t-distribution and \code{qt(alpha/2, df)}. For a MARSS model, the degrees of freedom can be different for different x_t in X_t (different rows) and just depends on your model structure. You could construct bootstrap CIs but that is incredibly slow since you have to bootstrap your bootstraps. This type of Normal approximation for confidence intervals is commonly used for more complex models, but just aware that it does mean that the confidence intervals are a bit narrow. 
}