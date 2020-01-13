\name{residuals.marssMLE}
\alias{residuals}
\alias{residuals.marssMLE}
\alias{residuals.MARSS}
\alias{residuals.marss}
\title{ MARSS Standardized Residuals }
\description{
  Calculates the standardized (or auxiliary) residuals sensu Harvey, Koopman and Penzer (1998).
}
\usage{
\method{residuals}{marssMLE}(object,..., Harvey=FALSE, normalize=FALSE)
}
\arguments{
  \item{object}{ An object of class \code{\link{marssMLE}}.}
  \item{...}{ Not used.}
  \item{Harvey}{ TRUE/FALSE }
  \item{normalize}{ TRUE/FALSE }
}
\value{
A list with the following components  
  \item{model.residuals}{ The smoothed model residuals E(V(t)|y(1:T),Theta), where Theta is the set of model parameters.  Sometimes called the smoothations. This is different than the Kalman filter innovations which are E(V(t)|y(1:t-1),Theta).}
  \item{state.residuals}{ The smoothed state residuals E(X(t+1)|y(1:T))-E(X(t+1)|E(x(t)|y(1:T))).}
  \item{residuals}{ The model residuals as a (n+m) x TT matrix with \code{model.residuals} on top and \code{state.residuals} below. \code{model.residuals} is \eqn{\hat(\eta_t)} on page 112 of Harvey, Koopman and Penzer (1998).  }
  \item{var.residuals}{ The variance of the model residuals and state residuals as a (n+m) x (n+m) x TT matrix with the model residuals in rows 1 to n.}
  \item{std.residuals}{ The standardized model residuals as a (n+m) x TT matrix. This is \code{residuals} divided by the square root of \code{var.residuals}---although the matrix equivalent of that equation.}
}
\details{

This function returns the conditional expected value (mean) and variance of the model and state residuals.  'conditional' means in this context, conditioned on the observed data and a set of parameters.  If there are no missing values in the data, E(Y|y)=y, and the model residual is the difference between y and the model prediction.  If the data are missing, then the model residual is 0 because E(Y)=model prediction and the difference between these is 0.  The variance however is non-zero.

If Harvey=TRUE, the function uses the algorithm on page 112 of Harvey, Koopman and Penzer (1998) to compute the conditional residuals and variance of the residuals.  If Harvey=FALSE, the function uses the equations in the technical report Holmes (2014).  The difference in the algorithm only concerns model residuals and their variance for missing values.  The residual variance for missing data is not normally used but is needed for leave-one-out analyses when a non-missing y(t) is removed and then its residual is compared the model residual variance at time t.

The residuals matrix (and se and variance) has a value for each time step.  The residuals in column t is the residual associated with the data at time t and the process (x) transition from x(t) to x(t+1). In the MARSS package, the process equation is written \eqn{x(t)=f(x(t-1))+w(t)}, and the residual in column t is then {v(t) w(t+1)}.  In many other texts, the process equation is written is written \eqn{x(t)=f(x(t-1))+w(t-1)} and the residual in column t is denoted {v(t) w(t)}.  Regardless, the meaning is the same; the state residual in column t is associated with the transition from x(t) to x(t+1) not the transition from x(t-1) to x(t). Thus the last state residual (and associated variance) is NA since there is no x(t+1).

\strong{model.residuals}

The model residuals \eqn{v_t} are the difference between the data and the predicted data at time t:
\deqn{ v_t = y_t - Z x_t - a}
In a state-space model, \eqn{x_t} is stochastic and the model residuals are a random variable.  \eqn{y_t} is also stochastic, though often observed unlike \eqn{x_t}.  The model residual random variable is:
\deqn{ V_t= Y_t - Z X_t - a}
The unconditional mean and variance of \eqn{V_t} is 0 and R.  This function (\code{residuals(MLEobj)}) returns the \emph{conditional} mean and variance of \eqn{V_t}.

\code{model.residuals} is the expected value of \eqn{V_t} conditioned on the data and parameter set \eqn{\Theta} (all the parameters including Z, a, B, u, R and Q).
\deqn{ E(V_t | y_1^T,\Theta) = E(Y_t | y_1^T,\Theta) - Z E(X_t | y_1^T, \Theta) - a}
If there are no missing data, this becomes
\deqn{ E(V_t | y_1^T,\Theta) = y_t - Z E(X_t | y_1^T, \Theta) - a}
\eqn{y_t} are data at time t and \eqn{E(X_t|y_1^T, \Theta)} is the Kalman smoother estimate of the states \eqn{x_t} at time t, i.e. the expected value of the states conditioned on all the data and the parameter set \eqn{\Theta}. Thus res1 and res2 in the code below will be the same.
\preformatted{dat = t(harborSeal)[2:3,]
MLEobj = MARSS(dat)
Z = coef(MLEobj, type="matrix")$Z
A = coef(MLEobj, type="matrix")$A
res1 = dat - Z \%*\% MLEobj$states - A \%*\% matrix(1,1,ncol(dat))
res2 = residuals(MLEobj)$model.residuals
}
The model.residuals for the missing data is 0 (in res2) because that is the expected value of \eqn{V_t} when data are missing.

\code{var.residuals} returned by the function is the conditional variance of the residuals. Rows 1 to n are the conditional variance of the model residuals. This is the variance of \eqn{V_t} conditioned on the data and the parameter set \eqn{\Theta}.  The unconditional variance (no data) would just be \eqn{R}.  See Holmes 2014.

If \code{Harvey=TRUE}, there will be no variance calculation for the missing y.  If \code{Harvey=FALSE}, the variance of the model residuals at the t with missing y are computed (via Holmes 2014). The interpretation is that although the data are missing, you can still imagine that these data exist and you need the variance of residuals. For example, if you are doing a leave-one-out cross-validation, the data exist and you need their variance because you are going to compute some diagnostics using the left-out data.  For outlier diagnostics and shock detection, the variances for the missing values are not needed.


\strong{state.residuals}

The state residuals \eqn{w_t} are the difference between the state at time t and the expected value of the state at time t given the state at time t-1:
\deqn{ w_t = x_t - B x_{t-1} - u}
Like the model residual, the state residual \eqn{w_t} is a random variable since \eqn{x_t} is a random variable:
\deqn{ W_t = X_t - B X_{t-1} - u}
The unconditional mean and variance of \eqn{W_t} is 0 and Q.  \code{residuals(MLEobj)} returns the \emph{conditional} mean and variance of \eqn{W_t}.

\code{state.residuals} is the expected value of \eqn{W_t} conditioned on the data (all the data 1 to T) and parameter set \eqn{\Theta}.
\deqn{ E(W_t|y_1^T,\Theta) = E(X_t|y_1^T,\Theta) - B E(X_{t-1}|y_1^T, \Theta) - u}
Thus res1 and res2 in the code below will be the same.
\preformatted{dat = t(harborSeal)[2:3,]
TT = ncol(dat)
MLEobj = MARSS(dat)
B = coef(MLEobj, type="matrix")$B
U = coef(MLEobj, type="matrix")$U
statest = MLEobj$states[,2:TT]
statestm1 = MLEobj$states[,1:(TT-1)]
res1 = statest - B \%*\% statestm1 - U \%*\% matrix(1,1,TT-1)
res2 = residuals(MLEobj)$state.residuals
}
The state residuals always exist since the expected value of the states exist without data and will be identical with \code{Harvey=TRUE} or \code{Harvey=FALSE}.  Generally speaking, \eqn{E(W_t|y_1^T)} is not 0 even if there are missing data.  Note that the state residual at the last time step (T) will be NA because it is the residual associated with x(T) to x(T+1) and T+1 is beyond the data.  Similarly, the variance matrix at the last time step will have NAs for the same reason.

\strong{standardized residuals}

\code{residuals.marssMLE} will return the standardized residuals sensu Harvey et al. (1998) for you in \code{std.residuals} for outlier and shock detection.  These are the model and state residuals scaled by the inverse square root of the missing values corrected variance of the residuals.  Note the standardized model residuals are set to NA when there are missing data (if there is no data point, there is no model residual). The standardized state residuals however always exist since the expected value of the states exist without data. The calculation of the standardized residuals for both the observations and states requires the full residuals variance matrix. Since the state residuals variance is NA at the last time step (since that is x(t) to x(t+1) and x(t+1) does not exist), the standarized residual in the last time step will be all NA.

The interpretation of the standardized residuals is not straight-forward when the Q and R variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in MVN(0,I) space.  For example, if v is 2x2 correlated errors with variance-covariance matrix R. The transformed residuals (from this function) for the i-th row of v is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the transformed residuals can look rather non-intuitive.

\strong{normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{V_t} and \eqn{W_t} are R and Q and the model is assumed to be written as
\deqn{ y_t = Z x_t + a + v_t}
\deqn{ x_t = B x_{t-1} + u + w_t}
Harvey et al (1998) writes the model as
\deqn{ y_t = Z x_t + a + Hv_t}
\deqn{ x_t = B x_{t-1} + u + Gw_t}
with the variance of \eqn{V_t} and \eqn{W_t} equal to I (identity).

\code{residuals.marssMLE} returns the residuals defined as in the first equations. To get the residuals defined as Harvey et al. (1998) define them (second equations), then use \code{normalize=TRUE}.  In that case the unconditional variance of residuals will be I instead of R and Q.  Note, that the `normalized' residuals are not the same as the `standardized' residuals.  In former, the unconditional residuals have a variance of I while in the latter it is the conditional residuals that have a variance of I.

}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSkem}} \code{\link{marssMLE}} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  MLEobj <- MARSS(dat)
  
  #state residuals
  state.resids1 <- residuals(MLEobj)$state.residuals
  #this is the same as
  states <- MLEobj$states
  Q <- coef(MLEobj,type="matrix")$Q
  state.resids2 <- states[,2:30]-states[,1:29]-matrix(coef(MLEobj,type="matrix")$U,2,29)
  #compare the two
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #normalize to variance of 1
  state.resids1 <- residuals(MLEobj, normalize=TRUE)$state.residuals
  state.resids2 <- (solve(t(chol(Q))) \%*\% state.resids2)
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #standardized (by joint variance) model & state residuals
  residuals(MLEobj)$std.residuals
}
\references{
Harvey, A., S. J. Koopman, and J. Penzer. 1998. Messy time series: a unified approach. Advances in Econometrics 13: 103-144  (see page 112-113).  Eqn 21 is the Kalman eqns.  Eqn 23 and 24 is the backward recursion to compute the smoothations.  This function uses the MARSSkf output for eqn 21 and then implements the backwards recursion in eqn 23 and eqn 24.  Pages 120-134 discuss the use of standardized residuals for outlier and structural break detection.

de Jong, P. and J. Penzer. 1998. Diagnosing shocks in time series. Journal of the American Statistical Association 93: 796-806.  This one shows the same equations; see eqn 6.  This paper mentions the scaling based on the inverse of the sqrt (chol) of the variance-covariance matrix for the residuals (model and state together).  This is in the right column, half-way down on page 800.

Koopman, S. J., N. Shephard, and J. A. Doornik. 1999. Statistical algorithms for models in state space using SsfPack 2.2. Econometrics Journal 2: 113-166. (see pages 147-148).

Harvey, A. and S. J. Koopman. 1992. Diagnostic checking of unobserved-components time series models. Journal of Business & Economic Statistics 4: 377-389.

Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045.  This report shows how to compute the residuals from the Kalman smoother output without using the Harvey et al. (1998) algorithm.  The variance of the hatvt for the missing values can then be computed.
}

