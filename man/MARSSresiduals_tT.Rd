\name{MARSSresiduals.tT}
\alias{MARSSresiduals.tT}

\title{ MARSS Smoothed Residuals }
\description{
  Calculates the standardized (or auxiliary) smoothed residuals sensu Harvey, Koopman and Penzer (1998). The expected values and variance for missing (or left-out) data are also returned (Holmes 2014). Not exported. Access this function with \code{residuals(object, conditioning="T")}.
}
\usage{
MARSSresiduals.tT(object, Harvey=FALSE, normalize=FALSE)
}
\arguments{
  \item{object}{ An object of class \code{\link{marssMLE}}.}
  \item{Harvey}{ TRUE/FALSE. Use the Harvey et al. (1998) algorithm or use the Holmes (2014) algorithm. The values are the same except for missing values. }
  \item{normalize}{ TRUE/FALSE }
}
\value{
A list with the following components  
  \item{model.residuals}{ The the observed smoothed model residuals: data minus the model predictions conditioned on all observed data. This is different than the Kalman filter innovations which use on the data up to time \eqn{t-1} for the predictions. See details. }
  \item{state.residuals}{ The smoothed state residuals \eqn{\tilde{x}_{t+1}^T - \mathbf{Z}\tilde{x}_t^T - \mathbf{u}}{E(X_{t+1}|y(1:T))-E(X_{t+1}|E(x_t|y(1:T)))}.}
  \item{residuals}{ The residuals conditioned on the observed data. Returned as a (n+m) x TT matrix with \code{model.residuals} in rows 1 to n and \code{state.residuals} in rows n+1 to n+m.  NAs will appear in rows 1 to n is the places where data are missing. }
  \item{var.residuals}{ The joint variance of the model and state residuals conditioned on observed data (with the same location of missing values as the observed data). Returned as a (n+m) x (n+m) x TT matrix. For Harvey=FALSE, this is Holmes (2014) equation 57. For Harvey=TRUE, this is the residual variance in eqn. 24, page 113, in Harvey et al. (1998). They are identical except for missing values, for those Harvey=TRUE returns 0s.}
  \item{std.residuals}{ The Cholesky standardized residuals as a (n+m) x TT matrix. This is \code{residuals} multiplied by the inverse of the Cholesky decomposition of \code{var.residuals}. The model standardized residuals associated with the missing data are replaced with NA. This for convenience for residuals diagnostics. }
  \item{mar.residuals}{ The marginal standardized residuals as a (n+m) x TT matrix. This is \code{residuals} multiplied by the inverse of the diagonal matrix formed by the square-root of the diagonal of \code{var.residuals}.  The model marginal residuals associated with the missing data are replaced with NA. This for convenience for residuals diagnostics. }
    \item{E.obs.residuals}{ The expected value of the model residuals conditioned on the observed data. Returned as a n x TT matrix.  For observed data, this will be the observed residuals. For unobserved data, this will be 0 if \eqn{\mathbf{R}}{R} is diagonal but non-zero if \eqn{\mathbf{R}}{R} is non-diagonal. See details. }
    \item{var.obs.residuals}{ The variance value of the model residuals conditioned on the observed data. Returned as a n x n x TT matrix.  For observed data, this will be 0. See details. }
  \item{msg}{ Any warning messages. This will be printed unless Object$control$trace = -1 (suppress all error messages). }

}
\details{

This function returns the raw, Cholesky standardized and marginal standardized smoothed model and state residuals.  'smoothed' means conditioned on all the observed data and a set of parameters. These are the residuals presented in Harvey, Koopman and Penzer (1998) pages 112-113, with the addition of the values for unobserved data (Holmes 2014).  If Harvey=TRUE, the function uses the algorithm on page 112 of Harvey, Koopman and Penzer (1998) to compute the conditional residuals and variance of the residuals.  If Harvey=FALSE, the function uses the equations in the technical report Holmes (2014).  

The residuals matrix has a value for each time step.  The residuals in column \eqn{t} rows 1 to n is the model residual associated with the data at time \eqn{t}. The residuals in rows n+1 to n+m are the state residuals associated with the transition from \eqn{\mathbf{x}_t}{x(t)} to \eqn{\mathbf{x}_{t+1}}{x(t+1)}, not the transition from \eqn{\mathbf{x}_t}{x(t-1)} to \eqn{\mathbf{x}_{t+1}}{x(t)}. Because \eqn{\mathbf{x}_{t+1}}{x(t+1)} does not exist at time \eqn{T}, the state residuals and associated variances at time \eqn{T} are NA.

Below the conditional residuals and their variance are discussed. The random variables are capitalized and the realizations from the random variables are lower case. The random variables are \eqn{\mathbf{X}}{X}, \eqn{\mathbf{Y}}{Y}, \eqn{\mathbf{V}}{V} and \eqn{\mathbf{W}}{W}. There are two types of \eqn{\mathbf{Y}}{Y}. The observed \eqn{\mathbf{Y}}{Y} that are used to estimate the states \eqn{\mathbf{x}}{x}. These are termed \eqn{\mathbf{Y}^{(1)}}{Y(1)}. The unobserved  \eqn{\mathbf{Y}}{Y} are termed \eqn{\mathbf{Y}^{(2)}}{Y(2)}. These are not used to estimate the states \eqn{\mathbf{x}}{x} and we may or may not know the values of \eqn{\mathbf{y}^{(2)}}{y(2)}. Typically we treat \eqn{\mathbf{Y}^{(2)}}{Y(2)} as unknown but it may be known but we did not include it in our model fitting.  Note that the model parameters \eqn{\Theta}{Theta} are treated as fixed or known. The 'fitting' does not involve estimating \eqn{\Theta}{Theta}; it involves estimating \eqn{\mathbf{x}}{x}. All MARSS parameters can be time-varying but the \eqn{t} subscripts are left off parameters to reduce clutter.

\strong{Model residuals}

\eqn{\mathbf{v}_t}{v_t} is the difference between the data and the predicted data at time \eqn{t} given \eqn{\mathbf{x}_t}{x_t}:
\deqn{ \mathbf{v}_t = \mathbf{y}_t - \mathbf{Z} \mathbf{x}_t - \mathbf{a} }{ v_t = y_t - Z x_t - a}
The observed model residuals \eqn{\hat{\mathbf{v}}_t}{hatv_t} are the difference between the observed data and the predicted data at time \eqn{t} using the fitted model. \code{MARSSresiduals.tT} fits the model using all the data. So
\deqn{ \hat{\mathbf{v}}_t = \mathbf{y}_t - \mathbf{Z}\tilde{x}_t^T - \mathbf{a} }{ hatv_t = y_t - Z xtT - a}
where \eqn{\tilde{x}_t^T}{xtT} is the expected value of \eqn{\mathbf{X}_t}{X_t} conditioned on the data from 1 to \eqn{T} (all the data), i.e. the Kalman smoother estimate of the states at time \eqn{t}. \eqn{\mathbf{y}_t}{y_t} are your data and missing values will appear as NA in the observed model residuals. These are returned as \code{model.residuals} and rows 1 to \eqn{n} of \code{residuals}.

\code{res1} and \code{res2} in the code below will be the same.
\preformatted{dat = t(harborSeal)[2:3,]
MLEobj = MARSS(dat)
Z = coef(MLEobj, type="matrix")$Z
A = coef(MLEobj, type="matrix")$A
res1 = dat - Z \%*\% MLEobj$states - A \%*\% matrix(1,1,ncol(dat))
res2 = residuals(MLEobj)$model.residuals
}

\strong{state.residuals}

The state residuals \eqn{\mathbf{w}_t}{w_t} are the difference between the state at time \eqn{t} and the expected value of the state at time \eqn{t} given the state at time \eqn{t-1}:
\deqn{ \mathbf{w}_t = \mathbf{x}_t - \mathbf{B} \mathbf{x}_{t-1} - \mathbf{u} }{ w_t = x_t - B x_{t-1} - u}
The estimated state residuals \eqn{\hat{\mathbf{w}}_t}{hatw_t} are the difference between estimate of \eqn{\mathbf{x}_t}{x_t} minus the estimate using \eqn{\mathbf{x}_{t-1}}{x_{t-1}}. 
\deqn{ \hat{\mathbf{w}}_t = \tilde{x}_t^T - \mathbf{B}\tilde{x}_{t-1}^T - \mathbf{u} }{ hatw_t = xtT - B xt1T - u}
where \eqn{\tilde{x}_t^T}{xtT} is the Kalman smoother estimate of the states at time \eqn{t} and \eqn{\tilde{x}_{t-1}^T}{xt1T} is the Kalman smoother estimate of the states at time \eqn{t-1}.
The estimated state residuals are returned in \code{state.residuals} and rows n+1 to \eqn{n+m} of \code{residuals}. There are no NAs in the estimated state residuals as an estimate of the state exists whether or not there are associated data.

\code{res1} and \code{res2} in the code below will be the same.
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
Note that the state residual at the last time step (\eqn{T}) will be NA because it is the residual associated with \eqn{\mathbf{x}_T}{x_T} to \eqn{\mathbf{x}_{T+1}}{x_{T+1}} and \eqn{T+1} is beyond the data.  Similarly, the variance matrix at the last time step will have NAs for the same reason.

\strong{Variance of the residuals}

In a state-space model, \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} are stochastic, and the model and state residuals are random variables \eqn{\hat{\mathbf{V}}_t}{hatV_t} and \eqn{\hat{\mathbf{W}}_{t+1}}{hatW_{t+1}}. To evaluate the residuals we observed (with \eqn{\mathbf{y}^{(1)}}), we use the joint distribution of \eqn{\hat{\mathbf{V}}_t, \hat{\mathbf{W}}_{t+1}}{hatV_t, hatW_{t+1}} across all the different possible data sets that our MARSS equations with parameters \eqn{\Theta}{Theta} might generate. Denote the matrix of \eqn{\hat{\mathbf{V}}_t, \hat{\mathbf{W}}_{t+1}}{hatV_t, hatW_{t+1}}, as as \eqn{\widehat{\mathcal{E}}_t}{Epsilon_t}. That distribution has a mean and variance:
\deqn{ \textrm{E}[\widehat{\mathcal{E}}_t] = 0; \textrm{var}[\widehat{\mathcal{E}}_t] = \hat{\Sigma}_t }{ E[Epsilon_t] = 0; var[Epsilon_t] = hatSigma_t}
Our observed residuals \code{residuals} are one sample from this distribution.
To standardize the observed residuals, we will use \eqn{ \hat{\Sigma}_t }{ hatSigma_t }. The joint distribution of the residuals is returned in \code{var.residuals}. Rows/columns 1 to n are the conditional variances of the model residuals and rows/columns n+1 to n+m are the conditional variances of the state residuals.

\strong{Standardized residuals}

\code{residuals.marssMLE} will return the Cholesky standardized residuals sensu Harvey et al. (1998) in \code{std.residuals} for outlier and shock detection.  These are the model and state residuals multiplied by the inverse of the Cholesky decomposition of \code{var.residuals}. The standardized model residuals are set to NA when there are missing data. The standardized state residuals however always exist since the expected value of the states exist without data. The calculation of the standardized residuals for both the observations and states requires the full residuals variance matrix. Since the state residuals variance is NA at the last time step, the standarized residual in the last time step will be all NA.

The interpretation of the Cholesky standardized residuals is not straight-forward when the \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in MVN(0,I) space.  For example, if v is 2x2 correlated errors with variance-covariance matrix R. The transformed residuals (from this function) for the i-th row of v is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the transformed residuals can look rather non-intuitive.

The marginal standardized residuals are returned in \code{mar.residuals}. These are the model and state residuals multiplied by the inverse of the diagonal matrix formed by the square root of the diagonal of \code{var.residuals}. These residuals will be correlated (across the residuals at time \eqn{t}) but are easier to interpret when \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} are non-diagonal.


\strong{Normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{W_t} and \eqn{V_t} are \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} and the model is assumed to be written as
\deqn{ y_t = Z x_t + a + v_t}
\deqn{ x_t = B x_{t-1} + u + w_t}
Harvey et al (1998) writes the model as
\deqn{ y_t = Z x_t + a + Hv_t}
\deqn{ x_t = B x_{t-1} + u + Gw_t}
with the variance of \eqn{V_t} and \eqn{W_t} equal to I (identity).

\code{MARSSresiduals.tT} returns the residuals defined as in the first equations. To get the residuals defined as Harvey et al. (1998) define them (second equations), then use \code{normalize=TRUE}.  In that case the unconditional variance of residuals will be I instead of \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R}.

\strong{Missing or left-out data}

\eqn{ \textrm{E}[\widehat{\mathcal{E}}_t] }{ E[Epsilon_t]} and \eqn{ \textrm{var}[\widehat{\mathcal{E}}_t] = \hat{\Sigma}_t }{ var[Epsilon_t] = hatSigma_t} are for the distribution across all possible  \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y}. We can also compute the expected value and variance conditioned on a specific value of \eqn{\mathbf{Y}}{Y}, the one we obseved \eqn{\mathbf{y^{(1)}}}{y(1)} (Holmes 2014). If there are no missing values, this is not very interesting as  \eqn{\textrm{E}\[\hat{\mathbf{V}}_t|\mathbf{y}^{(1)}]=\hat{\mathbf{v}}_t}{E[hatV_t|y(1)] = hatv_t} and \eqn{\textrm{var}\[\hat{\mathbf{V}}_t|\mathbf{y}^{(1)}] = 0}{var[hatV_t|y(1)] = 0}. The variance is zero because \eqn{ \hat{\mathbf{V}}_t|\mathbf{y}^{(1)} }{ hatV_t|y(1) } is a fixed value. If we have data that are missing because we left them out, however, \eqn{\textrm{E}\[\hat{\mathbf{V}}_t|\mathbf{y}^{(1)}]}{E[hatV_t|y(1)]} and \eqn{\textrm{var}\[\hat{\mathbf{V}}_t|\mathbf{y}^{(1)}]}{var[hatV_t|y(1)]} are the values we need to evaluate whether the left-out data are unusual relative to what you expect given the data you did collect. 

\code{E.obs.residuals} is the conditional expected value of \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV_t|y(1)}. It is 
\deqn{\textrm{E}[\mathbf{Y}_t|\mathbf{y}^{(1)}] - \mathbf{Z}\tilde{x}_t^T - \mathbf{a} }{ E[Y_t|y(1)] - Z xtT - a}
It is similar to \eqn{\hat{\mathbf{v}}_t}{hatv_t}. The difference is the \eqn{\mathbf{y}}{y} term. \eqn{\textrm{E}[\mathbf{Y}^{(1)}_t|\mathbf{y}^{(1)}] }{ E[Y(1)_t|y(1)] } is \eqn{\mathbf{y}^{(1)}_t}{y(1)_t} for the non-missing values. For the missing values, the value depends on \eqn{\mathbf{R}}{R}. If \eqn{\mathbf{R}}{R} is diagonal, \eqn{\textrm{E}[\mathbf{Y}^{(2)}_t|\mathbf{y}^{(1)}] }{ E[Y(2)_t|y(1)] } is \eqn{\mathbf{Z}\tilde{x}_t^T + \mathbf{a}}{Z xtT + a} and the expected residual value is 0. If \eqn{\mathbf{R}}{R} is non-diagonal however, it will be non-zero.

\code{var.obs.residuals} is the conditional variance of \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV_t|y(1)} (eqn 24 in Holmes (2014)).  For the non-missing values, this variance is 0 since \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV_t|y(1)} is a fixed value. For the missing values, \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV_t|y(1)} is not fixed because \eqn{\mathbf{Y}^{(2)}}{Y(2)} is a random variable. For these values, the variance of \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV_t|y(1)} is determined by the variance of \eqn{\mathbf{Y}^{(2)}}{Y(2)} conditioned on \eqn{\mathbf{Y}^{(1)}=\mathbf{y}^{(1)}}{Y(1)=y(1)}. This variance matrix is returned in \code{var.obs.residuals}. The variance of \eqn{\hat{\mathbf{W}}|\mathbf{y}^{(1)}}{hatW_t|y(1)} is 0 and thus is not included.

The variance of \eqn{\textrm{var}[\hat{\mathbf{V}}_t|\mathbf{Y}^{(1)}] }{ var[hatV_t|Y(1)] } (notice uppercase \eqn{ \mathbf{Y} }{Y}) returned in the 1 to n rows/columns of \code{var.residuals} may also be of interest depending on what your are investigating with regards to missing values. For example, it may be of interest in a simulation study or cases where you have multiple replicated \eqn{\mathbf{Y}}{Y} data sets. \code{var.residuals} would allow you to determine if the left-out residuals are unusual with regards to what you would expect for left-out data in that location of the \eqn{\mathbf{Y}}{Y} matrix but not relative to the data you did collect. If \eqn{\mathbf{R}}{R} is non-diagonal and the \eqn{\mathbf{y}^{(1)}}{y(1)} and \eqn{\mathbf{y}^{(2)}}{y(2)} are highly correlated, the variance of \eqn{\textrm{var}[\hat{\mathbf{V}}_t|\mathbf{Y}^{(1)}] }{ var[hatV_t|Y(1)] } and variance of \eqn{\textrm{var}[\hat{\mathbf{V}}_t|\mathbf{y}^{(1)}] }{ var[hatV_t|y(1)] } would be quite different. In the latter, the variance is low because \eqn{\mathbf{y}^{(1)} }{ y(1) } has strong information about \eqn{\mathbf{y}^{(2)} }{ y(2) }. In the former, we integrate over \eqn{\mathbf{Y}^{(1)} }{ Y(1) } and the variance could be high (depending on the parameters).
}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{residuals.marssMLE}}, \code{\link{MARSSresiduals.tt1}}, \code{\link{fitted.marssMLE}}, \code{\link{plot.marssMLE}} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  MLEobj <- MARSS(dat)
  
  #state residuals
  state.resids1 <- residuals(MLEobj, conditioning="T")$state.residuals
  #this is the same as
  states <- MLEobj$states
  Q <- coef(MLEobj,type="matrix")$Q
  state.resids2 <- states[,2:30]-states[,1:29]-matrix(coef(MLEobj,type="matrix")$U,2,29)
  #compare the two
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #normalize to variance of 1
  state.resids1 <- residuals(MLEobj, normalize=TRUE, conditioning="T")$state.residuals
  state.resids2 <- (solve(t(chol(Q))) \%*\% state.resids2)
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #Cholesky standardized (by joint variance) model & state residuals
  residuals(MLEobj)$std.residuals
}
\references{
Harvey, A., S. J. Koopman, and J. Penzer. 1998. Messy time series: a unified approach. Advances in Econometrics 13: 103-144  (see page 112-113).  Eqn 21 is the Kalman eqns.  Eqn 23 and 24 is the backward recursion to compute the smoothations.  This function uses the MARSSkf output for eqn 21 and then implements the backwards recursion in eqn 23 and eqn 24.  Pages 120-134 discuss the use of standardized residuals for outlier and structural break detection.

de Jong, P. and J. Penzer. 1998. Diagnosing shocks in time series. Journal of the American Statistical Association 93: 796-806.  This one shows the same equations; see eqn 6.  This paper mentions the scaling based on the inverse of the sqrt (chol) of the variance-covariance matrix for the residuals (model and state together).  This is in the right column, half-way down on page 800.

Koopman, S. J., N. Shephard, and J. A. Doornik. 1999. Statistical algorithms for models in state space using SsfPack 2.2. Econometrics Journal 2: 113-166. (see pages 147-148).

Harvey, A. and S. J. Koopman. 1992. Diagnostic checking of unobserved-components time series models. Journal of Business & Economic Statistics 4: 377-389.

Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 
}

