\name{MARSSresiduals.tT}
\alias{MARSSresiduals.tT}

\title{ MARSS Smoothed Residuals }
\description{
  Calculates the standardized (or auxiliary) smoothed residuals sensu Harvey, Koopman and Penzer (1998). The expected values and variance for missing (or left-out) data are also returned (Holmes 2014). Not exported. Access this function with \code{MARSSresiduals(object, type="tT")}. At time \eqn{t} (in the returned matrices), the model residuals are for time \eqn{t}, while the state residuals are for the transition from \eqn{t} to \eqn{t+1} following the convention in Harvey, Koopman and Penzer (1998).
}
\usage{
MARSSresiduals.tT(object, Harvey = FALSE, normalize = FALSE, 
    silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss"))
}
\arguments{
  \item{object}{ An object of class \code{\link{marssMLE}}.}
  \item{Harvey}{ TRUE/FALSE. Use the Harvey et al. (1998) algorithm or use the Holmes (2014) algorithm. The values are the same except for missing values. }
  \item{normalize}{ TRUE/FALSE See details.}
  \item{silent}{ If TRUE, don't print inversion warnings. }
  \item{fun.kf}{ Kalman filter function to use. Can be ignored. }
}
\value{
A list with the following components  
  \item{model.residuals}{ The the observed smoothed model residuals: data minus the model predictions conditioned on all observed data. This is different than the Kalman filter innovations which use on the data up to time \eqn{t-1} for the predictions. See details. }
  \item{state.residuals}{ The smoothed state residuals \eqn{\mathbf{x}_{t+1}^T - \mathbf{Z} \mathbf{x}_{t}^T - \mathbf{u}}{xtT(t+1) - Z xtT(t) - u}. The last time step will be NA because the last step would be for T to T+1 (past the end of the data).}
  \item{residuals}{ The residuals conditioned on the observed data. Returned as a (n+m) x T matrix with \code{model.residuals} in rows 1 to n and \code{state.residuals} in rows n+1 to n+m.  NAs will appear in rows 1 to n in the places where data are missing. }
  \item{var.residuals}{ The joint variance of the model and state residuals conditioned on observed data. Returned as a (n+m) x (n+m) x T matrix. For Harvey=FALSE, this is Holmes (2014) equation 57. For Harvey=TRUE, this is the residual variance in eqn. 24, page 113, in Harvey et al. (1998). They are identical except for missing values, for those Harvey=TRUE returns 0s. For the state residual variance, the last time step will be all NA because the last step would be for T to T+1 (past the end of the data).}
  \item{std.residuals}{ The Cholesky standardized residuals as a (n+m) x T matrix. This is \code{residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals}. The model standardized residuals associated with the missing data are replaced with NA. }
  \item{mar.residuals}{ The marginal standardized residuals as a (n+m) x T matrix. This is \code{residuals} multiplied by the inverse of the diagonal matrix formed by the square-root of the diagonal of \code{var.residuals}.  The model marginal residuals associated with the missing data are replaced with NA. }
  \item{bchol.residuals}{ The Block Cholesky standardized residuals as a (n+m) x T matrix. This is \code{model.residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals[1:n,1:n,]} and \code{state.residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals[(n+1):(n+m),(n+1):(n+m),]}. }
  \item{E.obs.residuals}{ The expected value of the model residuals conditioned on the observed data. Returned as a n x T matrix.  For observed data, this will be the observed residuals (values in \code{model.residuals}). For unobserved data, this will be 0 if \eqn{\mathbf{R}}{R} is diagonal but non-zero if \eqn{\mathbf{R}}{R} is non-diagonal. See details. }
    \item{var.obs.residuals}{ The variance of the model residuals conditioned on the observed data. Returned as a n x n x T matrix.  For observed data, this will be 0. See details. }
  \item{msg}{ Any warning messages. This will be printed unless Object$control$trace = -1 (suppress all error messages). }

}
\details{

This function returns the raw, the Cholesky standardized and the marginal standardized smoothed model and state residuals.  'smoothed' means conditioned on all the observed data and a set of parameters. These are the residuals presented in Harvey, Koopman and Penzer (1998) pages 112-113, with the addition of the values for unobserved data (Holmes 2014).  If Harvey=TRUE, the function uses the algorithm on page 112 of Harvey, Koopman and Penzer (1998) to compute the conditional residuals and variance of the residuals.  If Harvey=FALSE, the function uses the equations in the technical report (Holmes 2014).  Unlike the innovations residuals, the smoothed residuals are autocorrelated (section 4.1 in Harvey and Koopman 1992) and thus an ACF test on these residuals would not reveal model inadequacy.

The residuals matrix has a value for each time step.  The residuals in column \eqn{t} rows 1 to n are the model residuals associated with the data at time \eqn{t}. The residuals in rows n+1 to n+m are the state residuals associated with the transition from \eqn{\mathbf{x}_{t}}{x(t)} to \eqn{\mathbf{x}_{t+1}}{x(t+1)}, not the transition from \eqn{\mathbf{x}_{t-1}}{x(t-1)} to \eqn{\mathbf{x}_{t}}{x(t)}. Because \eqn{\mathbf{x}_{t+1}}{x(t+1)} does not exist at time \eqn{T}, the state residuals and associated variances at time \eqn{T} are NA.

Below the conditional residuals and their variance are discussed. The random variables are capitalized and the realizations from the random variables are lower case. The random variables are \eqn{\mathbf{X}}{X}, \eqn{\mathbf{Y}}{Y}, \eqn{\mathbf{V}}{V} and \eqn{\mathbf{W}}{W}. There are two types of \eqn{\mathbf{Y}}{Y}. The observed \eqn{\mathbf{Y}}{Y} that are used to estimate the states \eqn{\mathbf{x}}{x}. These are termed \eqn{\mathbf{Y}^{(1)}}{Y(1)}. The unobserved  \eqn{\mathbf{Y}}{Y} are termed \eqn{\mathbf{Y}^{(2)}}{Y(2)}. These are not used to estimate the states \eqn{\mathbf{x}}{x} and we may or may not know the values of \eqn{\mathbf{y}^{(2)}}{y(2)}. Typically we treat \eqn{\mathbf{y}^{(2)}}{y(2)} as unknown but it may be known but we did not include it in our model fitting.  Note that the model parameters \eqn{\Theta}{Theta} are treated as fixed or known. The 'fitting' does not involve estimating \eqn{\Theta}{Theta}; it involves estimating \eqn{\mathbf{x}}{x}. All MARSS parameters can be time varying but the \eqn{t} subscripts are left off parameters to reduce clutter.

\strong{Model residuals}

\eqn{\mathbf{v}_{t}}{v(t)} is the difference between the data and the predicted data at time \eqn{t} given \eqn{\mathbf{x}_{t}}{x(t)}:
\deqn{ \mathbf{v}_{t} = \mathbf{y}_{t} - \mathbf{Z} \mathbf{x}_{t} - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ v(t) = y(t) - Z x(t) - a - D d(t)}
\eqn{\mathbf{x}_{t}}{x(t)} is unknown (hidden) and our data are one realization of \eqn{\mathbf{y}_{t}}{y(t)}. The observed model residuals \eqn{\hat{\mathbf{v}}_{t}}{hatv(t)} are the difference between the observed data and the predicted data at time \eqn{t} using the fitted model. \code{MARSSresiduals.tT} fits the model using all the data, thus
\deqn{ \hat{\mathbf{v}}_{t} = \mathbf{y}_{t} - \mathbf{Z}\mathbf{x}_{t}^T - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ hatv(t) = y(t) - Z xtT(t) - a - D d(t)}
where \eqn{\mathbf{x}_{t}^T}{xtT(t)} is the expected value of \eqn{\mathbf{X}_{t}}{X(t)} conditioned on the data from 1 to \eqn{T} (all the data), i.e. the Kalman smoother estimate of the states at time \eqn{t}. \eqn{\mathbf{y}_{t}}{y(t)} are your data and missing values will appear as NA in the observed model residuals. These are returned as \code{model.residuals} and rows 1 to \eqn{n} of \code{residuals}.

\code{res1} and \code{res2} in the code below will be the same.
\preformatted{dat = t(harborSeal)[2:3,]
fit = MARSS(dat)
Z = coef(fit, type="matrix")$Z
A = coef(fit, type="matrix")$A
res1 = dat - Z \%*\% fit$states - A \%*\% matrix(1,1,ncol(dat))
res2 = MARSSresiduals(fit, type="tT")$model.residuals
}

\strong{State residuals}

\eqn{\mathbf{w}_{t+1}}{w(t+1)} are the difference between the state at time \eqn{t+1} and the expected value of the state at time \eqn{t+1} given the state at time \eqn{t}:
\deqn{ \mathbf{w}_{t+1} = \mathbf{x}_{t+1} - \mathbf{B} \mathbf{x}_{t} - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ w(t+1) = x(t+1) - B x(t) - u - C c(t+1)}
The estimated state residuals \eqn{\hat{\mathbf{w}}_{t+1}}{hatw(t+1)} are the difference between estimate of \eqn{\mathbf{x}_{t+1}}{x(t+1)} minus the estimate using \eqn{\mathbf{x}_{t}}{x(t)}. 
\deqn{ \hat{\mathbf{w}}_{t+1} = \mathbf{x}_{t+1}^T - \mathbf{B}\mathbf{x}_{t}^T - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ hatw(t+1) = xtT(t+1) - B xtT(t) - u - C c(t+1)}
where \eqn{\mathbf{x}_{t+1}^T}{xtT(t+1)} is the Kalman smoother estimate of the states at time \eqn{t+1} and \eqn{\mathbf{x}_{t}^T}{xtT(t)} is the Kalman smoother estimate of the states at time \eqn{t}.
The estimated state residuals \eqn{\mathbf{w}_{t+1}}{w(t+1)} are returned in \code{state.residuals} and rows \eqn{n+1} to \eqn{n+m} of \code{residuals}. \code{state.residuals[,t]} is \eqn{\mathbf{w}_{t+1}}{w(t+1)} (notice time subscript difference). There are no NAs in the estimated state residuals as an estimate of the state exists whether or not there are associated data.

\code{res1} and \code{res2} in the code below will be the same.
\preformatted{dat <- t(harborSeal)[2:3,]
TT <- ncol(dat)
fit <- MARSS(dat)
B <- coef(fit, type="matrix")$B
U <- coef(fit, type="matrix")$U
statestp1 <- MARSSkf(fit)$xtT[,2:TT]
statest <- MARSSkf(fit)$xtT[,1:(TT-1)]
res1 <- statestp1 - B \%*\% statest - U \%*\% matrix(1,1,TT-1)
res2 <- MARSSresiduals(fit, type="tT")$state.residuals[,1:(TT-1)]
}
Note that the state residual at the last time step (not shown) will be NA because it is the residual associated with \eqn{\mathbf{x}_T}{x_T} to \eqn{\mathbf{x}_{T+1}}{x_{T+1}} and \eqn{T+1} is beyond the data.  Similarly, the variance matrix at the last time step will have NAs for the same reason.

\strong{Variance of the residuals}

In a state-space model, \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} are stochastic, and the model and state residuals are random variables \eqn{\hat{\mathbf{V}}_{t}}{hatV(t)} and \eqn{\hat{\mathbf{W}}_{t+1}}{hatW(t+1)}. To evaluate the residuals we observed (with \eqn{\mathbf{y}^{(1)}}{y(1)}), we use the joint distribution of \eqn{\hat{\mathbf{V}}_{t}, \hat{\mathbf{W}}_{t+1}}{hatV(t), hatW(t+1)} across all the different possible data sets that our MARSS equations with parameters \eqn{\Theta}{Theta} might generate. Denote the matrix of \eqn{\hat{\mathbf{V}}_{t}, \hat{\mathbf{W}}_{t+1}}{hatV(t), hatW(t+1)}, as \eqn{\widehat{\mathcal{E}}_{t}}{Epsilon(t)}. That distribution has an expected value (mean) and variance:
\deqn{ \textrm{E}[\widehat{\mathcal{E}}_{t}] = 0; \textrm{var}[\widehat{\mathcal{E}}_{t}] = \hat{\Sigma}_{t} }{ E[Epsilon(t)] = 0; var[Epsilon(t)] = hatSigma(t)}
Our observed residuals (returned in \code{residuals}) are one sample from this distribution.
To standardize the observed residuals, we will use \eqn{ \hat{\Sigma}_{t} }{ hatSigma(t) }. \eqn{ \hat{\Sigma}_{t} }{ hatSigma(t) } is returned in \code{var.residuals}. Rows/columns 1 to \eqn{n} are the conditional variances of the model residuals and rows/columns \eqn{n+1} to \eqn{n+m} are the conditional variances of the state residuals. The off-diagonal blocks are the covariances between the two types of residuals.

\strong{Standardized residuals}

\code{MARSSresiduals} will return the Cholesky standardized residuals sensu Harvey et al. (1998) in \code{std.residuals} for outlier and shock detection.  These are the model and state residuals multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals} (note \code{chol()} in R returns the upper triangle thus a transpose is needed). The standardized model residuals are set to NA when there are missing data. The standardized state residuals however always exist since the expected value of the states exist without data. The calculation of the standardized residuals for both the observations and states requires the full residuals variance matrix. Since the state residuals variance is NA at the last time step, the standardized residual in the last time step will be all NA (for both model and state residuals).

The interpretation of the Cholesky standardized residuals is not straight-forward when the \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in MVN(0,I) space.  For example, if v is 2x2 correlated errors with variance-covariance matrix R. The transformed residuals (from this function) for the i-th row of \eqn{\mathbf{v}}{v} is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the transformed residuals can look rather non-intuitive.

The marginal standardized residuals are returned in \code{mar.residuals}. These are the model and state residuals multiplied by the inverse of the diagonal matrix formed by the square root of the diagonal of \code{var.residuals}. These residuals will be correlated (across the residuals at time \eqn{t}) but are easier to interpret when \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} are non-diagonal.

The Block Cholesky standardized residuals are like the Cholesky standardized residuals except that the full variance-covariance matrix is not used, only the variance-covariance matrix for the model or state residuals (respectively) is used for standardization. For the model residuals, the Block Cholesky standardized residuals will be the same as the Cholesky standardized residuals because the upper triangle of the lower triangle of the Cholesky decomposition (which is what we standardize by) is all zero. For the state residuals, the Block Cholesky standardization will be different because Block Cholesky standardization treats the model and state residuals as independent (which they are not in the smoothations case). 


\strong{Normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} are \eqn{\mathbf{R}}{R} and \eqn{\mathbf{Q}}{Q} and the model is assumed to be written as
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{v}_t}{ y(t) = Z x(t) + a + v(t)}
\deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{w}_t}{ x(t) = B x(t-1) + u + w(t)}
If normalize=TRUE, the model is assumed to be written
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{H}\mathbf{v}_t}{ y(t) = Z x(t) + a + Hv(t)}
\deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{G}\mathbf{w}_t}{ x(t) = B x(t-1) + u + Gw(t)}
with the variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} equal to \eqn{\mathbf{I}}{I} (identity).

\code{MARSSresiduals.tT} returns the residuals defined as in the first equations. To get the residuals defined as Harvey et al. (1998) define them (second equations), then use \code{normalize=TRUE}.  In that case the unconditional variance of residuals will be \eqn{\mathbf{I}}{I} instead of \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R}.

\strong{Missing or left-out data}

\eqn{ \textrm{E}[\widehat{\mathcal{E}}_{t}] }{ E[Epsilon(t)] } and \eqn{ \textrm{var}[\widehat{\mathcal{E}}_{t}] }{ var[Epsilon(t)] } are for the distribution across all possible  \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y}. We can also compute the expected value and variance conditioned on a specific value of \eqn{\mathbf{Y}}{Y}, the one we observed \eqn{\mathbf{y}^{(1)}}{y(1)} (Holmes 2014). If there are no missing values, this is not very interesting as  \eqn{\textrm{E}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}]=\hat{\mathbf{v}}_{t}}{E[hatV(t)|y(1)] = hatv(t)} and \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}] = 0}{var[hatV(t)|y(1)] = 0}. If we have data that are missing because we left them out, however, \eqn{\textrm{E}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}]}{E[hatV(t)|y(1)]} and \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}]}{var[hatV(t)|y(1)]} are the values we need to evaluate whether the left-out data are unusual relative to what you expect given the data you did collect. 

\code{E.obs.residuals} is the conditional expected value \eqn{\textrm{E}[\hat{\mathbf{V}}|\mathbf{y}^{(1)}]}{E[hatV(t)|y(1)]} (notice small \eqn{\mathbf{y}}{y}). It is 
\deqn{\textrm{E}[\mathbf{Y}_{t}|\mathbf{y}^{(1)}] - \mathbf{Z}\mathbf{x}_t^T - \mathbf{a} }{ E[Y(t)|y(1)] - Z xtT(t) - a}
It is similar to \eqn{\hat{\mathbf{v}}_{t}}{hatv(t)}. The difference is the \eqn{\mathbf{y}}{y} term. \eqn{\textrm{E}[\mathbf{Y}^{(1)}_{t}|\mathbf{y}^{(1)}] }{ E[Y(1)(t)|y(1)] } is \eqn{\mathbf{y}^{(1)}_{t}}{y(1)(t)} for the non-missing values. For the missing values, the value depends on \eqn{\mathbf{R}}{R}. If \eqn{\mathbf{R}}{R} is diagonal, \eqn{\textrm{E}[\mathbf{Y}^{(2)}_{t}|\mathbf{y}^{(1)}] }{ E[Y(2)(t)|y(1)] } is \eqn{\mathbf{Z}\mathbf{x}_t^T + \mathbf{a}}{Z xtT(t) + a} and the expected residual value is 0. If \eqn{\mathbf{R}}{R} is non-diagonal however, it will be non-zero.

\code{var.obs.residuals} is the conditional variance  \eqn{\textrm{var}[\hat{\mathbf{V}}|\mathbf{y}^{(1)}]}{var[hatV(t)|y(1)]} (eqn 24 in Holmes (2014)).  For the non-missing values, this variance is 0 since \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV(t)|y(1)} is a fixed value. For the missing values, \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV(t)|y(1)} is not fixed because \eqn{\mathbf{Y}^{(2)}}{Y(2)} is a random variable. For these values, the variance of \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV(t)|y(1)} is determined by the variance of \eqn{\mathbf{Y}^{(2)}}{Y(2)} conditioned on \eqn{\mathbf{Y}^{(1)}=\mathbf{y}^{(1)}}{Y(1)=y(1)}. This variance matrix is returned in \code{var.obs.residuals}. The variance of \eqn{\hat{\mathbf{W}}|\mathbf{y}^{(1)}}{hatW(t)|y(1)} is 0 and thus is not included.

The variance \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{Y}^{(1)}] }{ var[hatV(t)|Y(1)] } (uppercase \eqn{ \mathbf{Y} }{Y}) returned in the 1 to \eqn{n} rows/columns of \code{var.residuals} may also be of interest depending on what you are investigating with regards to missing values. For example, it may be of interest in a simulation study or cases where you have multiple replicated \eqn{\mathbf{Y}}{Y} data sets. \code{var.residuals} would allow you to determine if the left-out residuals are unusual with regards to what you would expect for left-out data in that location of the \eqn{\mathbf{Y}}{Y} matrix but not specifically relative to the data you did collect. If \eqn{\mathbf{R}}{R} is non-diagonal and the \eqn{\mathbf{y}^{(1)}}{y(1)} and \eqn{\mathbf{y}^{(2)}}{y(2)} are highly correlated, the variance of \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{Y}^{(1)}] }{ var[hatV(t)|Y(1)] } and variance of \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}] }{ var[hatV(t)|y(1)] } for the left-out data would be quite different. In the latter, the variance is low because \eqn{\mathbf{y}^{(1)} }{ y(1) } has strong information about \eqn{\mathbf{y}^{(2)} }{ y(2) }. In the former, we integrate over \eqn{\mathbf{Y}^{(1)} }{ Y(1) } and the variance could be high (depending on the parameters).

Note, if \code{Harvey=TRUE} then the rows and columns of \code{var.residuals} corresponding to missing values will be NA. This is because the Harvey et al. algorithm does not compute the residual variance for missing values.
}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
} 
\seealso{ \code{\link{MARSSresiduals}}, \code{\link{MARSSresiduals.tt1}}, \code{\link{fitted.marssMLE}}, \code{\link{plot.marssMLE}} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  fit <- MARSS(dat)
  
  #state residuals
  state.resids1 <- MARSSresiduals(fit, type="tT")$state.residuals
  #this is the same as hatx_t-(hatx_{t-1}+u)
  states <- fit$states
  state.resids2 <- states[,2:30]-states[,1:29]-matrix(coef(fit,type="matrix")$U,2,29)
  #compare the two
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #normalize the state residuals to a variance of 1
  Q <- coef(fit,type="matrix")$Q
  state.resids1 <- MARSSresiduals(fit, type="tT", normalize=TRUE)$state.residuals
  state.resids2 <- (solve(t(chol(Q))) \%*\% state.resids2)
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #Cholesky standardized (by joint variance) model & state residuals
  MARSSresiduals(fit, type="tT")$std.residuals
  
  # Returns residuals in a data frame in long form
  residuals(fit, type="tT")
}
\references{
Harvey, A., S. J. Koopman, and J. Penzer. 1998. Messy time series: a unified approach. Advances in Econometrics 13: 103-144  (see page 112-113).  Equation 21 is the Kalman eqns.  Eqn 23 and 24 is the backward recursion to compute the smoothations.  This function uses the MARSSkf output for eqn 21 and then implements the backwards recursion in equation 23 and equation 24.  Pages 120-134 discuss the use of standardized residuals for outlier and structural break detection.

de Jong, P. and J. Penzer. 1998. Diagnosing shocks in time series. Journal of the American Statistical Association 93: 796-806.  This one shows the same equations; see eqn 6.  This paper mentions the scaling based on the inverse of the sqrt (Cholesky decomposition) of the variance-covariance matrix for the residuals (model and state together).  This is in the right column, half-way down on page 800.

Koopman, S. J., N. Shephard, and J. A. Doornik. 1999. Statistical algorithms for models in state space using SsfPack 2.2. Econometrics Journal 2: 113-166. (see pages 147-148).

Harvey, A. and S. J. Koopman. 1992. Diagnostic checking of unobserved-components time series models. Journal of Business & Economic Statistics 4: 377-389.

Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 
}

