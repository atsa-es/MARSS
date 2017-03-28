\name{predict.marssMLE}
\alias{predict}
\alias{predict.MARSS}
\alias{predict.marssMLE}
\title{ Compute the Prediction Intervals, Expected Values, and Standard Errors for States (X) and Observation (Y) from MARSS fits }
\description{
This function is still in beta development.  Not exported in version 3.5.  Computes the expected values of states (X component) and observations (Y component) for MARSS models fit with \code{\link{MARSS}}. A list of \code{newdata}, which includes the inputs (e.g. c and d in a marxss model) must be passed in if they are included in the model.  Optionally data (Y) can be passed in in \code{newdata}.  This can be the original data or new data and can have missing values.  To simulate data, use \code{\link{simulate.marssMLE}}.
}
\usage{
\method{predict}{marssMLE}(object, ..., 
                           n.ahead=1, 
                           t.start=NULL, 
                           newdata=list(), 
                           se.fit=TRUE, 
                           nboot=1000, 
                           param.gen="hessian", 
                           verbose=FALSE, 
                           prediction.intervals=TRUE)
}
\arguments{
  \item{object}{ A marssMLE object \code{\link{marssMLE}}. }
  \item{...}{ Not used. }
  \item{n.ahead}{ Number of time steps ahead to predict (starting from \code{t.start}).}
  \item{t.start}{ A what t to start the prediction relative to the start (t=1) of the data to which the model was fit.  The default is to start at t=T+1, where T is the length of the data to which the model was fit. For example, if TT=10, t.start=1 and n.ahead = 5, prediction is for time steps 1, 2, 3, 4, 5 so is not past the end of the data.}
  \item{newdata}{ Data (explanatory variables) to be used for prediction if included in the model.}
  \item{se.fit}{ Logical: should standard errors of prediction be returned? }
    \item{prediction.intervals}{ Logical. Whether to return prediction intervals.  Calculation of prediction intervals can be slow (if \code{param.gen="MLE"}) and will exit with an error if the Hessian cannot be computed (if \code{param.gen="hessian"}).}
      \item{param.gen}{ "MLE" or "hessian". Default="hessian". The perdiction interval calculation requires samples of the parameters from their estimated distribution.  You can generate these via the estimated Hessian matrix (fast but assumes normality) or via parametric bootstrapped (terribly slow). }
      \item{nboot}{ Number of bootstrapped parameter sequences to generate for the prediction interval calculation. }
      \item{verbose}{ If TRUE, the boot straps used for the prediction interval calculation are returned. }
}
\details{
  \code{newdata} needs to be consistent with the model fitted.  Thus if the model included inputs (in X or Y component) then these must also be included in \code{newdata} and have the proper size and form.  If \code{newdata$data} is not passed in, then data (Y component) are assumed to be all missing and the expected values of the data (Ys) will be returned.  The data is assumed to be all missing even if \code{t.start} is less than T and thus the prediction period overlaps the data to which the model was fit.  If that data should be included in the predictions, then it must be passed in in \code{newdata}.  The number of time steps in components of \code{newdata} must be consistent with \code{n.ahead}.
  
  Prediction of the state (X) component is via the Kalman smoother (i.e. conditioned on all the data, if there is any) using \code{\link{MARSSkf}}.  Prediction of the observation (Y) component is via  \code{\link{MARSShatyt}}, and is also smoothed (i.e. conditioned on all the data, if there is any).  The standard errors of prediction exclude the uncertainty in the estimation of the model parameters. You could use bootstrapping via \code{\link{MARSSboot}} if you wanted to include that.
}
\value{
  \item{E.x}{ A matrix of the expected values of the states (X) given \code{newdata}. This is the Kalman smoother \code{xtT} output using \code{newdata} as the input for data and any inputs. See \code{\link{MARSSkf}}. }
  \item{E.y}{  A matrix of the expected values of the observations (Y) given \code{newdata}. This is \code{ytT} output using \code{newdata} as the input for data and any inputs. See \code{\link{MARSShatyt}}.}
  \item{x.se}{ A matrix of the estimated standard errors of the X predictions. }
  \item{y.se}{ A matrix of the estimated standard errors of the Y predictions. }
  \item{newdata}{ The \code{newdata} list that was constructed from the inputs. }
  \item{pred.MLEobj}{ The \code{\link{marssMLE}} object that was constructed with \code{newdata}.  This can be passed to \code{\link{MARSSkf}} if the full Kalman filter and smoother output is needed. }
  \item{pred.MLEobj}{ The \code{\link{marssMLE}} object that was constructed with \code{newdata}.  This can be passed to \code{\link{MARSSkf}} if the full Kalman filter and smoother output is needed. }
  \item{pred.MLEobj}{ The \code{\link{marssMLE}} object that was constructed with \code{newdata}.  This can be passed to \code{\link{MARSSkf}} if the full Kalman filter and smoother output is needed. }
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{marssMODEL}} \code{\link{marssMLE}} \code{\link{MARSSboot}} 
}


