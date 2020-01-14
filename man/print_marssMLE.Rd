\name{print.marssMLE}
\alias{print.marssMLE}

\title{ Printing functions for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs \code{\link{marssMLE}} objects.  \code{print(marssMLE)}, where \code{marssMLE} is a \code{\link{marssMLE}} object, will print out information on the fit.  However, \code{print} can be used to print a variety of information (residuals, smoothed states, imputed missing values, etc) from a \code{\link{marssMLE}} object using the \code{what} argument in the print call.
}
\usage{
\method{print}{marssMLE}(x, digits = max(3, getOption("digits")-4), ..., what="fit", form=NULL, silent=FALSE)
}
\arguments{
  \item{x}{ A marssMLE object.  }
  \item{digits}{ Number of digits for printing.  }
  \item{...}{ Other arguments for print. }
  \item{what}{ What to print.  Default is "fit".  If you input what as a vector, print returns a list. See examples.
  \itemize{
    \item{ "model"}{ The model parameters with names for the estimted parameters.  The output is customized by the form of the model that was fit.  This info is in \code{ attr(x$model, "form") }. }
    \item{ "par" }{ A list of only the estimated values in each matrix. Each model matrix has it's own list element. Standard function: \code{coef(x)}}
    \item{ "start" or "inits" }{ The values that the optimization algorithm was started at. Note, \code{x$start} shows this in form="marss" while \code{print} shows it in whatever form is in \code{ attr(x$model, "form") }.}
    \item{ "paramvector" }{ A vector of all the estimated values in each matrix. Standard function: \code{coef(x, type="vector")}.  See \code{\link{coef.marssMLE}}.}
    \item{ "par.se","par.bias","par.lowCIs","par.upCIs" }{ A vector the estimated parameter standard errors, parameter bias, lower and upper confidence intervals. Standard function: \code{MARSSparamCIs(x)}  See \code{\link{MARSSparamCIs}}.}
    \item{ "xtT" or "states" }{ The estimated states conditioned on all the data. \code{x$states}}
    \item{ "data" }{ The data. This is in \code{x$model$data}}
    \item{ "logLik" }{ The log-likelihood. Standard function: \code{x$logLik}.  See \code{\link{MARSSkf}} for a discussion of the computation of the log-likelihood for MARSS models. }
    \item{ "ytT" }{ The expected value of the data conditioned on all the data.  Returns the data if present and the expected value if missing. This is in \code{x$ytT} (ytT is analogous to xtT).}
    \item{ "states.se" }{ The state standard errors. \code{x$states.se} }
    \item{ "states.cis" }{ Approximate confidence intervals for the states. See \code{\link{MARSSparamCIs}}.}
    \item{ "model.residuals" }{ The smoothed model residuals. y(t)-E(y(t)|xtT(t)), aka actual data at time t minus the expected value of the data conditioned on the smoothed states estimate at time t. Standard function: \code{residuals(x)$model.residuals}}  See \code{\link{residuals.marssMLE}} for a discussion of residuals in the context of MARSS models.
    \item{ "state.residuals" }{ The smoothed state residuals. E(xtT(t))-E(x(t)|xtT(t-1)), aka the expected value of x at t conditioned on all the data minus the expected value of x at t conditioned on (x(t-1) conditioned on all the data). Standard function: \code{residuals(x)$state.residuals}}  See \code{\link{residuals.marssMLE}}.
    \item{ parameter name }{ Returns the parameter matrix for that parameter with fixed values at their fixed values and the estimated values at their estimated values. Standard function: \code{coef(x, type="matrix")$elem} }
    \item{ "kfs" }{ The Kalman filter and smoother output.  See \code{\link{MARSSkf}} for a description of the output. The full kf output is not normally attached to the output from a \code{\link{MARSS}()} call.  This will run the filter/smoother if needed and return the list INVISIBLY.  So assign the output as \code{foo=print(x,what="kfs")}}
    \item{ "Ey" }{ The expectations involving y conditioned on all the data.  See \code{\link{MARSShatyt}} for a discussion of these expectations.  This output is not normally attached to the output from a \code{\link{MARSS}()} call--except \code{ytT} which is the predicted value of any missing y. The list is returned INVISIBLY.  So assign the output as \code{foo=print(x,what="Ey")}.}
      } }
  \item{form}{ By default, print uses the model form specified in the call to \code{\link{MARSS}}.  This information is in \code{ attr(marssMLE$model, "form") }, however you can specify a different form.  \code{ form="marss" } should always work since this is the model form in which the model objects are stored (in \code{marssMLE$marss}).}
  \item{silent}{ If TRUE, do not print just return the object.  If print call is assigned, nothing will be printed.  See examples.  If \code{what="fit"}, there is always output printed.}
}

\value{
  A print out of information.  If you assign the print call to a value, then you can reference the output.  See the examples.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{ 
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  
  print(MLEobj)
  
  print(MLEobj,what="model")
  
  print(MLEobj,what="par")
  
  #silent doesn't mean silent unless the print output is assigned
  print(MLEobj,what="paramvector", silent=TRUE)
  tmp=print(MLEobj,what="paramvector", silent=TRUE)
  #silent means some info on what you are printing is shown whether
  #or not the print output is assigned
  print(MLEobj,what="paramvector", silent=FALSE)
  tmp=print(MLEobj,what="paramvector", silent=FALSE)
  
  cis=print(MLEobj,what="states.cis")
  cis$up95CI
  
  vars=print(MLEobj, what=c("R","Q"))
}