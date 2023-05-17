\name{checkMARSSInputs}
\alias{checkMARSSInputs}
\title{ Check inputs to MARSS call }
\keyword{internal}

\description{
  This is a helper function to check the inputs to a \code{\link{MARSS}()} call for any errors. Not exported.
}
\usage{
checkMARSSInputs( MARSS.inputs, silent = FALSE )
}
\arguments{
  \item{MARSS.inputs}{ A list comprised of the needed inputs to a MARSS call: 
    data, 
    inits, 
    model,
    control,
    method,
    form) }
  \item{silent}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }
}
\details{
This is a helper function to check that all the inputs to a \code{\link{MARSS}()} function call are properly specified.  

If arguments \code{inits} or \code{control} are not provided by the user, they will be set by the \code{alldefaults[[method]]} object specified in \code{MARSSsettings}.  Argument \code{model} specifies the model structure using a list of matrices; see \code{\link{MARSS}} or the \href{https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf}{User Guide} for instructions on how to specify model structure. If \code{model} is left off, then the function \code{MARSS.form()} is used to determine the default model structure.
}
\value{
  If the function does not stop due to errors, it returns an updated list with elements
  \item{data}{ Data supplied by user. }
  \item{model}{ Not changed.  Will be updated by the \code{MARSS.form} function (e.g. \code{\link{MARSS.marxss}}). }
  \item{inits}{ A list specifying initial values for parameters to be used at iteration 1 in iterative maximum-likelihood algorithms. }
  \item{method}{ The method used for estimation. }
  \item{form}{ The equation form used to convert wrapper object to a \code{\link{marssMODEL}} object. }
  \item{control}{ See Arguments. }
}

\author{ 
  Kellie Wills, NOAA, Seattle, USA.
}
\seealso{ 
  \code{\link{MARSS}()}, \code{\link{marssMODEL}}, \code{\link{checkModelList}()} 
}


