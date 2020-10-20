\name{marssMODEL-class}
\Rdversion{1.1}
\docType{class}
\alias{marssMODEL}
%\alias{print,marssMODEL-method}
%\alias{summary,marssMODEL-method}

\title{Class "marssMODEL"}
\description{
\code{marssMODEL} objects describe a vectorized form for the multivariate autoregressive state-space models used in the package \code{\link{MARSS-package}}.
}

\details{

The object has the following attributes:
\itemize{
  \item{form}{ The form that the model object is in. }
  \item{par.names}{ The names of each parameter matrix in the model. }
  \item{model.dims}{ A list with the dimensions of all the matrices in non-vectorized form. }
  \item{X.names}{ Names for the X rows. } 
  \item{Y.names}{ Names for the Y rows. } 
  \item{equation}{ The model equation.  Used to write the model in LaTeX. } 
}

These attributes are set in the MARSS_form.R file, in the \code{MARSS.form()} function and must be internally consistent with the elements of the model.  These attributes are used for internal error checking.

Each parameter matrix in a MARSS equation can be written in vectorized form: vec(P) = f + Dp, where f is the fixed part, p are the estimated parameters, and D is the matrix that transforms the p into a vector to be added to f.

An object of class \code{marssMODEL} is a list with elements:
\itemize{
  \item{data}{ Data supplied by user. }
  \item{fixed}{ A list with the f row vectors for each parameter matrix.  }
  \item{free}{ A list with the D matrices for each parameter matrix. }
  \item{tinitx}{ At what t, 0 or 1, is the initial x defined at? } 
  \item{diffuse}{ Whether a diffuse initial prior is used.  TRUE or FALSE.  Not used unless \code{method="BFGS"} was used. } 
  }

For the marss form, the matrices are called:   Z, A, R, B, U, Q, x0, V0.  This is the form used by all internal algorithms, thus alternate forms must be transformed to marss form before fitting.  For the marxss form (the default form in a \code{\link{MARSS}()} call), the matrices are called: Z,  A,  R,  B,  U, Q, D, C, d, c, x0, V0.

Each form, should have a file called MARSS_form.R, with the following functions.  Let foo be some form.
\itemize{
  \item{MARSS.foo(MARSS.call)}{ This is called in \code{\link{MARSS}()} and takes the input from the \code{\link{MARSS}()} call (a list called MARSS.call) and returns that list with two model objects added. First is a model object in marss form in the $marss element and a model object in the form foo.}
  \item{marss_to_foo(marssMLE or marssMODEL)}{ If called with marssMODEL (in form marss), marss_to_foo returns a model in form foo.  If marss_to_foo is called with a \code{\link{marssMLE}} object (which must have a $marss element by definition), it returns a $model element in form foo and all if the marssMLE object has par, par.se, par.CI, par.bias, start elements, these are also converted to foo form. The function is mainly used by print.foo which needs the par (and related) elements of a marssMLE object to be in foo form for printing.}
  \item{foo_to_marss(marssMODEL or marssMLE)}{ This converts marssMODEL(form=foo) to marssMODEL(form=marss).  This transformation is always possible since MARSS only works for models for which this is possible. If called with marssMODEL, it returns only a \code{\link{marssMODEL}} object.  If called with a \code{\link{marssMLE}} object, it adds the \code{$marss} element with a \code{\link{marssMODEL}} in "marss" form and if the par (or related) elements exists, these are converted to "marss" form. }
  \item{print_foo(marssMLE or marssMODEL)}{ print.marssMLE prints the par (and par.se and start) element of a marssMLE object but does not make assumptions about its form. Normally this element is in form=marss.  print.marssMLE checks for a print_foo function and runs that on the marssMLE object first.  This allows one to call foo_to_marss() to covert the par (and related) element to foo form so they look familiar to the user (the marss form will look strange). If called with marssMLE, print_foo returns a marssMLE object with the par (and related) elements in foo form.  If called with a marssMODEL, print_foo returns a marssMODEL in foo form.} 
  \item{coef_foo(marssMLE)}{ See print_foo.  coef.marssMLE also uses the par (and related) elements. }
  \item{predict_foo(marssMLE)}{ Called by predict.marssMLE to do any needed conversions.  Typically a form will want the newdata element in a particular format and this will need to be converted to marss form. This returns an updated marssMLE object and newdata.}
 \item{describe_foo(marssMODEL)}{ Called by describe.marssMODEL to do allow custom description output.} 
 \item{is.marssMODEL_foo(marssMODEL)}{ Check that the model object in foo form has all the parts it needs and that these have the proper size and form.}
  \item{MARSSinits_foo(marssMLE, inits.list)}{ Allows customization of the inits used by the form. Returns an inits list in marss form.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "marssMODEL")}: ... }
    \item{summary}{\code{signature(object = "marssMODEL")}: ... }
    \item{toLatex}{\code{signature(object = "marssMODEL")}: ... }
    \item{model.frame}{\code{signature(object = "marssMODEL")}: ... }
	 }
}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\keyword{classes}
