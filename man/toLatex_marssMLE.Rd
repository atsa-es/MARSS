\name{toLatex.marssMODEL}
\alias{toLatex.marssMODEL}
\alias{toLatex.marssMLE}

\title{ Create a LaTeX Version of the Model }
\description{
  Creates LaTex and a PDF (if LaTeX compiler available) using the tools in the Hmisc package. The files are saved in the working directory.
}
\usage{
\method{toLatex}{marssMODEL}(object, ..., file = NULL, digits = 2, greek = TRUE, orientation = "landscape", 
math.sty = "amsmath", output = c("pdf", "tex"), replace = TRUE, simplify = TRUE)
\method{toLatex}{marssMLE}(object, ..., file = NULL, digits = 2, greek = TRUE, orientation = "landscape", 
math.sty = "amsmath", output = c("pdf", "tex"), replace = TRUE, simplify = TRUE)
}
\arguments{
  \item{object}{ A \code{\link{marssMODEL}} or \code{\link{marssMLE}} object.  }
  \item{...}{ Other arguments. Not used. }
  \item{file}{ Name of file to save to. Optional. }
  \item{digits}{ Number of digits to display for numerical values (if real). }
  \item{greek}{ Use greek symbols. }
  \item{orientation}{ Orientation to use. landscape or portrait. }
  \item{math.sty}{ LaTeX math styling to use. }
  \item{output}{ pdf or tex. If blank, both are output. }
  \item{replace}{ Replace existing file if present. }
  \item{simplify}{ If TRUE, then if B or Z are identity, they do not appear. Any zero-ed out elements also do not appear. }
 }

\value{
  A LaTeX and or PDF file of the model.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{ 
# Example with linear constraints
dat <- t(harborSeal)
dat <- dat[c(2:4), ]
Z1 <- matrix(list("1*z1+-1*z2",0,"z2","2*z1","z1",0),3,2)
A1 <- matrix(list("a1",0,0),3,1)
MLEobj <- MARSS(dat, model=list(Z=Z1, A=A1, Q=diag(0.01,2)))
 \dontrun{
 toLatex(MLEobj)
 toLatex(MLEobj$model)
 }

}