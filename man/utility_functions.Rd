\name{utility.functions}
\alias{utility.functions}
\alias{makediag}
\alias{takediag}
\alias{is.blockdiag}
\alias{is.identity}
\alias{is.diagonal}
\alias{is.equaltri}
\alias{is.validvarcov}
\alias{is.design}
\alias{is.fixed}
\alias{is.zero}
\alias{ is.solvable }
\alias{vec}
\alias{unvec}
\alias{is.wholenumber}
\alias{Imat}
\alias{rwishart}
\alias{ mystrsplit }
\alias{ convert.model.mat }
\alias{ parmat }
\alias{fixed.free.to.formula }
\alias{ matrix.power }
\alias{ sub3D }
\alias{ pinv }
\alias{ pcholinv }
\alias{ pchol }
\alias{ all.equal.vector }
\keyword{internal}

\title{ Matrix Utilities }
\description{
  Matrix utilities for MARSS functions in the \code{\link{MARSS-package}}.  These are not exported but can be accessed using the \code{MARSS:::} prefix.
}
\usage{
is.blockdiag(x)
is.validvarcov(x, method="kem")
is.identity(x, dim=NULL)
is.diagonal(x, na.rm=FALSE)
is.equaltri(x)
makediag(x, nrow=NA)
takediag(x)
is.design(x, strict=TRUE, dim=NULL, zero.rows.ok=FALSE, zero.cols.ok=FALSE)
is.fixed(x, by.row=FALSE)
is.identity(x, dim=NULL)
is.zero(x)
vec(x)
unvec(x, dim=NULL)
is.wholenumber(x, tol = .Machine$double.eps^0.5)
Imat(x)
rwishart(nu, V)
mystrsplit(x)
convert.model.mat(param.matrix)
fixed.free.to.formula(fixed,free,dim)
matrix.power(x, n)
sub3D(x,t=1)
pinv(x)
pcholinv(x)
pchol(x) 
is.solvable(A,y=NULL)
all.equal.vector(x)
parmat(MLEobj, elem = c("B", "U", "Q", "Z", "A", "R", "x0", "V0", "G", "H", "L"), 
       t = 1, dims = NULL, model.loc = "marss")
}
\arguments{
  \item{x, A, y}{ A matrix (or vector for '\code{makediag}' or string for '\code{mystrsplit}'.). }
  \item{na.rm}{ How to treat NAs in the block diag test. }
  \item{dim, dims}{ Matrix dimensions. Some functions will take the vec of a matrix.  In this case, the optional dim arg specifies the matrix dimensions. }
  \item{fixed}{ A fixed matrix per the MARSS specification for fixed matrix syntax. }
  \item{free}{ A free matrix per the MARSS specification for free matrix syntax. }
  \item{nrow}{ Number of rows. }
  \item{tol}{ Tolerance. }
  \item{method}{ kem or BFGS.  Used to add extra test for MARSSoptim(). }
  \item{t}{ The time index or third dimension of a 3D matrix }
  \item{nu, V}{ Parameters of a Wishart distribution. }
  \item{param.matrix}{ The list matrix version of a time-invariant MARSS model. }
  \item{n}{ An interger for the power function. }
  \item{ zero.rows.ok, zero.cols.ok }{ Means the design matrix can have all zero rows or columns. }
  \item{ strict }{ Specifies whether the design matrix must be only 0s and 1s. }
  \item{ by.row }{ For is.fixed, reports whether is.fixed by row rather than for the whole matrix. }
  \item{ MLEobj}{ A \code{\link{marssMLE}} object. }
  \item{ elem }{ The parameter matrix of a marss model to return. }
  \item{ model.loc }{ Whether to use the marss or model marssMODEL in the \code{\link{marssMLE}} object. }
}
\details{
\itemize{
  \item \code{is...} tests for various matrix properties.  isDiagonal() from the Matrix package is used to test numeric matrices for diagonality.  is.diagonal() is only used to determine if list matrices (that combine numeric and character values) are diagonal. is.zero tests for near zeroness and give TRUE for is.zero((.5-.3)-(.3-.1)) unlike ==0.
  \item \code{vec(x)} creates a column vector from a matrix per the standard vec math function.
  \item \code{unvec(c,dim)} takes the vector c and creates a matrix with the specified dimensions.
  \item \code{Imat(nrow)} returns the identity matrix of dimension nrow.
  \item \code{ fixed.free.to.formula } takes a fixed and free pair and constructs a list matrix (or array if time-varying) with formulas in each matrix element.
  \item \code{ convert.model.mat } takes a list matrix with formulas in each element and converts to a fixed/free pair.
  \item \code{ sub3D } returns a 2D matrix after subsetting a 3D matrix on the third (time) dimension. Ensures that R always returns a matrix.
  \item \code{ mystrsplit } is a customized string splitter used by \code{convert.model.mat}.
  \item \code{ rwishart } generates random draws from a wishart distribution.
  \item \code{ matrix.power } is a faster way to get the n-th power of a matrix.
  \item \code{ pinv } is the pseudoinverse based on singular value decomposition PInv=UD^+V' where a diagonal matrix with non-zero diagonal values of D (from svd) replaced with 1/D.
  \item \code{ pcholinv } is the inverse based on the Cholsky decomposition but modified to allow 0s on the diagonal of x (with corresponding 0 row/column).  These appear as 0 row/columns in the returned inverse.
  \item \code{ pchol } returns the Cholsky decomposition but modified to allow 0s on the diagonal of x (with corresponding 0 row/column).
  \item \code{ is.solvable } returns information on the solvability of the linear system y=Ax using the SVD decomposition.
  \item \code{ all.equal.vector } tests if the all the elements in a vector, matrix, or array are all equal.  Works on list matrices too.
  \item \code{ parmat } constructs the parameter matrix with both the fixed and free values from the vectorized form in a \code{\link{marssMLE}} object.  Users should use \code{coef(MLEobj)} (See \code{\link{coef.marssMLE}}).
} 
}
\value{
See above.
}

\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}

