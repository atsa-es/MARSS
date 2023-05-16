\name{ldiag}
\alias{ldiag}
\keyword{user-helper}

\title{ Return a diagonal list matrix }
\description{
  Creates a list diagonal matrix where the diagonal can be a combination of numbers and characters. Characters are names of parameters to be estimated.
}
\usage{
ldiag(x, nrow = NA)
}
\arguments{
  \item{x}{ A vector or list of single values }
  \item{nrow}{ Rows in square matrix }
}
\details{
  A diagonal list matrix is returned. The off-diagonals will be 0 and the diagonal will be \code{x
  }. \code{x} can be a combination of numbers and characters.  If \code{x} is numeric, the diagonal will still be list type so that later the diagonal can be replace with characters. See examples. 
}
\value{
  a square list matrix
}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{

ldiag(list(0, "b"))
ldiag("a", nrow=3)

# This works
a <- ldiag(1:3)
diag(a) <- "a"
diag(a) <- list("a", 0, 0)

# ldiag() is a convenience function to replace having to 
# write code like this
a <- matrix(list(0), 3, 3)
diag(a) <- list("a", 0, 0)

# diag() alone won't work because it cannot handle mixed number/char lists

# This turns the off-diagonals to character "0" 
a <- diag(1:3)
diag(a) <- "a"

# This would turn our matrix into a list only (not a matrix anymore)
a <- diag(1:3)
diag(a) <- list("a", 0, 0)

# This would return NA on the diagonal
a <- diag("a", 3)

}

