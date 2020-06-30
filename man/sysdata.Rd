\name{sysdata}
\alias{hcl_palette_52}
\alias{hcl_palette_100}
\keyword{internal}

\docType{data}
\title{ Palettes }
\description{
Palettes to use in the plotting functions. This is from \code{colorspace::sequential_hcl(n)} with n=100 or n=52.
}
\usage{
data(hcl_palette_52s)
data(hcl_palette_100)
}
\format{
A vector of colors.
}

\examples{
plot(1:10, col=hcl_palette_100[1])
}

