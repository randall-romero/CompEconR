% Generated by roxygen2 (4.0.2): do not edit by hand
\name{qnwsimp}
\alias{qnwsimp}
\title{Simpson's rule quadrature nodes and weights}
\usage{
qnwsimp(n, a = rep(0, length(n)), b = rep(1, length(n)))
}
\arguments{
\item{n}{1.d number of nodes per dimension (must be odd positive integers)}

\item{a}{1.d left endpoints}

\item{b}{1.d right endpoints}
}
\value{
List with fields
 \itemize{
     \item \code{xpoints} prod(n).d quadrature nodes
     \item \code{weights} prod(n).1 quadrature weights
     }
}
\description{
Generates Simpson's rule quadrature nodes and weights for computing the
 definite integral of a real-valued function defined on a hypercube [a,b] in R^d.
}
\examples{
# To compute definte integral of a real-valued function f defined on a hypercube
# [a,b] in R^d, write a function f that returns an m.1 vector when passed an
# m.d matrix, and write
q <- qnwsimp(n,a,b,type);
Intf <- crossprod(q$w, f(q$x))

# Alternatively, use the quadrature function
Intf <- quadrature(f,qnwnsimp,n,a,b)
}
\author{
Randall Romero-Aguilar, based on Miranda & Fackler's CompEcon toolbox
}
\references{
Miranda, Fackler 2002 Applied Computational Economics and Finance
}
\seealso{
Other quadrature functions: \code{\link{qnwequi}};
  \code{\link{qnwlege}}; \code{\link{qnwlogn}};
  \code{\link{qnwnorm}}; \code{\link{qnwtrap}};
  \code{\link{qnwunif}}; \code{\link{quadrature}}
}
\keyword{quadrature}

