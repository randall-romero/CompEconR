#==============================================================================
#                   QNWUNIF
#
#' Gaussian quadrature nodes and weights for uniform distribution
#'
#' Generates Gaussian quadrature nodes and probability weights for multivariate
#' Uniform distribution on hypercube [a,b] in R^d.
#'
#' @param n 1.d number of nodes per dimension
#' @param a 1.d left endpoints
#' @param b   1.d right endpoints
#' @return List with fields
#'  \itemize{
#'      \item \code{xpoints} prod(n).d quadrature nodes
#'      \item \code{weights} prod(n).1 quadrature weights
#'      }
#' @family quadrature functions
#' @keywords quadrature
#'
#' @author Randall Romero-Aguilar, based on Miranda & Fackler's CompEcon toolbox
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance
#'
#' @examples
#' # To compute Ef(X) when f is real-valued and X is Uniform(a,b) on R^d, write a
#' # function f that returns m.1 vector when passed an m.d matrix, and write
#' q <- qnwunif(n,a,b)
#' Ef <- crossprod(q$w,f(q$x))
#'
#' # Alternatively, use the quadrature function
#' Ef <- quadrature(f,qnwunif,n,mu,var)


qnwunif <- function(n,a=matrix(0,1,length(n)),b=matrix(1,1,length(n))){
  d <- length(n)
  if (any(a>b)) stop('In qnwunif: right endpoints must exceed left endpoints.')
  temp <- qnwlege(n,a,b)

  return(list(xpoints=temp$x,weights=as.matrix(temp$w/prod(b-a))))
}
