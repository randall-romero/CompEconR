#==============================================================================
#                   QNWLOGN
#
#' Gaussian quadrature nodes and weights for lognormal distribution
#'
#' Generates Gaussian quadrature nodes and probability weights for multivariate
#' Lognormal distribution with parameters \code{mu} and \code{var}, or,
#' equivalently, for a random vector whose natural logarithm is normally distributed
#'  with mean vector \code{mu} and variance matrix \code{var}
#'
#' @param n   1.d number of nodes per dimension
#' @param mu  1.d mean vector (default: zeros)
#' @param var d.d positive definite variance matrix (default: identity matrix)
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
#' # To compute Ef(X) when f is real-valued and X is lognormal(mu,var) on R^d,
#' # write a function f that returns m.1 vector when passed an m.d matrix, and write
#' q  <- qnwlogn(n,mu,var)
#' Ef  <- crossprod(q$w, f(q$x))
#'
#' # Alternatively, use the quadrature function
#' Ef <- quadrature(f,qnwlogn,n,mu,var)


qnwlogn <- function(n,mu=matrix(0,1,length(n)),var=diag(length(n))){
  # Generates Gaussian quadrature nodes and probability weights for
  # multivariate Lognormal distribution with parameters mu and var, or,
  # equivalently, for a random vector whose natural logarithm is normally
  # distributed with mean vector mu and variance matrix var.

  temp <- qnwnorm(n,mu,var)
  return(list(xpoints=exp(temp$x),weights=temp$w))
}

