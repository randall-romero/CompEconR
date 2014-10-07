#==============================================================================
#                   QUADRATURE
#
#
#' Compute numerical quadrature
#'
#' @param f A function, integrand.
#' @param method A function, specifying the quadrature method.
#' @param n A scalar, number of nodes.
#' @param ... Other parameters passed to the selected \code{method}.
#'
#'   Available \code{method} options are: \code{qnwequi}, \code{qnwlege},
#'   \code{qnwlogn}, \code{qnwnorm}, \code{qnwsimp}, \code{qnwtrap},
#'   \code{qnwunif}
#'
#' @return The integral of function \code{f}.
#' @family quadrature functions
#' @keywords quadrature
#'
#' @author Randall Romero-Aguilar
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance, ch.5
#'
#' @examples
#' quadrature(function(x) sin(x)^2,qnwtrap,50,0,pi)


quadrature <- function(f,method=qnwtrap,n,...){
  A <- method(n,...)
  return(crossprod(A$w,f(A$x)))
}
