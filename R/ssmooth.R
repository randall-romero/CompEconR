#==============================================================================
#                   SSMOOTH
#
#' Semi-smooth approximation
#'
#' Computes \eqn{\tilde f(x) = \phi^-(\phi^+(f(x),a-x),b-x)} and its derivative
#'
#' @param x   n-array of evaluation points
#' @param a   n-array of lower bounds
#' @param b   n-array of upper bounds
#' @param fx  n-array of function \code{f} evaluated at f
#' @param J   n.n matrix of Jacobian of \code{f} (optional)
#' @return If input \code{J} is provided, it returns a list with fields
#'  \itemize{
#'      \item \code{f} n.1 matrix, function  evaluated at x
#'      \item \code{J} n.n matrix, Jacobian  evaluated at x
#'      }
#'  Otherwise, it returns only the matrix \code{f}.

#' @family nonlinear equations
#'
#' @author Randall Romero-Aguilar
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance, section 3.8
#'


ssmooth <- function(x, a, b, fx, J){

  nx <- length(x)
  ainf <- is.infinite(a)
  binf <- is.infinite(b)

  if (hasArg(J)){
    # apply fischer plus
    fplus <- fischer(fx,a-x,J,-diag(nx),plus=TRUE)
    fplus$f[ainf] <- fx[ainf]
    fplus$J[ainf] <- J[ainf]

    # apply fischer minus
    fhat <- fischer(fplus$f,b-x,fplus$J,-diag(nx),plus=FALSE)
    fhat$f[binf] <- fplus$f[binf]
    fhat$J[binf] <- fplus$J[binf]
    return(fhat)
  } else {
    # apply fischer plus
    fplus <- fischer(fx,a-x,plus=TRUE)
    fplus[ainf] <- fx[ainf]

    # apply fischer minus
    fhat <- fischer(fplus,b-x,plus=FALSE)
    fhat[ainf] <- fplus[binf]
    return(fhat)
  }
}
