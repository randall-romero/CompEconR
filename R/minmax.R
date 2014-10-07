#==============================================================================
#                   MINMAX
#
#' min-max approximation
#'
#' Computes \eqn{\tilde f(x) = min (max ( f(x), a-x ), b-x)} and its derivative
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


minmax <- function(x,a,b,fx,J){
  nx <- length(x)

  if (length(a)==1) a <- rep(a,nx)
  if (length(b)==1) b <- rep(b,nx)

  da <- a - x
  db <- b - x
  fhatval <- pmin(pmax(fx,da),db) #transformed function value at x

  if (hasArg(J)){
    fhatjac <-  -diag(nx)
    i  <- which(fx>da & fx <db)
    fhatjac[i,] = J[i,] # transformed function Jacobian at x
    return(list(f=fhatval,J=fhatjac))
  } else {
    return(fhatval)
  }
}
