#==============================================================================
#                   FISCHER
#
#
#' Computes Fischer's function and (optionally) its Jacobian
#'
#' @param u    n-array or vector, first function evaluated at x
#' @param v    n-array or vector, second function evaluated at x
#' @param du   n-array or vector, derivative of first argument evaluated at x (optional)
#' @param dv   n-array or vector, derivative of second argument evaluated at x (optional)
#' @param plus logical, use 'plus' version if TRUE (default), 'minus' if FALSE (optional)
#' @return If inputs \code{du} and \code{dv} are provided, it returns a list with fields
#'  \itemize{
#'      \item \code{f} n.1 matrix, Fischer function  evaluated at x
#'      \item \code{J} n.n matrix, Jacobian evaluated at x
#'      }
#'  Otherwise, it returns only the matrix \code{f}.
#' @family nonlinear equations
#'
#' @author Randall Romero-Aguilar
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance, section 3.8
#' @details Fischer's function \eqn{\phi:R\times R \rightarrow R}  is specified by
#'   \deqn{\phi^{\pm}(u,v) = u + v \pm \sqrt{u^2 + v^2}, u,v\in R}
#'   In this implementation, \code{u,v} are assumed to be function of \code{x}.
#'
#' @examples
#' u   <-  c(1,2,3,4)
#' v   <-  c(1,2,2,1)
#' du  <-  c(0,1,1,0)
#' dv  <-  c(1,1,1,1)
#' fp  <- fischer(u,v,du,dv)
#' fm  <- fischer(u,v,plus=FALSE)


fischer <- function(u, v, du, dv, plus=TRUE){
  s <- if (plus) 1 else -1  # plus or minus sign
  n <- max(length(u),length(v))

  u <- matrix(u,n,1)
  v <- matrix(v,n,1)

  sq.uv <- sqrt(u*u+v*v)

  f <- u + v + s * sq.uv

  if (hasArg(du) && hasArg(dv)){
    u <- matrix(u,n,n)
    v <- matrix(v,n,n)
    sq.uv <- matrix(sq.uv,n,n)

    J <- du + dv + s * (u*du + v*dv) / sq.uv
    return(list(f=f,J=J))
  } else{
    return(f)
  }
}
