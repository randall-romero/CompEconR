#==============================================================================
#                   QNWTRAP
#
#' Trapezoid rule quadrature nodes and weights
#'
#' Generates Trapezoid rule quadrature nodes and weights for computing the
#'  definite integral of a real-valued function defined on a hypercube [a,b] in R^d.
#'
#' @param n 1.d number of nodes per dimension
#' @param a 1.d left endpoints
#' @param b 1.d right endpoints
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
#' # To compute definte integral of a real-valued function f defined on a hypercube
#' # [a,b] in R^d, write a function f that returns an m.1 vector when passed an
#' # m.d matrix, and write
#' q <- qnwtrap(n,a,b,type);
#' Intf <- crossprod(q$w, f(q$x))
#'
#' # Alternatively, use the quadrature function
#' Intf <- quadrature(f,qnwtrap,n,a,b)





qnwtrap <- function(n,a=rep(0,length(n)),b=rep(1,length(n))){

  qnwtrap1 <- function(ni,ai,bi){
    if (ni<=1) stop('In qnwtrap: n must be integer greater than one.')

    dxi <- (bi-ai)/(ni-1)
    x <- matrix(seq(from = ai, to = bi, by = dxi))
    w <- dxi * matrix(1,ni,1)
    w[c(1,nrow(w))] <- 0.5*w[c(1,nrow(w))]
    return(list(xpoints=x,weights=w))
  }

  if (any(a>b)) stop('In qnwtrap: right endpoints must exceed left endpoints.')

  d <- length(n)
  X <- list()
  W <- list()

  for (k in 1:d){
    temp <- qnwtrap1(n[k],a[k],b[k])
    X[[k]] <- temp$xpoints
    W[[k]] <- temp$weights
  }

  w  <- ckron(rev(W))
  x <- expand.grid(X)

  return(list(xpoints=x,weights=w))
}
