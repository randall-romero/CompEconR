#==============================================================================
#                   QNWLEGE
#
#' Gauss-Legendre quadrature nodes and weights
#'
#' Generates Guass-Legendre quadrature nodes and weights for computing the
#' definite integral of a real-valued function defined on a hypercube [a,b] in R^d.
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
#' #' # To compute definte integral of a real-valued function f defined on a hypercube
#' # [a,b] in R^d, write a function f that returns an m.1 vector when passed an
#' # m.d matrix, and write
#' q <- qnwlege(n,a,b);
#' Intf <- crossprod(q$w, f(q$x))
#'
#' # Alternatively, use the quadrature function
#' Intf <- quadrature(f,qnwlege,n,a,b)


qnwlege <- function(n,a=matrix(0,1,length(n)),b=matrix(1,1,length(n))){

  # Auxiliary function
  qnwlege1 <- function(n,a,b){
    if (n<=1) stop('In qnwlege: n must be integer greater than one.')

    maxit <- 100
    m <- (n+1) %/% 2
    xm <- 0.5*(b+a)
    xl <- 0.5*(b-a)
    x <- matrix(0,nrow=n,ncol=1)
    w <- x
    i <- 1:m
    z <- cos(pi*(i-0.25)/(n+0.5))
    for (its in 1:maxit){
      p1 <- 1
      p2 <- 0
      for (j in 1:n){
        p3 <- p2
        p2 <- p1
        p1 <- ((2*j-1)*z*p2-(j-1)*p3)/j
      }
      pp <- n*(z*p1-p2)/(z*z-1)
      z1 <- z
      z <- z1-p1/pp
      if (all(abs(z-z1)<1e-14)) break
    }
    if (its==maxit) stop('In qnwlege: failure to converge.')
    x[i] <- xm-xl*z
    x[n+1-i] <- xm+xl*z
    w[i] <- 2*xl/((1-z*z)*pp*pp)
    w[n+1-i] <- w[i]
    return(list(xpoints=x,weights=as.matrix(w)))
  }

  # MAIN FUNCTION
  if (any(a>b)) stop('In qnwlege: right endpoints must exceed left endpoints.')

  d <- length(n)
  X <- list()
  W <- list()

  for (k in 1:d){
    temp <- qnwlege1(n[k],a[k],b[k])
    X[[k]] <- temp$xpoints
    W[[k]] <- temp$weights
  }

  w  <- ckron(rev(W))
  x <- expand.grid(X)
  return(list(xpoints=x,weights=w))
}
