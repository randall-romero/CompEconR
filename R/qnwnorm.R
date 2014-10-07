#==============================================================================
#                   QNWNORM
#
#' Gaussian quadrature nodes and weights for normal distribution
#'
#' Generates Gaussian quadrature nodes and probability weights for multivariate
#' normal distribution with mean vector \code{mu} and variance matrix \code{var}
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
#' # To compute Ef(X) when f is real-valued and X is Normal(mu,var) on R^d,
#' # write a function f that returns m.1 vector when passed an m.d matrix, and write
#' q  <- qnwnorm(n,mu,var)
#' Ef  <- crossprod(q$w, f(q$x))
#'
#' # Alternatively, use the quadrature function
#' Ef <- quadrature(f,qnwnorm,n,mu,var)


qnwnorm <- function(n,mu=matrix(0,1,length(n)),var=diag(length(n))){

  # Auxiliary function
  qnwnorm1 <- function(n){
    maxit <- 100
    pim4  <- 1/pi^0.25
    m <- (n+1) %/% 2
    x <- matrix(0,n,1)
    w <- matrix(0,n,1)


    for (i in 1:m) {
      # Reasonable starting values
      icase <- min(i,5)
      z <- switch(icase,
                  sqrt(2*n+1) - 1.85575*((2*n+1)^(-1/6)), # i=1
                  z - 1.14*(n^0.426)/z,     # i=2
                  1.86*z + 0.86*x[1],       # i=3
                  1.91*z+0.91*x[2],         # i=4
                  2*z+x[i-2])               # i > 4

      # Rootfinding iterations
      it <- 0
      while (it<maxit){
        it <- it+1
        p1 <- pim4
        p2 <- 0
        for (j in 1:n){
          p3 <- p2
          p2 <- p1
          p1 <- z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3
        }
        pp <- sqrt(2*n)*p2
        z1 <- z
        z <- z1 - p1/pp
        if (abs(z-z1) < 1e-14) break
      }
      if (it >= maxit) stop('In qnwnorm: failure to converge.')
      x[n+1-i] <- z
      x[i] <- -z
      w[i] <- 2/(pp*pp)
      w[n+1-i] <- w[i]

  }
    w <- w/sqrt(pi)
    x <- x*sqrt(2)
    return(list(xpoints=x,weights=as.matrix(w)))
  }

  # Main function
  d <- length(n)
  X <- list()
  W <- list()

  for (k in 1:d){
    temp <- qnwnorm1(n[k])
    X[[k]] <- temp$xpoints
    W[[k]] <- temp$weights
  }

  w  <- ckron(rev(W))
  x <- expand.grid(X)

  x = as.matrix(x) %*% chol(var) + matrix(mu,nrow=prod(n),ncol=d)

  return(list(xpoints=x,weights=w))
}
