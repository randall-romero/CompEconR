#==============================================================================
#                   QNWEQUI
#
#
#' Quadrature nodes and weights using equidistributed nodes.
#'
#' Generates equidistributed nodes for computing the definite integral of
#'   a real-valued function defined on a hypercube [a,b] in R^d.
#'
#' @param n total number of nodes
#' @param a 1.d left endpoints
#' @param b 1.d right endpoints
#' @param type string, type of sequence. Choose from:
#' \itemize{
#'    \item 'N' - Neiderreiter (default)
#'    \item 'W' - Weyl
#'    \item 'H' - Haber
#'    \item 'R' - pseudo Random
#'    }
#'
#'  @return List with fields
#'  \itemize{
#'      \item \code{xpoints} n.d quadrature nodes
#'      \item \code{weights} n.1 quadrature weights
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
#' q <- qnwequi(n,a,b,type);
#' Intf <- crossprod(q$w, f(q$x))
#'
#' # Alternatively, use the quadrature function
#' Intf <- quadrature(f,qnwnequi,n,a,b,type)


qnwequi <- function(n,                          # number of nodes
                    a = matrix(0,1,length(n)),  # lower bound (zeros)
                    b = matrix(1,1,length(n)),  # upper bound (ones)
                    type = 'N',                 # type (Neiderreiter)
                    equidist_pp = sqrt(primes(7920))  # good for d<=1000
){
  if (any(a>b)) stop('In qnwequi: right endpoints must exceed left endpoints.')

  d <- max(length(n),max(length(a),length(b)))
  n <- prod(n)
  i = matrix(1:n,ncol=1)

  decimalPart <- function(x) x - floor(x)

  x <- switch(toupper(substr(type,1,1)),
              N = decimalPart(i %*% 2^((1:d)/(d+1))),           # Neiderreiter
              W = decimalPart(i %*% equidist_pp[1:d]),          # Weyl
              H = decimalPart((i*(i+1)/2) * equidist_pp[1:d]),  # Haber
              R = matrix(runif(n*d),n,d),                       # pseudo-random
              error('In qnwequi: unknown sequence requested.')  # otherwise
  )

  r <- b-a
  x <- matrix(rep(a,n),nrow=n,byrow=TRUE) + x * matrix(rep(r,n),nrow=n,byrow=TRUE)
  w <- rep((prod(r)/n),n)
  return(list(xpoints=x,weights=as.matrix(w)))
}
