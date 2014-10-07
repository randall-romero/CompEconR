#==============================================================================
#                   PRIMES
#
#
#' List of prime numbers
#'
#' @param n integer scalar
#' @return Array with all prime numbers up to \code{n} (inclusive)
#'
#' @author Randall Romero-Aguilar
#'
#' @examples
#' primes(100)


primes <- function(n)
{
  n <- as.integer(n)
  if(n > 1e8) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  fsqr <- floor(sqrt(n))
  while (last.prime <= fsqr)
  {
    primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
    sel <- which(primes[(last.prime+1):(fsqr+1)])
    if(any(sel)){
      last.prime <- last.prime + min(sel)
    }else last.prime <- fsqr+1
  }
  which(primes)
}
