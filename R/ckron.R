#==============================================================================
#                   CKRON
#
#
#' Kronecker product
#'
#' Repeated Kronecker products on a list of matrices.
#'
#' @param b List of matrices
#' @param invert Logical, take inverse if TRUE
#' @return Matrix with Kronecker product
#'
#' @author Randall Romero-Aguilar, based on Miranda & Fackler's CompEcon toolbox
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance


ckron <- function(b,invert=FALSE){
  nargin <- length(as.list(match.call())) -1

  if (nargin<1) stop('At least one parameter must be passed')

  b <- lapply(b,as.matrix)
  bdims  <- csize(b)

  if (invert && any(bdims$nrows != bdims$ncols)) stop('Matrix elements must be square to invert')


  if (is.null(bdims$dim) || bdims$dim==1){
    z <- if(invert) solve(b[[1]]) else b[[1]]
    return(z)
  } else {
    nRows <- cumprod(bdims$nrows)
    nCols <- cumprod(bdims$ncols)
  }

  if (invert) b <- lapply(b,inverse) # invert all input matrices

  z <- matrix(0,nrow=nRows[bdims$d],ncol=nCols[bdims$d])  #preallocate memory
  z[1:nRows[1],1:nCols[1]] <-  b[[1]]  #initialize

  for (k in 2:bdims$dim){
    z[1:nRows[k],1:nCols[k]] <- kronecker(z[1:nRows[k-1],1:nCols[k-1]],b[[k]])
  }
  return(z)
}
