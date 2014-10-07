#==============================================================================
#                   CSIZE
#
#
#' Size on list of matrices
#'
#' Returns dimension information for list of arrays.
#'
#' @param b List of matrices
#' @return List with fields:
#'    \itemize{
#'      \item \code{dimension} scalar, number of matrices
#'      \item \code{nrows} array, number  of rows per dimension
#'      \item \code{ncols} array, number  of columns per dimension
#'    }
#'
#' @author Randall Romero-Aguilar, based on Miranda & Fackler's CompEcon toolbox
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance


csize <- function(b){
  if (is.list(b)){
    d <- length(b)
    m <- t(sapply(b,dim))
    return(list(dimension=d,nrows=m[,1],ncols=m[,2]))
    } else{
      d <- NULL
      m <- dim(b)
      return(list(dimension=d,nrows=m[1],ncols=m[2]))
    }

}
