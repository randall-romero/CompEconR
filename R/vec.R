#==============================================================================
#                   VEC
#
#
#' VEC operator: stacks matrix columns into a column vector
#'
#' @param x A m.n matrix
#' @return A mn.1 matrix (column vector)
#'
#' @author Randall Romero-Aguilar
#'
#' @examples
#' y <- matrix(c(1:6),2,3)
#' vec(y)


vec <- function(
  x # matrix
){
  return(matrix(x,length(x),1))
}
