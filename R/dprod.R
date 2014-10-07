dprod <- function(
  a,   # k.p matrix
  b    # k.q matrix
  ){
  if (nrow(a)!=nrow(b)) stop('In drprod: arguments should be two matrices with same number of rows')

  a1 <- matrix(1,nrow=1,ncol=ncol(a))
  b1 <- matrix(1,nrow=1,ncol=ncol(b))
  
  c <- kronecker(a,b1) * kronecker(a1,b)
  return(c)
}
