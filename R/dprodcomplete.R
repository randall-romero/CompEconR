dprodcomplete <- function(X,expandedNodes=TRUE){
  # X  is a list of matrices, whose n columns represent evaluated polynomials of degree 0 to n-1
  
  d <- length(X)      # dimension of base
  degs <- sapply(X,ncol) - 1  # degree of polynomials
  degmax <- max(degs)    # maximum degree of product of polynomials (for output)
  
  ldeg <- list()
  for (i in 1:d) ldeg[[i]] <- 0:degs[i]  
  
  deg.grid <- gridmake(ldeg)
  deg <- apply(deg.grid,MARGIN=1,sum) #degrees of all possible column products
  validDeg <- (deg <= degmax)
  deg <- apply(deg.grid[validDeg,],MARGIN=1,paste,sep='',collapse='')
  deg.all <- apply(deg.grid+1,MARGIN=1,paste,sep='',collapse='')
  
  nx <- sapply(X,nrow)
  if (expandedNodes && any(nx!=min(nx))) stop ('In dprodcomplete: to take direct product, all matrices must have same number of rows')
  
  nobs <- if(expandedNodes) nx[1] else prod(nx)
  ncols <- length(deg)
  
  XX <- matrix(0,nobs,ncols)  # preallocate memory

  index <- if (expandedNodes) matrix(1:nobs,nobs,ncols) else gridmake(ldeg) + 1 # row numbers 
  
  vX <- list()  # list of vectors, one from each basis
  for (ii in 1:nobs){
    
    for (k in 1:d) vX[[k]] <- X[[d+1-k]][index[ii,d+1-k],]    # make list in reverse order
  XX[ii,] <- dprodcomplete1(vX,validDeg)
  }
  colnames(XX) <- paste('Phi',deg,sep='')
  rownames(XX) <- paste('x',if (expandedNodes) 1:nobs else deg.all,sep='')
 return(XX)
}



dprodcomplete1 <- function(
  V,  # a list on numeric vectors
  validDeg # a boolean vector with valid degree combinations
  ){
  vgrid <- gridmake(V)
  vgrid <- vgrid[validDeg,]
  vprod <- apply(vgrid,MARGIN=1,prod)
  return (vprod)
}









# OLD VERSIONS


# dprodcomplete <- function(X){
#   d <- length(X)
#   degs <- sapply(X,ncol) - 1
#   degmax <- max(degs)
#   
#   ldeg <- list()
#   for (i in 1:d) ldeg[[i]] <- 0:degs[i]
#   
#   tempdeg <- ldeg[[d]]
#   tempProd <- X[[d]]
#   
#   for (i in (d-1):1){
#     P <- dprodtrunc(x=tempProd,y=X[[i]],dx=tempdeg,dy=ldeg[[i]],dmax=degmax)
#     tempdeg <- P$deg
#     tempProd <- P$prod
#   }
#   
#   return(list(prod=tempProd,deg=tempdeg))
#   
# }





# 
# dprodtrunc <- function(x,y,dx,dy,dmax){
#   nrows <- nrow(x)
#   if(nrows!=nrow(y)) stop('In dprodtrunc: x1 and x2 must have same number of rows')
#   
#   deg <- apply(gridmake(list(dx,dy)),MARGIN=1,sum)
#   validDeg <- (deg <= dmax)
#   deg <- deg[validDeg]
#   
#   ncols <- length(deg)
#   A <- matrix(0,nrows,ncols)
#   
#   for (i in 1:nrows){
#     xy <- gridmake(list(x[i,],y[i,]))
#     xy <- xy[validDeg,]
#     A[i,] <- xy[,1] * xy[,2]
#   }
#   return(list(prod=A,deg=deg))
# }
