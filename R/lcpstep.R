#' Newton step for Array Linear Complementarity Problem


lcpstep <- function(method,x,xl,xu,F,Fx){
  dim(x) <- dim(x)[1:2]

  if (!hasArg(Fx)){ # derivatives not needed
    F <- switch(method,
                minmax = pmin(pmax(F,xl-x),xu-x),
                ssmooth= arrayss(x,xl,xu,F),
                default=error(paste("In 'lcpstep': unknown method",method,"; use 'ssmooth' or 'minmax'"))
    )
    return(F)
  } else {
    switch(method,
           minmax = {
             F <- pmin(pmax(F,xl-x),xu-x)
             dx <- -arrayinvb(x,xl,xu,F,Fx)
             dx <- pmin(pmax(dx,xl-x),xu-x)
           },
           ssmooth= {
             FFx <- arrayss(x,xl,xu,F,Fx)
             dx <- -arrayinv(F,Fx);
             dx <- pmin(pmax(dx,xl-x),xu-x);

           },
           default = error(paste("In 'lcpstep': unknown method",method,"; use 'ssmooth' or 'minmax'"))
    )
    return(list(F=F,dx=dx))
  }
}



arrayinvb <- function(x,xl,xu,F,Fx){
  d <- dim(Fx)
  m <- d[1]
  p <- d[2]


  # Make sure inputs have three dimensions
  if (length(d)==2) dim(Fx) <- c(d,1)
  if (length(dim(F))==2) dim(F) <- c(dim(F),1)
  if (length(dim(x))==2) dim(x) <- c(dim(x),1)
  if (length(dim(xl))==2) dim(xl) <- c(dim(xl),1)
  if (length(dim(xu))==2) dim(xu) <- c(dim(xu),1)


  y <- matrix(0,m,p)
  AA <- -diag(p)
  for (i in 1:m){
    A <- Fx[i,,,drop=FALSE]
    dim(A) <- c(p,p)
    b <- F[i,,]
    bl <- xl[i,,] - x[i,,]
    bu <- xu[i,,] - x[i,,]
    ind1 <- b <= bl
    ind2 <- b >= bu
    b[ind1] <- bl[ind1]
    A[ind1,] <- AA[ind1,,drop=FALSE]
    b[ind2] <- bu[ind2]
    A[ind2,] <- AA[ind2,,drop=FALSE]
    y[i,] <- t(solve(A,b))
  }
  return(y)
}
