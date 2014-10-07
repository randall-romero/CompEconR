bisect <- function(f,a,b,...,opt=list(tol=1e-10,showiters=FALSE)){
  # Get options
  if (is.null(opt)) opt <- list()
  if (is.list(opt)){
    tol       <- if(is.null(opt$tol)) 1e-10        else opt$tol
    showiters <- if(is.null(opt$showiters)) FALSE  else as.logical(opt$showiters)
  } else {
    stop('opt must be a list object, with optional elements "tol"')
  }
  
  # Perform checks
  nargin <- length(as.list(match.call())) -1
  if (nargin<3) stop('In bisect: At least three parameters must be passed.')
  
  # lower point < upper point?
  if (a > b) stop('In bisect: Upper bound must be greater than lower bound.')
  
  # end points with different signs?
  sa  <-  sign(f(a,...))
  sb  <-  sign(f(b,...))
  if (sa == sb) stop('In bisect: Function has same sign at interval endpoints.')


  # Initializations
  dx   <-  (b-a)/2
  tol  <-  dx * tol
  x    <-  a + dx
  dx   <-  sb * dx
  
  
  if (showiters) {it=0; cat(sprintf('%4s %8s %8s\n','it','x','dx'))}
  
  # Iteration loop
  while (abs(dx) > tol) {
    if (showiters) {it <- it+1; cat(sprintf('%4i %8.2e %8.2e\n',it,x,dx))}
    dx <- 0.5*dx;
    x  <-  x - sign(f(x,...))*dx;
  }

  return(x)
}