#  FDJAC
# 
#   Computes two-sided finite difference Jacobian
# 
#   Usage
#     J = fdjac(f,x,varargin)
#   Input
#     f        : name of function of form fval = f(x,varargin)
#     x        : evaluation point
#     varargin : parameters passed to function f
#   Output
#     J        : finite difference Jacobian
#   Options
#     These options may be set by user with OPTSET (defaults in parentheses):
#     tol      : factor used to set step size (eps^(1/3))
# 
#   Copyright(c) 1997-2010
#    Mario J. Miranda - miranda.4@osu.edu
#    Paul L. Fackler  - paul_fackler@ncsu.edu

fdjac <- function(f,x,...,opt=list(tol=(.Machine$double.eps)^(1/3))){
  
  # Get options
  if (is.null(opt)) opt <- list()
  if (is.list(opt)){
    tol <- if(is.null(opt$tol)) (.Machine$double.eps)^(1/3) else opt$tol
    } else {
    stop('opt must be a list object, with optional elements "tol"')
  }
  
  # Compute stepsize
  h  <- tol * max(abs(x),1)
  xh1  <-  x+h
  xh0  <-  x-h
  h <- xh1-xh0 
  
  J <- matrix(0,nrow=length(f(x,...)),ncol=length(x))
  
  # Compute finite difference Jacobian
  for (j in 1:length(x)){
    xx <- x
    xx[j] <- xh1[j]; f1 <- f(xx,...)
    xx[j] <- xh0[j]; f0 <- f(xx,...) 
    J[,j] <- (f1-f0)/h[j]
  }
  return(J)
}







## ORIGINAL MATLAB CODE ########################## 
# function J = fdjac(f,x,varargin)
# 
# % Set options to defaults, if not set by user with OPTSET (see above)
# tol = optget(mfilename,'tol',eps.^(1/3));
# 
# % Compute stepsize
# h = tol.*max(abs(x),1);
# xh1 = x+h; 
# xh0 = x-h;
# h   = xh1-xh0;
# 
# % Compute finite difference Jacobian
# for j=1:length(x);
# xx = x;
# xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
# xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
# J(:,j) = (f1-f0)/h(j);
# end