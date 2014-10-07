##  FIXPOINT
# 
#   Computes fixed-point of function using function iteration
#  
#   Usage
#     x = fixpoint(g,x,varargin)
#   Input
#     g        : name of function of form gval=g(x)
#     x        : initial guess for fixed-point
#     varargin : parameters passed to function g
#   Output
#     x        : fixed-point of g
#     gval     : function value at x
#   Options
#     maxit    : maximum number of iterations (100)
#     tol      : convergence tolerance (sqrt(eps))
# 
#   Copyright(c) 1997-2010
#    Mario J. Miranda - miranda.4@osu.edu
#    Paul L. Fackler  - paul_fackler@ncsu.edu



fixpoint <- function(g,x,...,opt=list(tol = sqrt(.Machine$double.eps),maxit=100)){
  # Set options to defaults, if not set by user with opt
  if (is.null(opt)) opt <- list()
  if (is.list(opt)){
    tol <- if(is.null(opt$tol)) sqrt(.Machine$double.eps) else opt$tol
    maxit <- if(is.null(opt$maxit)) 100 else opt$maxit
  } else {
    stop('opt must be a list object, with optional elements "tol" and "maxit"')
  }
  
  norm_vec <- function(z){sqrt(sum(z^2))}
  
  for (it in 1:maxit){
    xold <- x
    x <- g(x,...)
    
    if (!is.finite(norm_vec(x-xold))) {it <- maxit; break} # break if NaN or INF
    if (norm_vec(x-xold)<tol) return(x) 
    #cat(norm_vec(x-xold),'\n')
  }
  
  if (it==maxit) warning('Failure to converge in fixpoint')
  return(x)
} 
