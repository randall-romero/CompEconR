#==============================================================================
#                   BROYDEN
#
#
#'   Computes root of function via Broyden's Inverse Update Method
#'
#' @param f    function, R^d --> R^d
#' @param x    d-array, initial guess for root
#' @param ...  additional arguments to pass to \code{f}
#' @param opt  list, options to control rootfinding process
#'   \itemize{
#'      \item \code{tol} convergence tolerance
#'      \item \code{maxit} maximum number of iterations
#'      \item \code{maxsteps} maximum number of backsteps
#'      \item \code{initb} an initial inverse Jacobian approximation matrix
#'      \item \code{initi} if initb is empty, use the identity matrix to initialize (
#'          if FALSE, a numerical Jacobian will be used)
#'      \item \code{showiters} display results of each iteration if TRUE
#'      }
#' @return A list with fields
#'    \itemize{
#'        \item \code{x}    d-array, root of f
#'        \item \code{fval} d-array, value of \code{f} at solution
#'        \item \code{fjacinv} d.d matrix, inverse of Jacobian estimate
#'    }
#'
#' @author Randall Romero-Aguilar, based on Miranda & Fackler's CompEcon toolbox
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance
#'
#' @examples
#' broyden(function(x) (x-2)^3, 4)
#' broyden(function(x) (x-2)^3, 3, opt = list(showiters=TRUE))
#'
#' # Compute fixedpoint of f(x1,x2)= [x1^2 + x2^3;   x1*x2 - 0.5]
#' # Initial values [x1,x2] = [-1,-1].  True fixedpoint is x1 = -(1/8)^(1/5),  x2 = 1/(2*x1).
#' broyden(function(x) c(x[1]^2 + x[2]^3, x[1]*x[2] - 0.5), c(-1,-1))



broyden <- function(f,x,...,
                    opt=list(
                      tol=sqrt(.Machine$double.eps),
                      maxit=100,
                      maxsteps=25,
                      initb=NULL,
                      initi=FALSE,
                      showiters=FALSE)){

  # Set options to defaults, if not set by user with OPTSET (see above)
  if (is.null(opt)) opt <- list()
  if (is.list(opt)){
    tol       <- if(is.null(opt$tol)) sqrt(.Machine$double.eps) else opt$tol
    maxit     <- if(is.null(opt$maxit)) 100                     else opt$maxit
    maxsteps  <- if(is.null(opt$maxsteps)) 25                   else opt$maxsteps
    initb     <- if(is.null(opt$initb)) NULL                    else opt$initb
    initi     <- if(is.null(opt$initi)) FALSE                   else as.logical(opt$initi)
    showiters <- if(is.null(opt$showiters)) FALSE                else as.logical(opt$showiters)
  } else {
    stop('opt must be a list object, with optional elements "tol", "maxit","maxsteps", "initb","initi","showiters"')
  }



  if (maxsteps<1) maxsteps <- 1  # fix maxsteps

  if (length(initb)==0){       # Set initial inverse Jacobian aprroximation matrix
    fjacinv <- if (initi) -diag(length(x))  else solve(fdjac(f,x,...))
  } else {
    fjacinv <- initb
  }

  norm_vec <- function(z){max(abs(z))}

  fval  <- f(x,...)
  fnorm <- norm_vec(fval)

  if (showiters) {it=0; cat(sprintf('%4s %4s %11s\n','it','back','fnorm'))}

  for (it in 1:maxit){
    if (fnorm<tol) return(list(x=x,fval=fval,fjacinv=fjacinv))
    dx  <- -(fjacinv %*% fval)
    fnormold <- Inf
    for (backstep  in 1:maxsteps){
      fvalnew <- f(x+dx,...)
      fnormnew <- norm_vec(fvalnew)
      if (fnormnew<fnorm) break
      if (fnormold<fnormnew){
        fvalnew <- fvalold
        fnormnew <- fnormold
        dx <- 2*dx
        break
      }
      fvalold <- fvalnew
      fnormold <- fnormnew
      dx <- dx/2
    }
    x <- x+dx
    if (any(!is.numeric(x))) stop('Infinites or NaNs encountered')
    if (fnormnew>fnorm){
      fjacinv  <- if (initi) diag(length(x)) else solve(fdjac(f,x,...))
    } else {
      temp <- fjacinv %*% (fvalnew-fval)
      fjacinv <- fjacinv + (dx-temp) %*% (crossprod(dx,fjacinv) / as.double(crossprod(dx,temp)))
    }
    fval <- fvalnew
    fnorm <- fnormnew
    if (showiters) cat(sprintf('%4i %4i %11.2e\n',it,backstep,fnorm))
  }
  warning('Failure to converge in broyden')
}
