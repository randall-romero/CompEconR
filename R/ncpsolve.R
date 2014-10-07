ncpsolve <- function(f,                               # name of user-supplied function
                     a,                               # d by 1 vector, left bound on x
                     b,                               # d by 1 vector, right bound on x
                     x = (a+b)/2,                     # d by 1 vector, initial guess for solution
                     ...,                             # optional parameters passed to f
                     maxit=100,                       # maximum number of iterations
                     tol = sqrt(.Machine$double.eps), # convergence tolerance
                     maxsteps = 10,                   # maximum number of backsteps
                     showiters = FALSE,               # display results of each iteration
                     type = 'ssmooth',                 # rootproblem transform, 'ssmooth' or 'minmax'
                     method = 'newton'                # method for rootfinding, 'newton' or 'funcit' 
){
  switch(method,
         newton = {
           for (it in 1:maxit){
             fx <- f(x,...)
             G <- switch(type,
                         ssmooth = ssmooth(x,a,b,fx$f,fx$J),
                         minmax  =  minmax(x,a,b,fx$f,fx$J)
             )
             
             fnorm <- max(abs(G$f))
             if (fnorm < tol) break
             dx <- - solve(G$J,G$f)
             fnormold <- Inf
             for (backstep in 1:maxsteps){
               xnew <- x + dx
               fxnew <- f(xnew,...)
               fnew <- switch(type,
                              ssmooth = ssmooth(xnew,a,b,fxnew$f,fxnew$J),
                              minmax  =  minmax(xnew,a,b,fxnew$f,fxnew$J)
               )
               fnormnew <- max(abs(fnew$f))
               if (fnormnew<fnorm) break
               if (fnormold<fnormnew) {dx <- 2*dx; break}
               fnormold <- fnormnew
               dx <- dx/2
             }
             x <- x + dx
             if (showiters) cat(sprintf('%4i %4i %6.2e\n',it,backstep,fnormnew)) 
           }
           
           if (it==maxit) warning('Failure to converge in ncpsolve')
           
           x <- Re(x)
           x <- pmax(a,x)
           x <- pmin(b,x)
           return(list(x=x,fx=fnew$f))
         },
         funcit={
           ff <- function(x){
             ffx <-  x + 
               switch(type,
                      ssmooth = ssmooth(x,a,b,f(x,...)$f),
                      minmax  =  minmax(x,a,b,f(x,...)$f)
               )
             return(ffx)
           }
           x <- fixpoint(ff,x)
           return(x)
         }
  )
}