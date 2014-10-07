# define auxiliary function
params <- list(r=0.04)

myfunc <- function(flag,s,x,i,j,e,params){
  R <- 1+params$r

  ns <- nrow(s)
  ds <- ncol(s)
  dx <- if(hasArg(x)) ncol(x)

  out <- switch(
    flag,
    f = list(
      f   =  log(s-x/R),
      fx  =  1/(x-R*s),
      fxx = array(-1/(x-R*s)^2,c(ns,dx,dx))
    ),
    g = list(
      g   = x,
      gx  = array(1,c(ns,ds,dx)),
      gxx = array(0,c(ns,ds,dx,dx))
    ),
    b = list(lower=matrix(0,ns),upper=R*s)
  )
  return(out)
}


test_that('create model',{

  mymodel <- new("dynamic.model",
                 horizon=Inf,
                 func=myfunc,
                 discount=0.95,
                 params = params,
                 X=matrix(seq(0,8,by=0.1)))

  expect_equal(mymodel@discount, 0.95)
})

  mymodel <- new("dynamic.model",
                 horizon=Inf,
                 func=myfunc,
                 discount=0.95,
                 params = params
  )



  mybasis <- fundefn(type='cheb',n=7,a=0,b=8)
  #mybasis2 <- fundefn(type='cheb',n=c(4,5),a=c(-1,0),b=c(1,3))

  myoptions <- new("dpsolve.options",algorithm="funcit")
  #myoptions <- new("dpsolve.options",algorithm="newton")


  #myfunc(flag='b',s=mybasis@nodes,i=1,j=1,params=params)
  solution <- dpsolve(model=mymodel,basis=as.ref(mybasis),options=myoptions)




#basis.getPhi(as.ref(mybasis),order=2)
#print(basis.getPhi(as.ref(mybasis),order=-2))

#sprintf('%10s %2d\n','hola',6)
