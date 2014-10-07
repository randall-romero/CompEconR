fundefn <- function(
  type,  # 1.d string indicating basis function type ('cheb','spli' or 'lin')
  n,     # 1.d vector of integers > 1 indicating order of approximation per dimension
  a,     # 1.d vector of left endpoints of approximation intervals per dimension
  b,     # 1.d vector of right endpoints of approximation intervals per dimension
  k = 3, # for 'spli' basistype, the order of the spline (default,   # 3 for cubic)
  format = 'tensor'  # if more than one dimension, how to combine them
){

  ## FUNDEFN
  #   Defines a basis function family.

  type <- tolower(substr(type,1,4))
  if (!all(type %in% c('cheb','spli','lin'))) stop("basis type must be 'cheb','spli', or 'lin'")


  d <- length(n)
  if (length(a) != d) stop('a must be same dimension as n')
  if (length(b) != d) stop('a must be same dimension as n')
  if (any(a>b)) stop('left point must be less than right point')
  if (any(n<2)) stop('n(i) must be greater than one')

  if(length(type==1)) type <- rep(type,d)  # if only one type provided, use it in all dimensions

  bas <- list()
  for (i in 1:d){
    bas[[i]] <-  # set a one-dimension basis
      switch(type[i],
             cheb = basis('cheb', n=n[i],a=a[i],b=b[i]),
             spli = basis('spli', breaks=c(a[i],b[i]),evennum=n[i]-k+1,k=k),
             lin  = basis('lin', breaks=c(a[i],b[i]), evennum=n[i])
      )
  }

  if (d==1){
    return(bas[[1]])
  } else{
    newBasis <- new("basisN")
    newBasis@d = d
    newBasis@n = n
    newBasis@a = a
    newBasis@b = b
    newBasis@type = type
    newBasis@bas = bas
    newBasis@format = format
    basis.setnodes(ref::as.ref(newBasis))
    return(newBasis)
  }
}
