#==========CLASS DEFINITION===============================================
setClass("basisCheb",
         contains="basis",
         representation(
           nodetype = "character"  # type of nodes
         )
)

#==========CLASS METHODS==================================================

#_________________________________________________________________________
#////////// DEFINE A NEW BASIS ///////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.define","basisCheb",
  function(OBJECT,n,a,b,nodetype){
    OBJECT@type <- 'cheb'
    OBJECT@n <- if (hasArg(n)) n else stop('Number of nodes missing!  n = ?')
    OBJECT@a <- if (hasArg(a)) a else stop('Lower bound missing!  a = ?')
    OBJECT@b <- if (hasArg(b)) b else stop('Upper bound missing!  b = ?')

    if (!hasArg(nodetype)){
      OBJECT@nodetype = 'canonical'
    } else {
      if(tolower(nodetype) %in% c('canonical','endpoint','lobatto')){
        OBJECT@nodetype <- tolower(nodetype)
      } else {
        warning(paste('Unknown nodetype',nodetype,": using 'canonical instead'"))
      }
    }
    return(OBJECT)
  }
)


#_________________________________________________________________________
#////////// COMPUTE NODES ////////////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  'basis.setnodes','basisCheb',
  function(
    OBJECT # basis or reference to a basis
  ){
    s <- (OBJECT@b - OBJECT@a)/2
    m <- (OBJECT@b + OBJECT@a)/2
    n <- OBJECT@n
    a <- OBJECT@a
    b <- OBJECT@b

    if (OBJECT@nodetype %in% c('canonical','endpoint')){
      k <- pi*(1:n - 0.5)
      x <- m - cos(k/n)*s
      if (OBJECT@nodetype=='endpoint')  {      # Extend nodes to endpoints
        aa <- x[1]
        bb <- x[n]
        x <- (bb*a - aa*b)/(bb-aa) + (b-a)/(bb-aa)*x
      }
    } else {  # Lobatto nodes
      k <- pi*(0:(n-1))
      x <- m - cos(k/(n-1))*s
    }
    OBJECT@nodes <- matrix(x,ncol=1)
    return(OBJECT)
  }
)



#_________________________________________________________________________
#////////// COMPUTE OPERATORS TO DIFFERENTIATE AND INTEGRATE//////////////
#_________________________________________________________________________
setMethod(
  "basis.diffoptr","basisCheb",
  function(
    OBJECT,        # basis or reference to a basis
    order = 1,     # order of derivative
    integrate = FALSE # integral instead of derivative
  ){

    derivative = !integrate

    # Use previously stored values for orders
    # this speeds up calculations considerably with repeated calls

    if(any(order<0)) stop('In basis.diffoptr: order must be nonnegative')

    if(order==0)   return(OBJECT)
    if(derivative && length(OBJECT@D) >= order) return(OBJECT)
    if(integrate  && length(OBJECT@I) >= order) return(OBJECT)

    n <- OBJECT@n
    a <- OBJECT@a
    b <- OBJECT@b

    if(n-order < 2 && derivative){
      warnmsg <- paste("Insufficient nodes: truncating order to k = ",n-2)
      warning(warnmsg)
      order <- n-2
    }

    if (derivative){
      #===============TAKE DERIVATIVE=========
      i <- matrix(rep(1:n,n),n)
      j <- t(i)
      rc <- which((i+j)%%2 & j>i,arr.ind=TRUE)

      d <- matrix(0,n-1,n)
      d[rc] <- (4/(b-a))*(j[rc]-1)
      d[1,] <- d[1,]/2

      d <- Matrix::Matrix(d,sparse=TRUE)

      if(length(OBJECT@D)==0) OBJECT@D[[1]] <- d

      ii <- length(OBJECT@D) + 1 # index of missing OBJECT@D
      while (ii<=order){ # OBJECT D still incomplete
        OBJECT@D[[ii]] <- Matrix::Matrix(d[1:(n-ii),1:(n-ii+1)] %*% OBJECT@D[[ii-1]],sparse=TRUE)
        ii <- ii+1
      }

    } else{
      #===============TAKE INTEGRAL=========
      nn <- n + order
      i <- (0.25*(b-a))/seq(1,nn)
      diags <- list(i,-i[1:(nn-2)])
      d <- bandSparse(n=nn,m=nn,k=c(0,2),diagonals=diags)
      d[1,1] <- 2*d[1,1]
      d0 <- ((-1)^(0:(nn-1)))*apply(d,MARGIN=2,sum)

      if(length(OBJECT@I)==0) OBJECT@I[[1]] <- rBind(d0[1:n],d[1:n,1:n])

      ii <- length(OBJECT@I)+1 # index of missing OBJECT@I

      while (ii <=order){
        temp <- rBind(d0[1:(n+ii-1)],d[1:(n+ii-1),1:(n+ii-1)])
        OBJECT@I[[ii]] <- Matrix::Matrix(temp %*% OBJECT@I[[ii-1]],sparse=TRUE)
        ii <- ii+1
      }
    }

    return(OBJECT)
  }
)


#_________________________________________________________________________
#////////// COMPUTE INTERPOLATION MATRICES ///////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.getPhi","basisCheb",
  function(
    OBJECT,              # basis OBJECT
    x = OBJECT@nodes,    # vector of the evaluations points
    order = 0,            # order of derivatives
    integrate = FALSE # integral instead of derivative
  ){
    derivative = !integrate
    maxorder <- max(order)
    nn  <- OBJECT@n + (if (integrate) maxorder else 0)

    # Compute operators for derivative/integral
    basis.diffoptr(ref::as.ref(OBJECT),maxorder,integrate)

    # Compute order 0 basis
    if (hasArg(x) || OBJECT@nodetype  !='canonical'){ # evaluate at arbitrary nodes
      x <- as.matrix(x)
      z <- (2/(OBJECT@b - OBJECT@a))* (x-(OBJECT@a + OBJECT@b)/2)
      m <- nrow(z)
      bas <- matrix(0,m,nn)
      bas[,1] <- matrix(1,m,1)
      bas[,2] <- z
      z <- 2*z
      for (i in 3:nn) bas[,i] <- z*bas[,i-1] - bas[,i-2]
    } else { # evaluate at standard Gaussian nodes
      temp <- matrix(seq(OBJECT@n-0.5,0.5,-1),ncol=1)
      bas  <- cos((pi/OBJECT@n) * temp %*% seq(0,nn-1))
    }


    B <- list()

    for (ii in 1:length(order)){

      sg.order <- if (integrate) order else -order  #signed order
      # Compute Phi ////////////////////////
      B[[ii]] <- if(order[[ii]]==0){
        bas[,1:OBJECT@n]
        } else {
          bas[,1:(OBJECT@n + sg.order[ii])] %*%
            (if (integrate) OBJECT@I[[order[ii]]] else OBJECT@D[[order[ii]]])
        }

      #  Label columns and rows in output matrices
      colPrefix <- if(order[ii]==0) {
        'phi'
      } else {
        paste(if(integrate) 'I' else 'D',
              if(order[ii]>1) order[ii] else '',
              'phi',sep='')
      }

      if (length(x)==1) B[[ii]] <- matrix(B[[ii]],nrow=1) #otherwise the dimension is dropped!

      colnames(B[[ii]]) <- paste(colPrefix,1:OBJECT@n,sep='')
      rownames(B[[ii]]) <- paste('x',1:nrow(x),sep='')
    }


    if (length(B)==1) return(B[[1]]) else return(B)
  }
)


#_________________________________________________________________________
#////////// EVALUATE FUNCTION WITH BASIS// ///////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.evaluate","basisCheb",
  function(
    OBJECT,            # basis OBJECT
    x = OBJECT@nodes,  # vector of the evaluations points
    coef,              # coefficient matrix
    order = 0,         # order of derivatives
    complete = FALSE   # compute complete polynomials or full tensor product
  ){


    # Get interpolation matrix
    B <- basis.getPhi(OBJECT,x,order,complete)


    y <- if (is.list(B)) lapply(B,function(b) b%*% coef) else B %*% coef
    return(y)
  }
)

