#===========================================================================
# DEFINE A d-DIMENSIONAL BASIS

setClass("basisN",
         representation(d       = "numeric",    # dimension of basis
                        type    = "character",  # type of basis 
                        n       = "numeric",    # number of nodes
                        a       = "numeric",    # lower limits
                        b       = "numeric",    # upper limits
                        bas     = "list",       # list of one-dimensional basis      
                        format  = "character",   # how to combine 1D basis: 'tensor','complete'
                        nodes   = "matrix"     # nodes         
         )
)


#=========================================================================
# METHODS FOR BASISN, TO BE DEFINE WITH SETMETHOD

basis.product  <- function(OBJECT,expandedNodes,completePoly) 'multiplies d interpolation bases'
  # NOTE: basis.product OBJECT is of class 'list', and it's called by the
  # basis.getPhi method for class 'basisN'
  

#_________________________________________________________________________
#////////// DISPLAY BASIS ////////////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  "show","basisN",
  function(object){
    strf <- paste('\t%5s %12s %15s %15s\n')
    cat("Basis of ", object@d, "dimensions\n")
    cat(sprintf(strf,'Dim.','type','# of nodes','interval'))
    
    for (d in 1:object@d){
      longtype <- switch(object@type[d],
                         cheb = 'Chebychev',
                         spli = paste('Spline(',object@bas[[d]]@k,')',sep=''),
                         lin = 'Linear')
      interval <- paste('[',object@a[d],', ',object@b[d],']',sep='')
      cat(sprintf(strf,d,longtype,object@n[d],interval))
    }
  }
)





#_________________________________________________________________________
#////////// COMPUTE NODES ////////////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  'basis.setnodes','basisN',
  function(
    OBJECT # basis or reference to a basis
  ){
    
    OBJECT@d <- length(OBJECT@bas)
    
    if(OBJECT@d < 1) stop("slot 'bas' must have at least one basis")
    
    theNodes = list()
    for(d in 1:OBJECT@d){
      if(length(OBJECT@bas[[d]]@nodes) != OBJECT@n[d]){
        OBJECT@bas[[d]]  <- basis.setnodes(OBJECT@bas[[d]])
      }
      theNodes[[d]] <- OBJECT@bas[[d]]@nodes
    }  
    OBJECT@nodes <- gridmake(theNodes)
    return(OBJECT)
  }
)

#_________________________________________________________________________
#////////// COMPUTE OPERATORS TO DIFFERENTIATE AND INTEGRATE//////////////
#_________________________________________________________________________
setMethod(
  "basis.diffoptr","basisN",
  function(
    OBJECT,    # basis or reference to a basis
    order = rep(1,OBJECT@d)  # order of derivative
  ){
    OBJECT@d <- length(OBJECT@bas)
    
    if(OBJECT@d != length(order)) stop("order should have same dimension as basis")
    
    for(d in 1:OBJECT@d){
      OBJECT@bas[[d]]  <- basis.diffoptr(OBJECT@bas[[d]],order=order[d])
    }
    return(OBJECT)
  }
)



#_________________________________________________________________________
#////////// COMPUTE INTERPOLATION MATRICES ///////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.getPhi","basisN",
  function( 
    OBJECT,           # basis OBJECT
    x = OBJECT@nodes, # vector of the evaluations points
    order = 0         # order of derivatives
  ){
    
    
    # Determine the problem dimension
    d <- OBJECT@d
    order <- if (is.matrix(order)) order else matrix(order,nrow=1)       # convert order to matrix
    order <- if (ncol(order)==1) matrix(order,nrow(order),d) else order  # Expand ORDER if it has a single columns  
    if(ncol(order)!=d) stop('In basis.getPhi, class basisN: order must have d columns (if matrix) or elements (if vector)')
    
    # Check validity of input x
    if (hasArg(x)){
      if (!(class(x)%in%c('list','matrix'))) stop('In basis.getPhi,x must be either a matrix or a list')
      if (is.matrix(x) && ncol(x)!=d) stop('In basis.getPhi, class basisN: x must have d columns')
      if (is.list(x) && length(x)!=d) stop('In basis.getPhi, class basisN: x must have d elements')
    }
   
    # Compute interpolation matrices for each dimension
    
    Blist <- list()
    orderj <- list()
    
    for (j in 1:d){
      orderj[[j]] <- sort(unique(order[,j]))
      
      temp <- if (!hasArg(x)){
        basis.getPhi(OBJECT@bas[[j]],order=orderj[[j]]) # use default nodes from 1D basis
      } else {
        if (is.list(x)){
          basis.getPhi(OBJECT@bas[[j]],x[[j]],orderj[[j]])
        } else {
          basis.getPhi(OBJECT@bas[[j]],x[,j],orderj[[j]])
        }
      }
      
      Blist[[j]]  <- if (is.list(temp)) temp else list(temp)   
    }
    
    xIsInputMatrix <- (hasArg(x) && is.matrix(x))
    
    B <- list()
    basList <- list()
    
    for (k in 1:nrow(order)){
      # MULTIPLY individual basis

      for (h in 1:d) {    #collect one basis per dimension, depending on order of derivative
        bpos <- which(orderj[[h]]==order[k,h])
        basList[[h]] <- Blist[[h]][[bpos]]
      }
      B[[k]] <- basis.product(basList,xIsInputMatrix,(OBJECT@format=='complete'))
      if (length(x)==1) B[[k]] <- matrix(B[[k]],nrow=1) #otherwise the dimension is dropped!    
      
      #  Label columns and rows in output matrices
      orderSign <- sign(order[k,])  
      if (max(orderSign) - min(orderSign)!=2){ # only label if derivative and integrals are not mixed
        if (all(orderSign==0)){
          colPrefix <- ''
        } else {
          colPrefix <- if (any(orderSign==1)){
            paste('D',paste(order[k,],sep='',collapse=''),'',sep='')
          } else {
            paste('I',paste(order[k,],sep='',collapse=''),'',sep='')
          }
        }
        colnames(B[[k]]) <- paste(colPrefix,colnames(B[[k]]),sep='')
      }
    }
    return(if (length(B)==1) B[[1]] else B)
  }
)

#_________________________________________________________________________
#////////// EVALUATE FUNCTION WITH BASIS// ///////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.evaluate","basisN",
  function( 
    OBJECT,           # basis OBJECT
    x,                # vector of the evaluations points
    coef,              # coefficient matrix
    order = 0        # order of derivatives
  ){
    
   
    # Get interpolation matrix
    B <- basis.getPhi(OBJECT,x,order)
    
    
    y <- if (is.list(B)) lapply(B,function(b) b%*% coef) else B %*% coef  
    return(y)
  }
)
