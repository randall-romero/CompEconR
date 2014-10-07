#==========CLASS DEFINITION===============================================
setClass("basisSpli",
         representation(
           breaks = "numeric", # break point for splines
           k = "numeric",      # order of spline polynomials
           evennum = "numeric"
         ),
         contains="basis"
)

#==========CLASS METHODS==================================================

#_________________________________________________________________________
#////////// DEFINE A NEW BASIS ///////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.define","basisSpli",
  function(OBJECT,breaks,evennum = 0,k=3){
    
    if (!hasArg(breaks))  stop("Parameter 'breaks' is required")
    if (length(k)>1)      stop('Spline order (k) should be a scalar')
    if (k<0)              stop('Spline order (k) is too small': ' num2str(k)]')
    if (length(breaks)<2) stop('breakpoint sequence must contain at least two elements')
    if (any(as.logical(diff(breaks)))<0) stop('Breakpoints must be non-decreasing')
    if (evennum==0){
      if (length(breaks)==2) evennum <- 2
    } else {
      if (length(breaks)==2) {
        breaks <- seq(breaks[1],breaks[2],len=evennum)
      } else {
        stop('Breakpoint sequence must contain 2 values when evennum>0')
      }
    }
    OBJECT@type <- 'spli'
    OBJECT@n <- length(breaks) + k - 1
    OBJECT@a <- breaks[1]
    OBJECT@b <- breaks[length(breaks)]
    OBJECT@breaks <- breaks
    OBJECT@k <- k
    OBJECT@evennum <- evennum
    return(OBJECT)
  }
)


#_________________________________________________________________________
#////////// COMPUTE NODES ////////////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  'basis.setnodes','basisSpli',
  function(OBJECT){
    
    n <- OBJECT@n
    a <- OBJECT@a
    b <- OBJECT@b
    k <- OBJECT@k
    
    x <- cumsum(rbind(matrix(a,k),matrix(OBJECT@breaks),matrix(b,k)))
    x <- (x[(1+k):(n+k)]-x[1:n])/k
    x[1] <- a
    x[length(x)] <- b
    
    OBJECT@nodes <- matrix(x,ncol=1)   
    return(OBJECT) 
  }
)



#_________________________________________________________________________
#////////// COMPUTE OPERATORS TO DIFFERENTIATE AND INTEGRATE//////////////
#_________________________________________________________________________
setMethod(
  "basis.diffoptr","basisSpli",
  function(
    OBJECT,    # basis or reference to a basis
    order = 1  # order of derivative
  ){
    
    
    
    
    
    
    return(OBJECT) 
  }
)


#_________________________________________________________________________
#////////// COMPUTE INTERPOLATION MATRICES ///////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.getPhi","basisSpli",
  function( 
    OBJECT,              # basis OBJECT
    x = OBJECT@nodes,    # vector of the evaluations points
    order = 0            # order of derivatives
  ){
    
    
    
    
    
    
    
    
    
    return(B)
  }
)