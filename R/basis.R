#===========================================================================
# DEFINE A GENERIC ONE DIMENSIONAL BASIS

setClass("basis",
         representation(type   = "character",  # type of basis
                        params = "list",       # parameters
                        n      = "numeric",    # number of nodes
                        a      = "numeric",    # lower limits
                        b      = "numeric",    # upper limits
                        nodes  = "matrix",     # nodes
                        D      = "list",       # operators to compute derivative
                        I      = "list",       # operators to compute integrals
                        varname = "character"  # optional name for variable
         ),
         validity = function(object){
           typeProvided <- substr(tolower(object@type),1,4)
           badTypes <- !(typeProvided %in% c("cheb","spli","lin"))

           if (any(badTypes)){
               stop(paste("Unknown type, valid options are",paste(validTypes,collapse = ", ")))
             } else {}
           return(TRUE)
           }

)


#=========================================================================
# METHODS FOR BASIS, TO BE DEFINE FOR EACH SPECIFIC TYPE
basis.define   <- function(OBJECT,...){} # defines n,a,b
basis.setnodes <- function(OBJECT,...){} # computes nodes
basis.diffoptr <- function(OBJECT,...){} # computes differential operator matrices
basis.getPhi   <- function(OBJECT,...){} # computes interpolation matrix
basis.evaluate <- function(OBJECT,...){} # evaluates multivariate functions with linear bases

#=========================================================================
# SAME METHODS, BUT BASED ON POINTERS (useful to keep basis updated)

setMethod(
  'basis.define','ref',
  function(OBJECT,...){
    ref2OBJECT <- OBJECT # save pointer in ref2OBJECT
    newObject <- basis.define(ref::deref(OBJECT),...)
    ref::deref(ref2OBJECT) <- newObject
    return(newObject)
  }
)

setMethod(
  'basis.setnodes','ref',
  function(OBJECT,...){
    ref2OBJECT <- OBJECT # save pointer in ref2OBJECT
    newObject <- basis.setnodes(ref::deref(OBJECT),...)
    ref::deref(ref2OBJECT) <- newObject
    return(newObject)
  }
)


setMethod(
  'basis.diffoptr','ref',
  function(OBJECT,...){
    ref2OBJECT <- OBJECT # save pointer in ref2OBJECT
    newObject <- basis.diffoptr(ref::deref(OBJECT),...)
    ref::deref(ref2OBJECT) <- newObject
    return(newObject)
  }
)


setMethod(
  'basis.getPhi','ref',
  function(OBJECT,...){
    ref2OBJECT <- OBJECT # save pointer in ref2OBJECT
    newObject <- basis.getPhi(ref::deref(OBJECT),...)
    ref::deref(ref2OBJECT) <- newObject
    return(newObject)
  }
)




#==============================================================================================
# FUNCTION TO INITIALIZE A ONE-DIMENSIONAL BASIS
basis <- function(type,...){
 # params <- as.list(match.call())
  type <- paste('basis',toupper(substr(type,1,1)),tolower(substr(type,2,4)),sep='')
  object <- new(type)             # create new basis object
  basis.define(ref::as.ref(object),...)        # defines the basis
  basis.setnodes(ref::as.ref(object))      # adds nodes to basis
  return(object)
}




setMethod(
  'basis.define','ref::ref',
  function(OBJECT,...){
    ref2OBJECT <- OBJECT # save pointer in ref2OBJECT
    newObject <- basis.define(ref::deref(OBJECT),...)
    ref::deref(ref2OBJECT) <- newObject
    return(newObject)
  }
)

setMethod(
  'basis.getPhi','ref::ref',
  function(OBJECT,...){
    ref2OBJECT <- OBJECT # save pointer in ref2OBJECT
    newObject <- basis.getPhi(ref::deref(OBJECT),...)
    ref::deref(ref2OBJECT) <- newObject
    return(newObject)
  }
)



#===============================================================================================
# METHOD TO SHOW BASIS
setMethod(
  "show","basis",
  function(object){
    strf <- '\t%10s : %1s\n'
    longtype <- switch(object@type,
                       cheb = 'Chebychev',
                       spli = paste('Spline(',object@k,')',sep=''),
                       lin = 'Linear')
    cat("Basis of 1 dimension\n")
    cat(sprintf(strf,'Type',longtype))
    cat(sprintf(strf,'# nodes',object@n))
    cat(sprintf(strf,'Interval',paste('[',object@a,', ',object@b,']',sep='')))
  }
)

