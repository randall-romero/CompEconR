#_________________________________________________________________________
#////////// MULTIPLY LIST OF BASES////////////////////////////////////////
#_________________________________________________________________________
setMethod(
  "basis.product","list",
  function(
    OBJECT,          # list of interpolation bases
    expandedNodes = TRUE, # whether rows of each basis have all desired node combinations
    completePoly = FALSE  # wheter resulting polynomials should be truncated at max(degree_i)
  ){
    
    if(completePoly){ # COMPUTE COMPLETE BASIS
      B <- dprodcomplete(OBJECT,expandedNodes)
      
    } else{  # COMPUTE TENSOR BASIS
      d <- length(OBJECT)  # basis dimension
      # multiply basis//////
      B <- OBJECT[[d]]
      for (h in (d-1):1) {        
        B <- if(expandedNodes){
          dprod(temp1,OBJECT[[h]])
        } else {
          kronecker(B,OBJECT[[h]])
        }
      }
      
      # get labels ///////
      degs <- list()  # degree of polynomials 
      for (h in 1:d) degs[[h]] <- 1:ncol(OBJECT[[h]])
      degs <- gridmake(degs)
      degs <- apply(degs,MARGIN=1,paste,sep='',collapse='')
      
      colnames(B) <- paste('Phi',degs,sep='')
      rownames(B) <- paste('x',if (expandedNodes) 1:nrow(B) else degs)      
    }
    return(B)
  }
  
)