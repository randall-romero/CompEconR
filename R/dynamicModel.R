#============================================================
#     CLASS TO CONTAIN MODEL FOR DPSOLVE
#============================================================


setClass("dynamic.model",
         representation(
           description = "character", # (optional) description of the model
           horizon = "numeric",    # time horizon (infinite)
           func    = 'function',   # name of function file (required)
           params  = "list",       # model parameters required by function file (empty) 
           discount = "numeric",   # discount factor (required)
           ds = "numeric",         # dimension ds of the continuous state s (1)
           dx = "numeric",         # dimension dx of the continuous action x (0)
           ni = "numeric",         # number ni of discrete states i (none)
           nj = "numeric",         # number nj of discrete actions j (none)
           e  = "matrix",         # ne.de array of discretized continuous state shocks (0)
           w  = "numeric",         # ne.1 vector of discretized continuous state shock probabilities (1)
           q  = "array",           # ni.ni.nj array of stochastic discrete state transition probabilities (empty)
           h  = "matrix",         # nj.ni array of deterministic discrete state transitions (empty)
           X  = "matrix"          # nx.dx array of discretized continuous actions (empty)
         ),
         prototype(
           horizon   = Inf,
           ds        = 1,
           dx        = 1,
           ni        = 1,
           nj        = 1,
           e         = matrix(0,nrow=1,ncol=1),
           w         = 1,
           q         = array(1,c(1,1,1)),
           h         = matrix(1),
           X         = matrix(0)
         )
         )

setMethod("show","dynamic.model",
          function(object){
            SPRINTF <- function(...) cat(sprintf(...))
            cat('Dynamic Model',object@description,'\n')
            format <- '\t%36s : %1s\n'
            
            SPRINTF(format,'time horizon',object@horizon)
            SPRINTF(format,'discount factor', object@discount)
            SPRINTF(format,'dim. of continuous state, ds',object@ds)
            SPRINTF(format,'dim. of continuous action, dx',object@dx)
            if (object@ni>1) SPRINTF(format,'# of discrete states, ni',object@ni)
            if (object@nj>1) SPRINTF(format,'# of discrete actions, nj',object@nj)
            if (nrow(object@e)>1) SPRINTF(format,'continuous state shocks, e',paste(nrow(object@e),'X',ncol(object@e),'matrix'))
            if (length(object@w)>1) SPRINTF(format,'continuous state shocks prob, w',paste(length(object@w),'X 1 vector'))
            if (nrow(object@X)>1) SPRINTF(format,'discretized continuous actions, X',paste(nrow(object@X),'X',ncol(object@X),'matrix'))
            if (nrow(object@h)>1) SPRINTF(format,'deterministic discrete state transitions, h',paste(nrow(object@h),'X',ncol(object@h),'matrix'))
            if (prod(dim(object@q))>1) SPRINTF(format,'stochastic discrete state transition prob, q',paste(dim(object@q)[3],'matrices size ',dim(object@q)[1],'X',dim(object@q)[2]))
            
            nparams <- length(object@params)
            if (nparams>0){
              parnames <- names(object@params)
              cat('\n\t Parameters for func:\n')
              for (k in 1:nparams) SPRINTF('\t %12s: %1s\n',parnames[k],object@params[[k]])
            }
#             
#             q  = "array",           # ni.ni.nj array of abilities (empty)


          }
)










#============================================================
#     CLASS TO SET OPTIONS FOR DPSOLVE
#============================================================

setClass("dpsolve.options",
         representation(
           algorithm = "character",  # collocation equation solution algorithm, either Newton's method (`newton') or function iteration `funcit'
           tol       = "numeric",    # convergence tolerance (square root of machine epsilon)
           maxit     = "numeric",    # maximum number of iterations, collocation equation solution algorithm (500)
           nr        = "numeric",    # continuous state array refinement factor (10)
           maxitncp  = "numeric",    # maximum number of iterations, nonlinear complementarity solver (50)
           ncpmethod = "character",  # nonlinear complementarity solution algorithm for continuous action maximization problem embedded in Bellman
                                     # equation, either semi-smooth formulation ('ssmooth') or min-max formulation 'minmax'
           output    = "logical"     # whether to generate printed output (TRUE)
           ),
         prototype(
           algorithm = "newton",
           tol       = sqrt(.Machine$double.eps),
           maxit     = 500,
           nr        = 10,
           maxitncp  = 50,
           ncpmethod = "minmax",
           output    = TRUE
           )
)


