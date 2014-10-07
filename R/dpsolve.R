#'  General Discrete-Time Bellman Equation Solver
#'
#' DPSOLVE uses the method of collocation to solve finite- and infinite-horizon
#' discrete-time stochastic dynamic optimization models with discrete, continuous,
#' or mixed states and actions of arbitrary dimension.  Usage of DPSOLVE is fully
#' documented in the pdf file:
#'
#' @author Randall Romero-Aguilar, based on Miranda & Fackler's CompEcon toolbox
#' @references Miranda, Fackler 2002 Applied Computational Economics and Finance
#'
#' @example demo/dummy_model_for_testing.r

dpsolve <- function(
  model,
  basis,
  guess = list(x = array(0,c(prod(basis@n),model@dx,model@ni,model@nj)),
               v = array(0,c(prod(basis@n),model@ni,model@nj))),
  options = new(Class="dpsolve.options")
){

  #  General Discrete-Time Bellman Equation Solver
  #
  # DPSOLVE uses the method of collocation to solve finite- and infinite-horizon
  # discrete-time stochastic dynamic optimization models with discrete, continuous,
  # or mixed states and actions of arbitrary dimension.  Usage of DPSOLVE is fully
  # documented in the pdf file:


  # '[' <- function(...) base::'['(...,drop=FALSE) # to prevent matrices turning into vectors



  vmax <- function(
    s, # continuous state
    x, # continuous action
    c  # approximation coefficients
  ){

    delta <- model@discount

    # Determine number of continuous state nodes
    ns <- nrow(s)             # number of continuous state nodes
    nx <- nrow(model@X)             # number of discretized continuous actions

    # Compute Bellman equation optimand - discrete or discretized continuous action
    if (dx==0 || nx>1){
      v <- array(0,c(ns,ni,nj))
      x <- array(0,c(ns,dx,ni,nj))
      for (i in 1:ni){
        for (j in 1:nj){
          vv <- matrix(0,ns,nx)
          xbound <- if (nx>1) func(flag='b',s=s,i=i,j=j,params=params)

          for (ix in 1:nx){
            vv[,ix] <- -Inf
            xx <- model@X[rep(ix,ns),,drop=FALSE]

            is <- if(nx>1) which(apply(xx>=xbound$lower & xx<=xbound$upper,MARGIN=1,all)) else 1:ns

            if (length(is)>0){
              # Initialize optimand and derivatives with reward terms
              vv[is,ix] <- func(flag='f',s=s[is,],x=xx[is,,drop=FALSE],i=i,j=j,params=params)$f
              for (k in 1:length(w)){
                # Compute states next period and derivatives
                ee <- e[rep(k,length(is)),]
                snext <- func(flag='g',s=s[is,],x=xx[is,,drop=FALSE],i=i,j=j,e=ee,params=params)$g
                snext <- Re(snext)
                for (iin in (1:ni)){
                  if (q[i,iin,j]==0) next
                  prob <- w[k]*q[i,iin,j]
                  vn <- basis.evaluate(basis,snext,c[,iin])
                  vv[is,ix] <- vv[is,ix] + delta*prob*vn
                }
              }
            }
          }
          v[,i,j] <- apply(vv,MARGIN=1,max)
          x[,,i,j] <- model@X[apply(vv,MARGIN=1,which.max),]
        }
      }

      #return(list(v=v,x=x))
    }

    # Compute Bellman equation optimand - continuous or mixed action
    if (dx>0 && nx<2){
      x <- array(x,dim=c(ns,dx,ni,nj))
      v <- array(0,c(ns,ni,nj))
      for (i in 1:ni){         #===for each discrete state
        for (j in 1:nj){       #===for each discrete action
          xbound <- func(flag='b',s=s,i=i,j=j,params=params)
          for (it in 1:options@maxitncp){    # while below maximum number or iterations
            # Initialize optimand and derivatives with reward terms
            xpolicy <- x[,,i,j,drop=FALSE]
            dim(xpolicy) <- dim(xpolicy)[1:2]

            reward <- func(flag='f',
                           s=s,
                           x=xpolicy,
                           i=i,j=j,
                           params=params)  # [vv,vx,vxx]
            vv  <- reward$f
            vx  <- reward$fx
            vxx <- reward$fxx



            for (k in 1:length(w)){  #===for each discretized shock
              # Compute states next period and derivatives
              ee <- e[rep(k,ns),]
              transition <- func(flag='g',
                                 s=s,
                                 x=x[,,i,j,drop=FALSE],
                                 i=i,j=j,
                                 e=ee,
                                 params=params)
              snext <- Re(transition$g)
              snx <- transition$gx
              snxx <- transition$gxx
              vnss <- array(0,dim=c(ns,ds,ds))   # Hessian matrix over all posible states
              for (iin in 1:ni){  #===for each FUTURE discrete state
                if (q[i,iin,j]==0) next
                prob <- w[k] * q[i,iin,j]

                vn  <- basis.evaluate(basis,coef=c[,iin,drop=FALSE])   #value function next period
                vns <- basis.evaluate(basis,coef=c[,iin,drop=FALSE],order=diag(ds)) #its first derivatives
                for (is in 1:ds){
                  for (js in is:ds){
                    order <- rep(0,ds)
                    order[is] <- order[is] + 1
                    order[js] <- order[js] + 1
                    vnss[,is,js] <- as.vector(basis.evaluate(basis,coef=c[,iin],order=order))
                    vnss[,js,is] <- vnss[,is,js]
                  }
                }
                vv <- vv + (delta*prob)*vn
                for (ix in 1:dx){
                  for (is in 1:ds){
                    vx[,ix] <- vx[,ix] + (delta*prob) * (vns[,is] * snx[,is,ix])
                    for (jx in 1:dx){
                      vxx[,ix,jx] <- vxx[,ix,jx] + (delta*prob) * (vns[,is] * snxx[,is,ix,jx])
                      for (js in 1:ds){
                        vxx[,ix,jx] <- vxx[,ix,jx] + (delta*prob) * vnss[,is,js] * snx[,is,ix] * snx[,js,jx]
                      }
                    }
                  }
                }
              }
            }
            # Compute Newton step, update continuous action, check convergence
            temp <-  lcpstep(options@ncpmethod,x[,,i,j,drop=FALSE],xbound$lower,xbound$upper,vx,vxx)
            tempx <- x[,,i,j,drop=FALSE]
            dim(tempx) <- dim(tempx)[1:2]
            x[,,i,j] <- tempx + temp$dx
            if (norm(temp$F,'I') < options@tol) break
          }
          v[,i,j] <- vv
        }
      }
    }



    if (options@algorithm=='funcit'){
      vc = NULL
    } else {

      #
      # if nargout<3, return, }
      #
      # Compute derivative of Bellman equation optimand with respect to basis
      # coefficients for Newton method, if requested

      if (ni*nj>1){
        vc <- array(0,c(ns,ni,ns,ni))
        jmax <- apply(v,c(1,2),which.max)
        for (i in 1:ni){
          for (j in 1:nj){
            is <- which(jmax[,i] == j)
            if (length(is)==0) next
            for (k in 1:length(w)){
              ee <- matrix(e[k,],ns,ncol(e),byrow=TRUE) #   e(k+zeros(ns,1),:);
              snext <- func('g',s[is,],x[is,,i,j,drop=FALSE],i,j,ee[is,],params=params);
              B <-  basis.getPhi(basis,snext)
              for (iin in 1:ni){
                if (q[i,iin,j]==0) next
                prob <- w[k]*q[i,iin,j]
                vc[is,i,,iin] <- vc[is,i,,iin] + delta*prob* matrix(B,c(length(is),1,ns))
              }
            }
          }
        }
        vc <- matrix(vc,c(ni*ns,ni*ns))
      } else {
        vc <- matrix(0,ns,ns)
        for (k in 1:length(w)){
          ee <- matrix(e[k,],ns,ncol(e),byrow=TRUE) #   e(k+zeros(ns,1),:);
          snext <- func('g',s,x,1,1,ee,params=params)$g
          vc <- vc + delta*w[k]*basis.getPhi(basis,snext)
        }
      }
    }
    return(list(v=v,x=x,vc=vc))

  }





















  inputByReference <- (class(basis) == 'ref')

  if (inputByReference){
    ref2basis <- basis
    basis <-  deref(basis)
  }



  # Unpack model structure
  T      <- model@horizon
  func   <- model@func
  params <- model@params
  ds     <- model@ds
  dx     <- model@dx
  ni     <- model@ni
  nj     <- model@nj
  e      <- model@e
  w      <- model@w
  q      <- array(model@q,c(ni,ni,nj))
  h      <- matrix(model@h,nj,ni)

  # Equivalent transition probability matrix for deterministic discrete state
  if (ni==1) h <- matrix(1,nj,1)
  if (is.nan(q)){
    if (is.nan(h)){
      if (ni==nj){
        warning('Neither q nor h specified; will default to h(j,i)=j. Please specify q or h if this is not correct.')
        h <- matrix(1:nj,ncol<-1) %*% matrix(1,1,ni)
      } else {
        stop('Either q or h must be specified.')
      }
    }
    for (i in 1:ni){
      for (j in 1:nj) {q[i,h[j,i],j] <- 1 }
    }
  }
  model@q <- q
  model@h <- h

  # Determine number of continuous state collocation nodes
  n  <- basis@n          # number of continuous state collocation node coordinates by dimension
  nc <- prod(n)          # number of continuous state collocation nodes


  # Set initial guess for values and actions if not provided
  x  <- guess$x
  v <- guess$v



  # Compute continuous state collocation nodes, collocation matrix and basis coefficients
  s   <- basis@nodes               # continuous state collocation nodes
  Phi <- basis.getPhi(basis)        # collocation matrix
  maxv <- apply(guess$v,c(1,2),max)        # max V over all discrete actions
  c <- solve(Phi,maxv)                     # basis coefficients



  # Solve finite horizon model by backward recursion.
  if (T<Inf){
    if (options@output) cat('Solving finite horizon model by backward recursion.\n')

    cc <- array(0,dim=c(nc,ni,T+2))
    vv <- array(0,dim=c(nc,ni,nj,T+2))
    xx <- array(0,dim=c(nc,dx,ni,nj,T+1))
    cc[,,T+2] <- c
    vv[,,,T+2] <- guess$v

    for (t in (T+1):1){
      VMAX <- vmax(s,x,c)        # update optimal values and actions
      x <- VMAX$x
      maxv <- apply(VMAX$v,c(1,2),max)       # max v over all discrete actions
      c <- solve(Phi,maxv)                    # update basis coefficients
      cc[,,t] <- c
      vv[,,,t] <- VMAX$v
      xx[,,,,t] <- x
    }
    return(list(c=drop(cc),s=s,v=drop(vv),x=drop(xx)))

  }


  # Solve infinite-horizon model by function iteration.
  if (T==Inf && options@algorithm=='funcit'){
    if (options@output) cat('Solving infinite-horizon model by function iteration.\n')
    for (it in 1:options@maxit){
      cold <- c                                           # store old basis coefficients
      VMAX <- vmax(s,x,c)                    # update optimal contingent values and actions
      x <- VMAX$x
      maxv <- apply(VMAX$v,c(1,2),max)       # max v over all discrete actions
      c <- solve(Phi,maxv)                    # update basis coefficients
      change <- norm(matrix(c-cold),"I")                    # compute change
      if (options@output) cat(sprintf ('%4i %10.1e\n',it,change))  # print change
      if (change<options@tol) break                         # convergence check
    }
    return(list(c=drop(c),s=s,v=drop(VMAX$v),x=drop(x)))
  }


  #Solve infinite horizon model by Newton's method.
  if (T==Inf && options@algorithm=='newton'){
    if (options@output) cat("Solving infinite horizon model by Newton's method.\n")
    Phik <- kronecker(diag(ni),Phi)
    for (it in 1:options@maxit){
      cold <- c                                           # store old basis coefficients
      VMAX  <- vmax(s,x,c)             # update optimal contingent values and actions
      x <- VMAX$x
      vm <- apply(VMAX$x,c(1,2),max)   # compute optimal value
      c <- vec(c) - solve(Phik - VMAX$vc, Phik %*% vec(c)-vec(vm))  # update basis coefficients
      c <- matrix(c,nc,ni)                               # reshape basis coefficient array
      change <- norm(vec(c- cold),'I');                    # compute change
      if (options@output) cat(sprintf('%4i %10.1e\n',it,change))  # print change
      if (change<options@tol) break                            # convergence check
    }
    return(list(c=drop(c),s=s,v=drop(VMAX$v),x=drop(x)))
  }































  if (inputByReference) deref(ref2basis) <- basis  # Update Object in calling environment
  return("SO FAR SO GOOD")
}


















































#
#                          # Solve finite horizon model by backward recursion.
#                          if T<inf
#                          if output, display('Solving finite horizon model by backward recursion.'), }
#                          cc <- zeros(nc,ni,T+2);
#                          vv <- zeros(nc,ni,nj,T+2);
#                          xx <- zeros(nc,dx,ni,nj,T+1);
#                          cc(:,:,T+2) <- c;
#                          vv(:,:,:,T+2) <- v;
#                          for t<-T+1:-1:1
#                          [v,x] <- vmax(model,basis,s,x,c);        # update optimal values and actions
#                          c <- Phi\max(v,[],3);                    # update basis coefficients
#                          cc(:,:,t) <- c;
#                          vv(:,:,:,t) <- v;
#                          xx(:,:,:,:,t) <- x;
#                          }
#                          sr <- s;
#                          c  <- squeeze(cc);
#                          vr <- squeeze(vv);
#                          xr <- squeeze(xx);
#                          resid <- [];
#                          return
#                          }
#
#                          tic
#
#                          # Solve infinite-horizon model by function iteration.
#                          if T==inf && strcmp(algorithm,'funcit')
#                          if output, display('Solving infinite-horizon model by function iteration.'), }
#                          for it<-1:maxit
#                          cold <- c;                                           # store old basis coefficients
#                          [v,x] <- vmax(model,basis,s,x,c);                    # update optimal contingent values and actions
#                          c <- Phi\max(v,[],3);                                # update basis coefficients
#                          change <- norm(c(:)-cold(:),inf);                    # compute change
#                          if output, fprintf ('#4i #10.1e\n',it,change), }	# print change
#                          if change<tol, break, };                          # convergence check
#                          }
#                          }
#
#                          # Solve infinite horizon model by Newton's method.
#                          if T==inf && strcmp(algorithm,'newton')
#                          if output, display('Solving infinite horizon model by Newton''s method.'), }
#                          Phik <- kron(eye(ni),Phi);
#                          for it<-1:maxit
#                          cold <- c;                                           # store old basis coefficients
#                          [v,x,vc] <- vmax(model,basis,s,x,c);                 # update optimal contingent values and actions
#                          vm <- max(v,[],3);                                   # compute optimal value
#                          c <- c(:)-(Phik-vc)\(Phik*c(:)-vm(:));               # update basis coefficients
#                          c <- reshape(c,nc,ni);                               # reshape basis coefficient array
#                          change <- norm(c(:)-cold(:),inf);                    # compute change
#                          if output, fprintf ('#4i #10.1e\n',it,change), }	# print change
#                          if change<tol, break, };                          # convergence check
#                          }
#                          }
#
#                          if output
#                          if it==maxit, display('Failure to converge in dpsolve.'), }
#                          fprintf('Elapsed Time <- #7.2f Seconds\n',toc)
#                          }
#
#                          # Check whether continuous state transitions remain within approximation interval
#                          snmin <-  inf;
#                          snmax <- -inf;
#                          for i<-1:ni
#                          [~,jmax] <- max(v(:,i,:),[],3);
#                          for j<-1:nj
#                          is <- find(jmax==j);
#                          if isempty(is), continue, }
#                          ns <- length(is);
#                          for k<-1:size(e,1);
#                          ee <- e(k+zeros(ns,1),:);
#                          sn <- feval(func,'g',s(is,:),x(is,:,i,j),i,j,ee,params{:});
#                          snmin <- min(snmin,min(sn));
#                          snmax <- max(snmax,max(sn));
#                          }
#                          }
#                          }
#                          if output
#                          if any(snmin<basis.a-eps), display('Extrapolating below smin. Decrease smin.'), };
#                          if any(snmax>basis.b+eps), display('Extrapolating above smax. Increase smax.'), };
#                          disp([basis.a' snmin' snmax' basis.b'])
#                          }
#
#                          # Compute values and actions on refined continuous state array, if requested.
#                          if nr>0
#                          n  <- nr*n;
#                          ns <- prod(n);
#                          for is<-1:ds
#                          srcoord{is} <- linspace(basis.a(is),basis.b(is),n(is))';
#                          }
#                          sr <- gridmake(srcoord);
#                          if dx>0
#                          xr <- zeros(ns,dx,ni,nj);
#                          for i<-1:ni
#                          for j<-1:nj
#                          xr(:,:,i,j) <- funeval(Phi\x(:,:,i,j),basis,sr);
#                          }
#                          }
#                          else
#                            xr <- [];
#                          }
#                          [vr,xr] <- vmax(model,basis,sr,xr,c);
#                          else
#                            sr <- s;
#                          vr <- v;
#                          xr <- x;
#                          }
#
#                          # Compute Bellman equation residual on refined continuous state array, if requested
#                          if nargout>4
#                          vrproj <- funeval(c,basis,sr);
#                          resid <- vrproj - max(vr,[],3);
#                          resid <- squeeze(resid);
#                          }
#
#                          # Eliminate singleton dimensions to facilitate analysis in calling program
#                          c  <- squeeze(c);
#                          vr <- squeeze(vr);
#                          xr <- squeeze(xr);
#
#
#                          ## VMAX
#                          #
#                          #   Function called by dpsolve to maximize the optimand embedded in the
#                          #   Bellman equation with respect to the continuous action x at ns input
#                          #   continuous state nodes s, given a single discrete state i and a single
#                          #   discrete action j.  Does so by solving associated Karush-Kuhn-Tucker
#                          #   necessary conditions as a nonlinear complementarity problem, using a
#                          #   derivative-based adapted Newton method (see Miranda and Fackler).
#                          #   Function optionally also computes derivative of the optimand with
#                          #   respect to the value function approximant basis function coefficients,
#                          #   which is required to solve the collocation equation using Newton's
#                          #   method.
#                          #
#                          # Usage
#                          #   [v,x,vc] <- vmax(model,basis,s,x,c)
#                          # Let
#                          #   n  <- number of basis functions and continuous state collocation nodes
#                          #   ns <- number of continuous state nodes on refined output array
#                          #   ds <- dimension of continuous state s
#                          #   dx <- dimension of continuous action x
#                          #   de <- dimension of continuous state transition shock e
#                          #   ni <- number of discrete states i
#                          #   nj <- number of discrete actions j
#                          #   nx <- number of discretized continuous actions, if applicable
#                          # Input
#                          #   model  : structured array containing model specifications
#                          #   basis  : n-dimensional basis for real-valued functions defined on the
#                          #            continuous state space
#                          #   s      : ns.ds array of continuous state nodes
#                          #   x      : n.dx.ni.nj array of initial guesses for optimal continuous
#                          #            actions at the ns continuous state nodes, per discrete state
#                          #            and discrete action
#                          #   c      : n.ni vector of value function approximant basis function
#                          #            coefficients, per discrete state
#                          # Output
#                          #   v      : n.ni.nj array of discrete-action-contingent values at the ns
#                          #            continuous state nodes, per discrete state and discrete action
#                          #   x      : ns.dx.ni.nj array of optimal continuous actions at the ns
#                          #            continuous state nodes, per discrete state and discrete action
#                          #   vc     : nn.nn (nn<-nc*ni*nj) array of derivatives of values with
#                          #            respect to basis coefficients
#
#


































# ## DPSOLVE
# #
# #  General Discrete-Time Bellman Equation Solver
# #
# #  DPSOLVE uses the method of collocation to solve finite- and infinite-
#   #  horizon discrete-time stochastic dynamic optimization models with
# #  discrete, continuous, or mixed states and actions of arbitrary
# #  dimension.  Usage of DPSOLVE is fully documented in the pdf file
# #  `DPSOLVE - A MATLAB Routine for Solving Discrete-Time Bellman Equations'
# #  included in the directory CompEcon2010 accompanying the distribution of
# #  this file.  A summary of its features is provided here:
# #
# #  Usage
# #    [c,sr,vr,xr,resid] = dpsolve(model,basis,v,x)
# #  Let
# #    n  = number of basis functions and continuous state collocation nodes
# #    ns = number of continuous state nodes on refined output array
# #    ds = dimension of continuous state s
# #    dx = dimension of continuous action x
# #    de = dimension of continuous state transition shock e
# #    ni = number of discrete states i
# #    nj = number of discrete actions j
# #    nx = number of discretized continuous actions, if applicable
# #  Input
# #    model  : structured array containing model specifications (see below)
# #    basis  : n-dimensional basis for real-valued functions defined on the
# #             continuous state space
# #    v      : n.ni.nj array of initial guesses for discrete-action-
# #             contingent value function values at the n continuous state
# #             collocation nodes, per discrete state and discrete action
# #             (optional, default is array of zeros)
# #    x      : n.dx.ni.nj array of initial guesses for optimal continuous
# #             actions at the n state collocation nodes, per discrete state
# #             and discrete action (optional, default is array of zeros)
# #  Output
# #    c      : n.ni vector of value function approximant basis function
# #             coefficients, per discrete state
# #    sr     : ns.ds refined array of continuous state nodes
# #    vr     : ns.ni.nj array of discrete-action-contingent values on the
# #             refined array of continuous state nodes, per discrete state
# #             and discrete action
# #    xr     : ns.dx.ni.nj array of discrete-action-contingent optimal
# #             continuous actions on the refined array of continuous state
# #             nodes, per discrete state and discrete action
# #    resid  : ns.ni array of Bellman equation residuals on the refined
# #             array of continuous state nodes, per discrete state
# #    Note: If the residual is requested, values, optimal continuous
# #    actions, and residuals are returned on a refined array of ns
# #    equally-spaced continuous state nodes. The degree of refinement is
# #    governed by nr, an optional parameter with default value of 10 that
# #    may be set by the user using optset (see below).  If nr>0, the refined
# #    continuous state array is created by forming the Cartesian product of
# #    equally-spaced coordinates along each dimension, with nr times the
# #    number of coordinates possessed by the continuous state collocation
# #    array along each dimension. If nr=0, the routine will return the
# #    values, optimal continuous actions, and residuals on the original
# #    array of continuous state collocation nodes.
# #    Note: The singleton dimensions of c, vr, xr, and resid are eliminated
# #    to facilitate analysis in the calling program.
# #  Model Structure
# #    The structured array "model" contains different fields that specify
# #    essential features of the model to be solved (default values for
# #    optional fields in parentheses):
# #      horizon   : time horizon (infinite)
# #      func      : name of function file (required)
# #      params    : model parameters required by function file (empty)
# #      discount  : discount factor (required)
# #      ds        : dimension ds of the continuous state s (1)
# #      dx        : dimension dx of the continuous action x (0)
# #      ni        : number ni of discrete states i (none)
# #      nj        : number nj of discrete actions j (none)
# #      e         : ne.de array of discretized continuous state shocks (0)
# #      w         : ne.1 vector of discretized continuous state shock
# #                  probabilities (1)
# #      q         : ni.ni.nj array of stochastic discrete state transition
# #                  probabilities (empty)
# #      h         : nj.ni array of deterministic discrete state transitions
# #                  (empty)
# #      X         : nx.dx array of discretized continuous actions (empty)
# #    Note: If X is empty (the default), the routine will attempt to solve
# #    the continuous action maximization problem embedded in the Bellman
# #    equation by solving the associated Karush-Kuhn-Tucker conditions as a
# #    nonlinear complementarity problem, using a derivative-based adapted
# #    Newton method (see Miranda and Fackler).
# #    Note: If the stochastic discrete state transition probabilities q are
# #    specified, then the deterministic state transitions h should not be
# #    specified, and vice versa.
# #  Function File
# #    User-supplied function that returns the bound, reward, and continuous
# #    state transition functions and needed derivatives with respect to the
# #    continuous action x at an arbitrary number ns of continuous states s
# #    and continuous actions x, given a discrete state i and a discrete
# #    action j, according to the format
# #      [out1,out2,out3] = func(flag,s,x,i,j,e,<params>)
# #    Function File Input
# #      flag      : flag that specifies the desired function
# #      s         : ns.ds array of continuous state nodes
# #      x         : ns.dx array of continuous actions
# #      i         : a single discrete state index, an integer between 1-ni
# #      j         : a single discrete action index, an integer between 1-nj
# #      e         : ns.de array of continuous state transition shocks
# #      params    : user-supplied list of function parameters
# #    Function File Output
# #      flag = 'b' returns lower and upper bounds on continuous action x
# #        out1    : ns.dx lower bounds on continuous action x
# #        out2    : ns.dx upper bounds on continuous action x
# #        out3    : empty
# #      flag = 'f' returns reward and derivatives
# #        out1    : ns.1 reward function f values
# #        out2    : ns.dx first derivative of f with respect to x
# #        out3    : ns.dx.dx second derivative of f with respect to x
# #      flag = 'g' returns continuous state transition and derivatives
# #        out1    : ns.ds continuous state transition function g values
# #        out2    : ns.ds.dx first derivative of g with respect to x
# #        out3    : ns.ds.dx.dx second derivative of g with respect to x
# #      Note: If the continuous action is discretized, then the function
# #      file need only provide the reward and transition, with the reward
# #      set to negative infinity at non-feasible values of x as in
# #      flag = 'f' returns reward and derivatives
# #        out     : ns.1 reward f
# #      flag = 'g' returns continuous state transition and derivatives
# #        out     : ns.ds continuous state transition g
# #  Options
# #    algorithm : collocation equation solution algorithm, either Newton's
# #                method (`newton') or function iteration `funcit'
#                          #    ncpmethod : nonlinear complementarity solution algorithm for
#                          #                continuous action maximization problem embedded in Bellman
#                          #                equation, either semi-smooth formulation ('ssmooth') or
#                          #                min-max formulation 'minmax'
#                          #    maxit     : maximum number of iterations, collocation equation
#                          #                solution algorithm (500)
#                          #    maxitncp  : maximum number of iterations, nonlinear complementarity
#                          #                solver (50)
#                          #    tol       : convergence tolerance (square root of machine epsilon)
#                          #    nr        : continuous state array refinement factor (10)
#                          #    output    : whether to generate printed output (1)
#
#                          #  Copyright(c) 1997-2011
#                          #    Mario J. Miranda - miranda.4@osu.edu
#                          #    Paul L. Fackler  - paul_fackler@ncsu.edu
#
#                          function [c,sr,vr,xr,resid] = dpsolve(model,basis,v,x)
#
#                          # Set user options to defaults, if not set by user with OPTSET (see above)
#                          global dpsolve_options
#                          if ~isfield(dpsolve_options,'algorithm'), dpsolve_options.algorithm = 'newton'; end
#                          if ~isfield(dpsolve_options,'tol'),       dpsolve_options.tol = sqrt(eps);      end
#                          if ~isfield(dpsolve_options,'maxit'),     dpsolve_options.maxit = 500;          end
#                          if ~isfield(dpsolve_options,'nr'),        dpsolve_options.nr = 10;              end
#                          if ~isfield(dpsolve_options,'maxitncp'),  dpsolve_options.maxitncp = 50;        end
#                          if ~isfield(dpsolve_options,'ncpmethod'), dpsolve_options.ncpmethod = 'minmax'; end
#                          if ~isfield(dpsolve_options,'output'),    dpsolve_options.output = 1;           end
#
#                          # Set model fields to default if nonexistent (see above)
#                          if ~isfield(model,'horizon'),  model.horizon = inf;              end
#                          if ~isfield(model,'func'),     error('Missing Function File');   end
#                          if ~isfield(model,'params'),   error('Missing Parameter List');  end
#                          if ~isfield(model,'discount'), error('Missing Discount Factor'); end
#                          if ~isfield(model,'ds'),       model.ds = 1;                     end
#                          if ~isfield(model,'dx'),       model.dx = 1;                     end
#                          if ~isfield(model,'ni'),       model.ni = 1;                     end
#                          if ~isfield(model,'nj'),       model.nj = 1;                     end
#                          if ~isfield(model,'e'),        model.e  = 0;                     end
#                          if ~isfield(model,'w'),        model.w  = 1;                     end
#                          if ~isfield(model,'q'),        model.q  = [];                    end
#                          if ~isfield(model,'h'),        model.h  = [];                    end
#                          if ~isfield(model,'X'),        model.X  = zeros(1,model.dx);     end
#                          model.ni = max(1,model.ni);
#                          model.nj = max(1,model.nj);
#
#                          # Unpack user options
#                          algorithm = dpsolve_options.algorithm;
#                          tol       = dpsolve_options.tol;
#                          maxit     = dpsolve_options.maxit;
#                          output    = dpsolve_options.output;
#                          if nargout<5
#                          nr = 0;
#                          else
#                          nr = dpsolve_options.nr;
#                          end
#
#                          # Unpack model structure
#                          T      = model.horizon;
#                          func   = model.func;
#                          params = model.params;
#                          ds     = model.ds;
#                          dx     = model.dx;
#                          ni     = model.ni;
#                          nj     = model.nj;
#                          e      = model.e;
#                          q      = model.q;
#                          h      = model.h;
#
#                          # Equivalent transition probability matrix for deterministic discrete state
#                          if ni==1; h = ones(nj,1); end
#                          if isempty(q)
#                          if isempty(h)
#                          if ni==nj
#                          'Neither q nor h specified; will default to h(j,i)=j.'
#                          'Please specify q or h if this is not correct.'
#                          h = (1:nj)'*ones(1,ni);
#                          else
#                            error('Either q or h must be specified.')
#                          end
#                          end
#                          for i=1:ni
#                          for j=1:nj
#                          q(i,h(j,i),j) = 1;
#                          end
#                          end
#                          end
#                          model.q = q;
#                          model.h = [];
#
#                          # Determine number of continuous state collocation nodes
#                          n  = basis.n;               # number of continuous state collocation
#                          # node coordinates by dimension
#                          nc = prod(n);               # number of continuous state collocation nodes
#
#                          # Set initial guess for values and actions if not provided
#                          if nargin<3, v = zeros(nc,ni,nj);    else v = reshape(v,nc,ni,nj); end
#                          if nargin<4, x = zeros(nc,dx,ni,nj); else x = reshape(x,nc,dx,ni,nj); end
#
#                          # Compute continuous state collocation nodes, collocation matrix and basis coefficients
#                          s   = funnode(basis);        # continuous state collocation nodes
#                          Phi = funbas(basis,s);      # collocation matrix
#                          c = Phi\max(v,[],3);        # basis coefficients
#
#                          # Solve finite horizon model by backward recursion.
#                          if T<inf
#                          if output, display('Solving finite horizon model by backward recursion.'), end
#                          cc = zeros(nc,ni,T+2);
#                          vv = zeros(nc,ni,nj,T+2);
#                          xx = zeros(nc,dx,ni,nj,T+1);
#                          cc(:,:,T+2) = c;
#                          vv(:,:,:,T+2) = v;
#                          for t=T+1:-1:1
#                          [v,x] = vmax(model,basis,s,x,c);        # update optimal values and actions
#                          c = Phi\max(v,[],3);                    # update basis coefficients
#                          cc(:,:,t) = c;
#                          vv(:,:,:,t) = v;
#                          xx(:,:,:,:,t) = x;
#                          end
#                          sr = s;
#                          c  = squeeze(cc);
#                          vr = squeeze(vv);
#                          xr = squeeze(xx);
#                          resid = [];
#                          return
#                          end
#
#                          tic
#
#                          # Solve infinite-horizon model by function iteration.
#                          if T==inf && strcmp(algorithm,'funcit')
#                          if output, display('Solving infinite-horizon model by function iteration.'), end
#                          for it=1:maxit
#                          cold = c;                                           # store old basis coefficients
#                          [v,x] = vmax(model,basis,s,x,c);                    # update optimal contingent values and actions
#                          c = Phi\max(v,[],3);                                # update basis coefficients
#                          change = norm(c(:)-cold(:),inf);                    # compute change
#                          if output, fprintf ('#4i #10.1e\n',it,change), end	# print change
#                          if change<tol, break, end;                          # convergence check
#                          end
#                          end
#
#                          # Solve infinite horizon model by Newton's method.
#                          if T==inf && strcmp(algorithm,'newton')
#                          if output, display('Solving infinite horizon model by Newton''s method.'), end
#                          Phik = kron(eye(ni),Phi);
#                          for it=1:maxit
#                          cold = c;                                           # store old basis coefficients
#                          [v,x,vc] = vmax(model,basis,s,x,c);                 # update optimal contingent values and actions
#                          vm = max(v,[],3);                                   # compute optimal value
#                          c = c(:)-(Phik-vc)\(Phik*c(:)-vm(:));               # update basis coefficients
#                          c = reshape(c,nc,ni);                               # reshape basis coefficient array
#                          change = norm(c(:)-cold(:),inf);                    # compute change
#                          if output, fprintf ('#4i #10.1e\n',it,change), end	# print change
#                          if change<tol, break, end;                          # convergence check
#                          end
#                          end
#
#                          if output
#                          if it==maxit, display('Failure to converge in dpsolve.'), end
#                          fprintf('Elapsed Time = #7.2f Seconds\n',toc)
#                          end
#
#                          # Check whether continuous state transitions remain within approximation interval
#                          snmin =  inf;
#                          snmax = -inf;
#                          for i=1:ni
#                          [~,jmax] = max(v(:,i,:),[],3);
#                          for j=1:nj
#                          is = find(jmax==j);
#                          if isempty(is), continue, end
#                          ns = length(is);
#                          for k=1:size(e,1);
#                          ee = e(k+zeros(ns,1),:);
#                          sn = feval(func,'g',s(is,:),x(is,:,i,j),i,j,ee,params{:});
#                          snmin = min(snmin,min(sn));
#                          snmax = max(snmax,max(sn));
#                          end
#                          end
#                          end
#                          if output
#                          if any(snmin<basis.a-eps), display('Extrapolating below smin. Decrease smin.'), end;
#                          if any(snmax>basis.b+eps), display('Extrapolating above smax. Increase smax.'), end;
#                          disp([basis.a' snmin' snmax' basis.b'])
#                          end
#
#                          # Compute values and actions on refined continuous state array, if requested.
#                          if nr>0
#                          n  = nr*n;
#                          ns = prod(n);
#                          for is=1:ds
#                          srcoord{is} = linspace(basis.a(is),basis.b(is),n(is))';
#                          end
#                          sr = gridmake(srcoord);
#                          if dx>0
#                          xr = zeros(ns,dx,ni,nj);
#                          for i=1:ni
#                          for j=1:nj
#                          xr(:,:,i,j) = funeval(Phi\x(:,:,i,j),basis,sr);
#                          end
#                          end
#                          else
#                            xr = [];
#                          end
#                          [vr,xr] = vmax(model,basis,sr,xr,c);
#                          else
#                            sr = s;
#                          vr = v;
#                          xr = x;
#                          end
#
#                          # Compute Bellman equation residual on refined continuous state array, if requested
#                          if nargout>4
#                          vrproj = funeval(c,basis,sr);
#                          resid = vrproj - max(vr,[],3);
#                          resid = squeeze(resid);
#                          end
#
#                          # Eliminate singleton dimensions to facilitate analysis in calling program
#                          c  = squeeze(c);
#                          vr = squeeze(vr);
#                          xr = squeeze(xr);
#
#
#                          ## VMAX
#                          #
#                          #   Function called by dpsolve to maximize the optimand embedded in the
#                          #   Bellman equation with respect to the continuous action x at ns input
#                          #   continuous state nodes s, given a single discrete state i and a single
#                          #   discrete action j.  Does so by solving associated Karush-Kuhn-Tucker
#                          #   necessary conditions as a nonlinear complementarity problem, using a
#                          #   derivative-based adapted Newton method (see Miranda and Fackler).
#                          #   Function optionally also computes derivative of the optimand with
#                          #   respect to the value function approximant basis function coefficients,
#                          #   which is required to solve the collocation equation using Newton's
#                          #   method.
#                          #
#                          # Usage
#                          #   [v,x,vc] = vmax(model,basis,s,x,c)
#                          # Let
#                          #   n  = number of basis functions and continuous state collocation nodes
#                          #   ns = number of continuous state nodes on refined output array
#                          #   ds = dimension of continuous state s
#                          #   dx = dimension of continuous action x
#                          #   de = dimension of continuous state transition shock e
#                          #   ni = number of discrete states i
#                          #   nj = number of discrete actions j
#                          #   nx = number of discretized continuous actions, if applicable
#                          # Input
#                          #   model  : structured array containing model specifications
#                          #   basis  : n-dimensional basis for real-valued functions defined on the
#                          #            continuous state space
#                          #   s      : ns.ds array of continuous state nodes
#                          #   x      : n.dx.ni.nj array of initial guesses for optimal continuous
#                          #            actions at the ns continuous state nodes, per discrete state
#                          #            and discrete action
#                          #   c      : n.ni vector of value function approximant basis function
#                          #            coefficients, per discrete state
#                          # Output
#                          #   v      : n.ni.nj array of discrete-action-contingent values at the ns
#                          #            continuous state nodes, per discrete state and discrete action
#                          #   x      : ns.dx.ni.nj array of optimal continuous actions at the ns
#                          #            continuous state nodes, per discrete state and discrete action
#                          #   vc     : nn.nn (nn=nc*ni*nj) array of derivatives of values with
#                          #            respect to basis coefficients
#
#
#                          function [v,x,vc] = vmax(model,basis,s,x,c)
#
#                          # Unpack user options
#                          global dpsolve_options
#                          tol       = dpsolve_options.tol;
#                          maxitncp  = dpsolve_options.maxitncp;
#                          ncpmethod = dpsolve_options.ncpmethod;
#
#                          # Unpack model structure
#                          func   = model.func;
#                          params = model.params;
#                          delta  = model.discount;
#                          ds     = model.ds;
#                          dx     = model.dx;
#                          ni     = model.ni;
#                          nj     = model.nj;
#                          e      = model.e;
#                          w      = model.w;
#                          q      = model.q;
#                          X      = model.X;
#
#                          # Determine number of continuous state nodes
#                          ns = size(s,1);             # number of continuous state nodes
#                          nx = size(X,1);             # number of discretized continuous actions
#
#                          # Compute Bellman equation optimand - discrete or discretized continuous action
#                          if dx==0||nx>1
#                          v = zeros(ns,ni,nj);
#                          x = zeros(ns,dx,ni,nj);
#                          for i=1:ni
#                          for j=1:nj
#                          vv = zeros(ns,nx);
#                          if nx>1
#                          [xl,xu] = feval(func,'b',s,[],i,j,[],params{:});
#                          end
#                          for ix=1:nx
#                          vv(:,ix) = -inf;
#                          xx = X(ix+zeros(ns,1),:);
#                          if nx>1
#                          is = find(all(xx>=xl,2).* all(xx<=xu,2)==1);
#                          else
#                          is = 1:ns;
#                          end
#                          if ~isempty(is)
#                          # Initialize optimand and derivatives with reward terms
#                          vv(is,ix) = feval(func,'f',s(is,:),xx(is,:),i,j,[],params{:});
#                          for k=1:length(w)
#                          # Compute states next period and derivatives
#                          ee = e(k+zeros(length(is),1),:);
#                          snext = feval(func,'g',s(is,:),xx(is,:),i,j,ee,params{:});
#                          snext = real(snext);
#                          for in=1:ni
#                          if q(i,in,j)==0, continue, end
#                          prob = w(k)*q(i,in,j);
#                          vn = funeval(c(:,in),basis,snext);
#                          vv(is,ix) = vv(is,ix) + delta*prob*vn;
#                          end
#                          end
#                          end
#                          end
#                          [v(:,i,j),ix] = max(vv,[],2);
#                          x(:,:,i,j) = X(ix,:);
#                          end
#                          end
#                          end
#
#                          # Compute Bellman equation optimand - continuous or mixed action
#                          if dx>0&&nx<2
#                          x = reshape(x,ns,dx,ni,nj);
#                          v = zeros(ns,ni,nj);
#                          for i=1:ni
#                          for j=1:nj
#                          [xl,xu] = feval(func,'b',s,[],i,j,[],params{:});
#                          for it=1:maxitncp
#                          # Initialize optimand and derivatives with reward terms
#                          [vv,vx,vxx] = feval(func,'f',s,x(:,:,i,j),i,j,[],params{:});
#                          for k=1:length(w)
#                          # Compute states next period and derivatives
#                          ee = e(k+zeros(ns,1),:);
#                          [snext,snx,snxx] = feval(func,'g',s,x(:,:,i,j),i,j,ee,params{:});
#                          snext = real(snext);
#                          B = funbasx(basis,snext,[0;1;2]);
#                          vnss = zeros(ns,nj,ds,ds);
#                          for in=1:ni
#                          if q(i,in,j)==0, continue, end
#                          prob = w(k)*q(i,in,j);
#                          vn  = funeval(c(:,in),basis,B);
#                          vns = funeval(c(:,in),basis,B,eye(ds));
#                          for is=1:ds
#                          for js=is:ds
#                          order = zeros(1,ds);
#                          order(is) = order(is)+1;
#                          order(js) = order(js)+1;
#                          vnss(:,is,js) = funeval(c(:,in),basis,B,order);
#                          vnss(:,js,is) = vnss(:,is,js);
#                          end
#                          end
#                          vv = vv + delta*prob*vn;
#                          for ix=1:dx
#                          for is=1:ds
#                          vx(:,ix) = vx(:,ix) + delta*prob*vns(:,is).*snx(:,is,ix);
#                          for jx=1:dx
#                          vxx(:,ix,jx) = vxx(:,ix,jx) + delta*prob*vns(:,is).*snxx(:,is,ix,jx);
#                          for js=1:ds
#                          vxx(:,ix,jx) = vxx(:,ix,jx) + delta*prob*vnss(:,is,js).*snx(:,is,ix).*snx(:,js,jx);
#                          end
#                          end
#                          end
#                          end
#                          end
#                          end
#                          # Compute Newton step, update continuous action, check convergence
#                          [vx,delx] = lcpstep(ncpmethod,x(:,:,i,j),xl,xu,vx,vxx);
#                          x(:,:,i,j) = x(:,:,i,j) + delx;
#                          if norm(vx(:),inf)<tol, break, end;
#                          end
#                          v(:,i,j) = vv;
#                          end
#                          end
#                          end
#
#                          if nargout<3, return, end
#
#                          # Compute derivative of Bellman equation optimand with respect to basis
#                          # coefficients for Newton method, if requested
#
#                          if ni*nj>1
#                          vc = zeros(ns,ni,ns,ni);
#                          [~,jmax] = max(v,[],3);
#                          for i=1:ni
#                          for j=1:nj
#                          is = find(jmax(:,i)==j);
#                          if isempty(is), continue, end
#                          for k=1:length(w)
#                          ee = e(k+zeros(ns,1),:);
#                          snext = feval(func,'g',s(is,:),x(is,:,i,j),i,j,ee(is,:),params{:});
#                          B = full(funbas(basis,snext));
#                          for in=1:ni
#                          if q(i,in,j)==0, continue, end
#                          prob = w(k)*q(i,in,j);
#                          vc(is,i,:,in) = vc(is,i,:,in) + delta*prob*reshape(B,length(is),1,ns);
#                          end
#                          end
#                          end
#                          end
#                          vc = reshape(vc,ni*ns,ni*ns);
#                          else
#                          vc = zeros(ns,ns);
#                          for k=1:length(w)
#                          ee = e(k+zeros(ns,1),:);
#                          snext = feval(func,'g',s,x,1,1,ee,params{:});
#                          vc = vc + delta*w(k)*funbas(basis,snext);
#                          end
#                          end
