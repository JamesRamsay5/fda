fRegress <- function(y, xfdlist, betalist, wt=NULL, y2cMap=NULL, SigmaE=NULL,
                     method, sep='.', ...)
{
  
  #  FREGRESS  Fits a functional linear model using multiple
  #  functional independent variables with the dependency being
  #  pointwise or concurrent.
  #  The case of a scalar independent variable is included by treating
  #  it as a functional independent variable with a constant basis
  #  and a unit coefficient.
  #
  #  Arguments:
  #  Y        ... an object for the dependent variable,
  #               which may be:
  #                   a functional data object,
  #                   a functional parameter (fdPar) object, or
  #                   a vector
  #  XFDLIST  ... a list object of length p with each list
  #               containing an object for an independent variable.
  #               the object may be:
  #                   a functional data object or
  #                   a vector
  #               if XFDLIST is a functional data object or a vector,
  #               it is converted to a list of length 1.
  #  BETALIST ... a list object of length p with each list
  #               containing a functional parameter object for
  #               the corresponding regression function.  If any of
  #               these objects is a functional data object, it is
  #               converted to the default functional parameter object.
  #               if BETALIST is a functional parameter object
  #               it is converted to a list of length 1.
  #  WT       ... a vector of nonnegative weights for observations
  #  Y2CMAP   ... the matrix mapping from the vector of observed values
  #               to the coefficients for the dependent variable.
  #               This is output by function SMOOTH_BASIS.  If this is
  #               supplied, confidence limits are computed, otherwise not.
  #  SIGMAE   ... Estimate of the covariances among the residuals.  This
  #               can only be estimated after a preliminary analysis
  #               with FREGRESS.
  #
  #  Returns FREGRESSLIST  ... A list containing seven members with names:
  #    yfdPar      ... first  argument of FREGRESS
  #    xfdlist     ... second argument of FREGRESS
  #    betalist    ... third  argument of FREGRESS
  #    betaestlist ... estimated regression functions
  #    yhatfdobj   ... functional data object containing fitted functions
  #    Cmat        ... coefficient matrix for the linear system defining
  #                    the regression coefficient basis vector
  #    Dmat        ... right side vector for the linear system defining
  #                    the regression coefficient basis vector
  #    Cmatinv     ... inverse of the coefficient matrix, needed for
  #                    function FREGRESS.STDERR that computes standard errors
  #    wt          ... weights for observations
  #    df          ... degrees of freedom for fit
  
  #  Last modified 7 January 2020 by Jim Ramsay
  
  UseMethod("fRegress")
  
  if (inherits(y, "formula")) {
    fRegressList <- fRegress.formula(y, xfdlist, betalist, wt=NULL, y2cMap=NULL, 
                                     SigmaE=NULL, method, sep, ...)
    return(fRegressList)
  }
  
  if (inherits(y, "numeric")) {
    data <- xfdlist
    fRegressList <- fRegress.numeric(y, xfdlist, betalist, wt=NULL, y2cMap=NULL, 
                                     SigmaE=NULL, method, sep, ...)
    return(fRegressList)
    
  } else {
    
    arglist <- fRegressArgCheck(y, xfdlist, betalist, wt)
    
    yfdPar   <- arglist$yfdPar
    xfdlist  <- arglist$xfdlist
    betalist <- arglist$betalist
    wt       <- arglist$wt
    
    p <- length(xfdlist)
    N <- dim(xfdlist[[1]]$coef)[1]
    
    wtconst <- var(wt) == 0
    
    #  extract dependent variable information
    
    yfdobj    <- yfdPar$fd
    ylambda   <- yfdPar$lambda
    yLfdobj   <- yfdPar$Lfd
    ycoef     <- yfdobj$coefs
    ycoefdim  <- dim(ycoef)
    N         <- ycoefdim[2]
    ybasisobj <- yfdobj$basis
    rangeval  <- ybasisobj$rangeval
    ynbasis   <- ybasisobj$nbasis
    onesbasis <- create.constant.basis(rangeval)
    onesfd    <- fd(1,onesbasis)
    
    if (length(ycoefdim) > 2) stop("YFDOBJ from YFDPAR is not univariate.")
    
    #  --------  set up the linear equations for the solution  -----------
    
    #  compute the total number of coefficients to be estimated
    
    ncoef <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      if (betafdParj$estimate) {
        ncoefj     <- betafdParj$fd$basis$nbasis
        ncoef      <- ncoef + ncoefj
      }
    }
    
    Cmat <- matrix(0,ncoef,ncoef)
    Dmat <- rep(0,ncoef)
    
    #  loop through rows of CMAT
    
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      if (betafdParj$estimate) {
        betafdj    <- betafdParj$fd
        betabasisj <- betafdj$basis
        ncoefj     <- betabasisj$nbasis
        #  row indices of CMAT and DMAT to fill
        mj1    <- mj2 + 1
        mj2    <- mj2 + ncoefj
        indexj <- mj1:mj2
        #  compute right side of equation DMAT
        xfdj <- xfdlist[[j]]
        if (wtconst) {
          xyfdj <- xfdj*yfdobj
        } else {
          xyfdj <- (xfdj*wt)*yfdobj
        }
        wtfdj <- sum(xyfdj)
        Dmatj <- inprod(betabasisj,onesfd,0,0,rangeval,wtfdj)
        Dmat[indexj] <- Dmatj
        #  loop through columns of CMAT
        mk2 <- 0
        for (k in 1:j) {
          betafdPark <- betalist[[k]]
          if (betafdPark$estimate) {
            betafdk    <- betafdPark$fd
            betabasisk <- betafdk$basis
            ncoefk     <- betabasisk$nbasis
            #  column indices of CMAT to fill
            mk1 <- mk2 + 1
            mk2 <- mk2 + ncoefk
            indexk <- mk1:mk2
            #  set up two weight functions
            xfdk <- xfdlist[[k]]
            if (wtconst) {
              xxfdjk <- xfdj*xfdk
            } else {
              xxfdjk <- (xfdj*wt)*xfdk
            }
            wtfdjk <- sum(xxfdjk)
            Cmatjk <- inprod(betabasisj, betabasisk, 0, 0,
                             rangeval, wtfdjk)
            Cmat[indexj,indexk] <- Cmatjk
            Cmat[indexk,indexj] <- t(Cmatjk)
          }
        }
        #  attach penalty term to diagonal block
        lambdaj <- betafdParj$lambda
        if (lambdaj > 0) {
          Rmatj <- betafdParj$penmat
          if (is.null(Rmatj)) {
            Lfdj  <- betafdParj$Lfd
            Rmatj <- eval.penalty(betabasisj, Lfdj)
          }
          Cmat[indexj,indexj] <- Cmat[indexj,indexj] +
            lambdaj*Rmatj
        }
      }
    }
    
    Cmat <- (Cmat+t(Cmat))/2
    
    #  check Cmat for singularity
    
    eigchk(Cmat)
    
    
    #  solve for coefficients defining BETA
    
    Lmat    <- chol(Cmat)
    Lmatinv <- solve(Lmat)
    Cmatinv <- Lmatinv %*% t(Lmatinv)
    
    betacoef <- Cmatinv %*% Dmat
    
    #  set up fdPar objects for reg. fns. in BETAESTLIST
    
    betaestlist <- betalist
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      if (betafdParj$estimate) {
        betafdj <- betafdParj$fd
        ncoefj  <- betafdj$basis$nbasis
        mj1     <- mj2 + 1
        mj2     <- mj2 + ncoefj
        indexj  <- mj1:mj2
        coefj   <- betacoef[indexj]
        betafdj$coefs <- as.matrix(coefj)
        betafdParj$fd <- betafdj
      }
      betaestlist[[j]] <- betafdParj
    }
    
    #  set up fd objects for predicted values in YHATFDOBJ
    
    nfine     <- max(501,10*ynbasis+1)
    tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
    yhatmat <- matrix(0,nfine,N)
    for (j in 1:p) {
      xfdj       <- xfdlist[[j]]
      xmatj      <- eval.fd(tfine, xfdj, 0)
      betafdParj <- betaestlist[[j]]
      betafdj    <- betafdParj$fd
      betavecj   <- eval.fd(tfine, betafdj, 0)
      yhatmat    <- yhatmat + xmatj*as.vector(betavecj)
    }
    yhatfd <- smooth.basis(tfine, yhatmat, ybasisobj)$fd
    
    df <- NA
    
    #  -----------------------------------------------------------------------
    #        Compute pointwise standard errors of regression coefficients
    #               if both y2cMap and SigmaE are supplied.
    #  -----------------------------------------------------------------------
    
    if (!(is.null(y2cMap) || is.null(SigmaE))) {
      
      #  check dimensions of y2cMap and SigmaE
      
      y2cdim = dim(y2cMap)
      if (y2cdim[1] != ynbasis || y2cdim[2] != dim(SigmaE)[1]) {
        stop("Dimensions of Y2CMAP not correct.")
      }
      
      ybasismat = eval.basis(tfine, ybasisobj, 0)
      
      deltat    = tfine[2] - tfine[1]
      
      #  compute BASISPRODMAT
      
      basisprodmat = matrix(0,ncoef,ynbasis*N)
      
      mj2 = 0
      for (j in 1:p) {
        betafdParj = betalist[[j]]
        betabasisj = betafdParj$fd$basis
        ncoefj     = betabasisj$nbasis
        bbasismatj = eval.basis(tfine, betabasisj, 0)
        xfdj       = xfdlist[[j]]
        tempj      = eval.fd(tfine, xfdj, 0)
        #  row indices of BASISPRODMAT to fill
        mj1    = mj2 + 1
        mj2    = mj2 + ncoefj
        indexj = mj1:mj2
        #  inner products of beta basis and response basis
        #    weighted by covariate basis functions
        mk2 = 0
        for (k in 1:ynbasis) {
          #  row indices of BASISPRODMAT to fill
          mk1    = mk2 + 1
          mk2    = mk2 + N
          indexk = mk1:mk2
          tempk  = bbasismatj*ybasismat[,k]
          basisprodmat[indexj,indexk] =
            deltat*crossprod(tempk,tempj)
        }
      }
      
      #  compute variances of regression coefficient function values
      
      c2bMap    = solve(Cmat,basisprodmat)
      VarCoef   = y2cMap %*% SigmaE %*% t(y2cMap)
      CVariance = kronecker(VarCoef,diag(rep(1,N)))
      bvar      = c2bMap %*% CVariance %*% t(c2bMap)
      betastderrlist = vector("list", p)
      mj2 = 0
      for (j in 1:p) {
        betafdParj = betalist[[j]]
        betabasisj = betafdParj$fd$basis
        ncoefj     = betabasisj$nbasis
        mj1 	     = mj2 + 1
        mj2 	     = mj2 + ncoefj
        indexj     = mj1:mj2
        bbasismat  = eval.basis(tfine, betabasisj, 0)
        bvarj      = bvar[indexj,indexj]
        bstderrj   = sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
        bstderrfdj = smooth.basis(tfine, bstderrj, betabasisj)$fd
        betastderrlist[[j]] = bstderrfdj
      }
    } else {
      betastderrlist = NULL
      bvar           = NULL
      c2bMap         = NULL
    }
    
    
    #  -------------------------------------------------------------------
    #                       Set up output list object
    #  -------------------------------------------------------------------
    
    fRegressList <-
      list(yfdPar         = yfdPar,
           xfdlist        = xfdlist,
           betalist       = betalist,
           betaestlist    = betaestlist,
           yhat           = yhatfd,
           betastderrlist = betastderrlist,
           df             = df,
           y2cMap         = y2cMap,
           SigmaE         = SigmaE,
           bvar           = bvar,
           c2bMap         = c2bMap)
    
  }
  
}

#  ------------------------------------------------------------------------

eigchk <- function(Cmat) {
  
  #  check Cmat for singularity
  
  eigval <- eigen(Cmat)$values
  ncoef  <- length(eigval)
  if (eigval[ncoef] < 0) {
    neig <- min(length(eigval),10)
    cat("\nSmallest eigenvalues:\n")
    print(eigval[(ncoef-neig+1):ncoef])
    cat("\nLargest  eigenvalues:\n")
    print(eigval[1:neig])
    stop("Negative eigenvalue of coefficient matrix.")
  }
  if (eigval[ncoef] == 0) stop("Zero eigenvalue of coefficient matrix.")
  logcondition <- log10(eigval[1]) - log10(eigval[ncoef])
  if (logcondition > 12) {
    warning("Near singularity in coefficient matrix.")
    cat(paste("\nLog10 Eigenvalues range from\n",
              log10(eigval[ncoef])," to ",log10(eigval[1]),"\n"))
  }
}

#  ------------------------------------------------------------------------

fRegressArgCheck <- function(yfdPar, xfdlist, betalist, wt=NULL) 
{
  #  FREGRESS_ARGCHECK checks the first four arguments for the functions
  #  for function regression, including FREGRESS.
  
  #  --------------------  Check classes of arguments  --------------------
  
  #  check YFDPAR and compute sample size N
  
  if (inherits(yfdPar, "fd")) yfdPar <- fdPar(yfdPar)
  
  if (!(inherits(yfdPar, "fdPar") || is.numeric(yfdPar))) stop(
    "First argument is not of class 'fd', 'fdPar' or 'numeric'.")
  
  if (inherits(yfdPar, "fdPar")) {
    yfd   <- yfdPar$fd
    ycoef <- yfd$coefs
    N     <- dim(ycoef)[2]
  }
  
  if (is.numeric(yfdPar)) {
    N <- length(yfdPar)
  }
  
  #  check that xfdlist is a list object
  
  #  check XFDLIST
  
  if (inherits(xfdlist, "fd") || inherits(xfdlist, "numeric")) 
    xfdlist <- list(xfdlist)
  
  if (!inherits(xfdlist, "list")) stop(
    "Argument XFDLIST is not a list object.")
  
  #  get number of independent variables p
  
  p <- length(xfdlist)
  
  #  check BETALIST
  
  if (inherits(betalist, "fd")) betalist <- list(betalist)
  
  if (!inherits(betalist, "list")) stop(
    "Argument BETALIST is not a list object.")
  
  if (length(betalist) != p)  {
    cat(paste("\nNumber of regression coefficients does not match\n",
              "number of independent variables."))
    stop("")
  }
  
  #  check that the regression is functional, and extract the range
  
  if (inherits(yfdPar, "fdPar")) {
    rangeval <- yfdPar$fd$basis$rangeval
  } else {
    allscalar <- TRUE
    for (j in 1:p) {
      if (inherits(xfdlist[[j]], "fd")) {
        rangeval <- xfdlist[[j]]$basis$rangeval            
        allscalar <- FALSE
        break
      }
    }
    if (allscalar) stop(
      paste("The dependent variable and all the independent",   
            "variables are scalar."))
  }
  
  
  #  --------------------  check contents of arguments  -------------------
  
  #  XFDLIST:
  
  #  If the object is a vector of length N,
  #  it is converted to a functional data object with a
  #  constant basis
  
  onebasis <- create.constant.basis(rangeval)
  onesfd   <- fd(1,onebasis)
  
  xerror <- FALSE
  for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
      xcoef <- xfdj$coefs
      if (length(dim(xcoef)) > 2) stop(
        paste("Covariate",j,"is not univariate."))
      #  check size of coefficient array
      Nj <- dim(xcoef)[2]
      if (Nj != N) {
        print(
          paste("Incorrect number of replications in XFDLIST",
                "for covariate",j))
        xerror = TRUE
      }
    } 
    if (inherits(xfdj, "numeric")) {
      if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
      Zdimj <- dim(xfdj)
      if (Zdimj[1] != N && Zdimj != 1) {
        print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
        xerror = TRUE 
      } 
      if (Zdimj[2] != 1) {
        print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
        xerror = TRUE 
      } 
      xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
    } 
    if (!(inherits(xfdlist[[j]], "fd") || 
          inherits(xfdlist[[j]], "numeric"))) {
      print(paste("XFDLIST[[",j,"]] is neither an FD object nor numeric."))
      xerror = TRUE
    }
  }
  
  #  BETALIST:
  
  berror <- FALSE
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (inherits(betafdParj, "fd") || inherits(betafdParj, "basisfd")) {
      betafdParj    <- fdPar(betafdParj)
      betalist[[j]] <- betafdParj
    }
    if (!inherits(betafdParj, "fdPar")) {
      print(paste("BETALIST[[",j,"]] is not a FDPAR object."))
      berror <- TRUE
    }
  }
  
  if (xerror || berror) stop(
    "An error has been found in either XFDLIST or BETALIST.")
  
  #  check weights
  
  if (is.null(wt)) wt = rep(1,N)
  if (length(wt) != N) stop("Number of weights not equal to N.")
  if (any(wt < 0))     stop("Negative weights found.")
  
  return(list(yfdPar=yfdPar, xfdlist=xfdlist, betalist=betalist, wt=wt))
  
}