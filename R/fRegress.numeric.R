fRegress.numeric <- function(y, xfdlist, betalist, wt=NULL, y2cMap=NULL, 
                             SigmaE=NULL, method, sep='.', ...)
{
  
  #  FREGRESS.NUMERIC  Fits a functional linear model using multiple
  #  functional independent variables with the dependency being
  #  pointwise or concurrent.
  #  The case of a scalar independent variable is included by treating
  #  it as a functional independent variable with a constant basis
  #  and a unit coefficient.
  #
  #  Arguments:
  #  Y        ... an vectpr or matrix object for the dependent variable
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
  #                    function FREGRESSStERR that computes standard errors
  #    wt          ... weights for observations
  #    df          ... degrees of freedom for fit
  
  #  Last modified 8 January 2020 by Jim Ramsay
  
  if (!inherits(y,"numeric")) stop("Argument y is not numeric.")
  
  arglist <- fRegressArgCheck(y, xfdlist, betalist, wt)
  
  ymat     <- arglist$y
  xfdlist  <- arglist$xfdlist
  betalist <- arglist$betalist
  wt       <- arglist$wt
  
  p <- length(xfdlist)
  N <- dim(ymat)[1]
  
  wtconst <- var(wt) == 0
  Zmat  <- NULL
  Rmat  <- NULL
  pjvec <- rep(0,p)
  ncoef <- 0
  for (j in 1:p) {
    xfdj       <- xfdlist[[j]]
    xcoef      <- xfdj$coefs
    xbasis     <- xfdj$basis
    betafdParj <- betalist[[j]]
    bbasis     <- betafdParj$fd$basis
    bnbasis    <- bbasis$nbasis
    pjvec[j]   <- bnbasis
    Jpsithetaj <- inprod(xbasis,bbasis)
    Zmat       <- cbind(Zmat,crossprod(xcoef,Jpsithetaj))
    if (betafdParj$estimate) {
      lambdaj    <- betafdParj$lambda
      if (lambdaj > 0) {
        Lfdj  <- betafdParj$Lfd
        Rmatj <- lambdaj*eval.penalty(bbasis,Lfdj)
      } else {
        Rmatj <- matrix(0,bnbasis,bnbasis)
      }
      if (ncoef > 0) {
        zeromat <- matrix(0,ncoef,bnbasis)
        Rmat    <- rbind(cbind(Rmat,       zeromat),
                         cbind(t(zeromat), Rmatj))
      } else {
        Rmat  <- Rmatj
        ncoef <- ncoef + bnbasis
      }
    }
  }
  
  #  -----------------------------------------------------------
  #          set up the linear equations for the solution
  #  -----------------------------------------------------------
  
  #  solve for coefficients defining BETA
  
  if (any(wt != 1)) {
    rtwt   <- sqrt(wt)
    Zmatwt <- Zmat*rtwt
    ymatwt <- ymat*rtwt
    Cmat   <- t(Zmatwt) %*% Zmatwt + Rmat
    Dmat   <- t(Zmatwt) %*% ymatwt
  } else {
    Cmat <- t(Zmat) %*% Zmat + Rmat
    Dmat <- t(Zmat) %*% ymat
  }
  
  eigchk(Cmat)
  
  Cmatinv  <- solve(Cmat)
  
  betacoef <- Cmatinv %*% Dmat
  
  #  compute and print degrees of freedom measure
  
  df <- sum(diag(Zmat %*% Cmatinv %*% t(Zmat)))
  
  #  set up fdPar object for BETAESTFDPAR
  
  betaestlist <- betalist
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    betafdj    <- betafdParj$fd
    ncoefj     <- betafdj$basis$nbasis
    mj1        <- mj2 + 1
    mj2        <- mj2 + ncoefj
    indexj     <- mj1:mj2
    betacoefj        <- betacoef[indexj]
    betaestfdj       <- betafdj
    betaestfdj$coefs <- as.matrix(betacoefj)
    betaestfdParj    <- betafdParj
    betaestfdParj$fd <- betaestfdj
    betaestlist[[j]] <- betaestfdParj
  }
  
  #  set up fd object for predicted values
  
  yhatmat <- matrix(0,N,1)
  for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
      xbasis     <- xfdj$basis
      xnbasis    <- xbasis$nbasis
      xrng       <- xbasis$rangeval
      nfine      <- max(501,10*xnbasis+1)
      tfine      <- seq(xrng[1], xrng[2], len=nfine)
      deltat     <- tfine[2]-tfine[1]
      xmat       <- eval.fd(tfine, xfdj, 0)
      betafdParj <- betaestlist[[j]]
      betafdj    <- betafdParj$fd
      betamat    <- eval.fd(tfine, betafdj, 0)
      fitj       <- deltat*(crossprod(xmat,betamat) -
                              0.5*as.matrix((xmat[1,    ]*betamat[1,    ]) +
                                            (xmat[nfine,]*betamat[nfine,])))
      yhatmat    <- yhatmat + fitj
    } else{
      betaestfdParj <- betaestlist[[j]]
      betavecj      <- betaestfdParj$fd$coefs
      yhatmat       <- yhatmat + xfdj %*% t(betavecj)
    }
  }
  yhat <- yhatmat
  
  #  -----------------------------------------------------------------------
  #        Compute pointwise standard errors of regression coefficients
  #               if both y2cMap and SigmaE are supplied.
  #  -----------------------------------------------------------------------
  
  if (!(is.null(y2cMap) || is.null(SigmaE))) {
    
    #  check dimensions of y2cMap and SigmaE
    
    y2cdim <- dim(y2cMap)
    if (y2cdim[2] != dim(SigmaE)[1])  stop(
      "Dimensions of Y2CMAP not correct.")
    
    
    #  compute linear mapping c2bMap takinging coefficients for
    #  response into coefficients for regression functions
    
    c2bMap <- Cmatinv %*% t(Zmat)
    y2bmap <- c2bMap
    bvar   <- y2bmap %*% as.matrix(SigmaE) %*% t(y2bmap)
    betastderrlist <- vector("list",p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj     <- betabasisj$nbasis
      mj1        <- mj2 + 1
      mj2        <- mj2 + ncoefj
      indexj     <- mj1:mj2
      bvarj      <- bvar[indexj,indexj]
      betarng    <- betabasisj$rangeval
      nfine      <- max(c(501,10*ncoefj+1))
      tfine      <- seq(betarng[1], betarng[2], len=nfine)
      bbasismat  <- eval.basis(tfine, betabasisj, 0)
      bstderrj   <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
      bstderrfdj <- smooth.basis(tfine, bstderrj, betabasisj)$fd
      betastderrlist[[j]] <- bstderrfdj
    }
  } else {
    betastderrlist = NULL
    bvar           = NULL
    c2bMap         = NULL
  }
  
  #  -----------------------------------------------------------------------
  #                  Set up output list object
  #  -----------------------------------------------------------------------
  
  fRegressList <-
    list(yfdPar         = y,
         xfdlist        = xfdlist,
         betalist       = betalist,
         betaestlist    = betaestlist,
         yhat           = yhat,
         df             = df,
         y2cMap         = y2cMap,
         SigmaE         = SigmaE,
         betastderrlist = betastderrlist,
         bvar           = bvar,
         c2bMap         = c2bMap)

  return(fRegressList)

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

fRegressArgCheck <- function(y, xfdlist, betalist, wt=NULL) 
{
  #  FREGRESS_ARGCHECK checks the first four arguments 
  
  #  --------------------  Check classes of arguments  --------------------
  
  if (is.numeric(y)) {
    N <- length(y)
    y = as.matrix(y)
  } else {
    stop("First argument is not numeric.")
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
  
  rangeval = betalist[[1]]$fd$basis$rangeval
  
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
  
  return(list(y=y, xfdlist=xfdlist, betalist=betalist, wt=wt))
  
}
