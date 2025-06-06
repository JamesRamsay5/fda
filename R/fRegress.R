
fRegress <- function(y, ...) {
  UseMethod("fRegress")
}

fRegress.fd <- function(y, xfdlist, betalist, wt=NULL,
                        y2cMap=NULL, SigmaE=NULL, returnMatrix=FALSE,
                        method=c('fRegress', 'model'),
                        sep='.', ...) {
  
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
  #                   a functional data object or a numerical vector
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
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  #
  #  Returns FREGRESSLIST  ... A list containing seven members with names:
  #    yfdobj      ... first  argument of FREGRESS
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
  #  This list object is converted to a class with the name "fRegress"
  #  function predict.fRegress is an example of a method that can be called simply
  #  as predict(fRegressList).  In this call fRegressList can be any object of the
  #  "fRegress".
  
  #  Last modified 5 November 2020 by Jim Ramsay
  
  if (is.fdPar(y)) y <- y$fd
  
  #  As of 2020, if yfd is an fdPar object, it is converted to an fd object.
  #  The added structure of the fdPar class is not used in any of the fRegress codes.
  #  The older versions of fda package used yfdPar as the name for the first member.
  
  arglist <- fRegressArgCheck(y, xfdlist, betalist, wt)
  
  yfdobj   <- arglist$yfd  # the older version used yfdPar as the name.
  xfdlist  <- arglist$xfdlist
  betalist <- arglist$betalist
  wt       <- arglist$wt
  
  p <- length(xfdlist)
  
  wtconst <- var(wt) == 0
  
  #  --------------------------------------------------------------------------
  #  branch depending on whether the dependent variable is functional or scalar
  #  --------------------------------------------------------------------------
  
  #  ----------------------------------------------------------------
  #           YFDOBJ is a functional data object
  #  ----------------------------------------------------------------
  
  #  extract dependent variable information
  ycoef     <- yfdobj$coefs
  ycoefdim  <- dim(ycoef)
  N         <- ycoefdim[2]
  ybasisobj <- yfdobj$basis
  rangeval  <- ybasisobj$rangeval
  ynbasis   <- ybasisobj$nbasis
  onesbasis <- create.constant.basis(rangeval)
  onesfd    <- fd(1,onesbasis)
  
  if (length(ycoefdim) > 2) stop("YFDOBJ from YFD is not univariate.")
  
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
  
  #  ------------------------------------------------------------------------
  #  Compute the symmetric positive definite matrix CMAT and
  #  the column matrix DMAT.  CMAT contains weighted inner products of
  #  bases for each pair of terms plus, for lambda > 0, a roughness penalty
  #  matrix to ensure that the estimated coefficient functions will be smooth
  #  The weight vector is the point-wise product of the associated functional
  #  covariates.  
  #  Dmat contains for each covariate the weighted integral of the basis 
  #  functions, with the weight function being the covariate function
  #  pointwise multiplied the dependent variate yobj.
  #  The estimated coefficients functions are defined by the solution
  #  CMAT %*% COEF = DMAT.
  #  ------------------------------------------------------------------------
  
  #  loop through rows of CMAT
  
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      #  get jth beta basis
      betafdj    <- betafdParj$fd
      betabasisj <- betafdj$basis
      ncoefj     <- betabasisj$nbasis
      #  row indices of CMAT and DMAT to fill
      mj1    <- mj2 + 1
      mj2    <- mj2 + ncoefj
      indexj <- mj1:mj2
      #  compute right side of equation DMAT
      #  compute weight function for DMAT
      xfdj <- xfdlist[[j]]
      if (wtconst) {
        xyfdj <- xfdj*yfdobj
      } else {
        xyfdj <- (xfdj*wt)*yfdobj
      }
      wtfdj <- sum(xyfdj)
      #  Compute jth component of DMAT
      Dmatj <- inprod(betabasisj,onesfd,0,0,rangeval,wtfdj)
      Dmat[indexj] <- Dmatj
      #  loop through columns of CMAT
      mk2 <- 0
      for (k in 1:j) {
        betafdPark <- betalist[[k]]
        if (betafdPark$estimate) {
          #  get the kth basis
          betafdk    <- betafdPark$fd
          betabasisk <- betafdk$basis
          ncoefk     <- betabasisk$nbasis
          #  column indices of CMAT to fill
          mk1 <- mk2 + 1
          mk2 <- mk2 + ncoefk
          indexk <- mk1:mk2
          #  set up weight function for CMAT component
          xfdk <- xfdlist[[k]]
          if (wtconst) {
            xxfdjk <- xfdj*xfdk
          } else {
            xxfdjk <- (xfdj*wt)*xfdk
          }
          wtfdjk <- sum(xxfdjk)
          #  compute the inner product
          Cmatjk <- inprod(betabasisj, betabasisk, 0, 0,
                           rangeval, wtfdjk)
          Cmat[indexj,indexk] <- Cmatjk
          Cmat[indexk,indexj] <- t(Cmatjk)
        }
      }
      #  attach penalty term to diagonal block if required
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
  
  #  ensure symmetry
  
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
  
  nfine   <- max(501,10*ynbasis+1)
  tfine   <- seq(rangeval[1], rangeval[2], len=nfine)
  yhatmat <- matrix(0,nfine,N)
  for (j in 1:p) {
    xfdj       <- xfdlist[[j]]
    xmatj      <- eval.fd(tfine, xfdj, 0, returnMatrix)
    betafdParj <- betaestlist[[j]]
    betafdj    <- betafdParj$fd
    betavecj   <- eval.fd(tfine, betafdj, 0, returnMatrix)
    yhatmat    <- yhatmat + xmatj*as.vector(betavecj)
  }
  yhatfdobj <- smooth.basis(tfine, yhatmat, ybasisobj)$fd
  
  df <- NA
  
  #  -----------------------------------------------------------------------
  #        Compute pointwise standard errors of regression coefficients
  #               if both y2cMap and SigmaE are supplied.
  #        y2cMap is supplied by the smoothing of the data that defined
  #        the dependent variable.
  #        SigmaE has to be computed from a previous analysis of the data.
  #  -----------------------------------------------------------------------
  
  if (!(is.null(y2cMap) || is.null(SigmaE))) {
    
    #  check dimensions of y2cMap and SigmaE
    
    y2cdim = dim(y2cMap)
    if (y2cdim[1] != ynbasis || y2cdim[2] != dim(SigmaE)[1]) {
      stop("Dimensions of Y2CMAP not correct.")
    }
    
    ybasismat = eval.basis(tfine, ybasisobj, 0, returnMatrix)
    
    deltat    = tfine[2] - tfine[1]
    
    #  compute BASISPRODMAT
    
    basisprodmat = matrix(0,ncoef,ynbasis*N)
    
    mj2 = 0
    for (j in 1:p) {
      betafdParj = betalist[[j]]
      betabasisj = betafdParj$fd$basis
      ncoefj     = betabasisj$nbasis
      bbasismatj = eval.basis(tfine, betabasisj, 0, returnMatrix)
      xfdj       = xfdlist[[j]]
      tempj      = eval.fd(tfine, xfdj, 0, returnMatrix)
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
      bbasismat  = eval.basis(tfine, betabasisj, 0, returnMatrix)
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
    list(yfdobj         = yfdobj,
         xfdlist        = xfdlist,
         betalist       = betalist,
         betaestlist    = betaestlist,
         yhatfdobj      = yhatfdobj,
         Cmat           = Cmat,
         Dmat           = Dmat,
         Cmatinv        = Cmatinv,
         wt             = wt,
         df             = df,
         y2cMap         = y2cMap,
         SigmaE         = SigmaE,
         betastderrlist = betastderrlist,
         bvar           = bvar,
         c2bMap         = c2bMap)
  
  
  return(fRegressList)
  
}

#  -------------------------------------------------------------------------------------

eigchk <- function(Cmat) {
  
  #  Last modified 25 August 2020 by Jim Ramsay
  
  #  Cmat for NA's
  
  if (anyNA(Cmat)) stop("Cmat has NA values.")
  
  #  check Cmat for Cmatmetry
  
  if (max(abs(Cmat-t(Cmat)))/max(abs(Cmat)) > 1e-10) {
    stop('CMAT is not symmetric.')
  } else {
    Cmat <- (Cmat + t(Cmat))/2
  }
  
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
