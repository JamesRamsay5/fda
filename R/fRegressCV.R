fRegressCV <- function(y, xfdlist, betalist, wt=NULL, CVobs=1:N, ...)
{
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
  
#  Note: In earlier version named fRegress.CV an argument CVobs was 
#  included.  Since now S3 methods for a particular object must all
#  share the same argument list, CVobs has been replaced by 1:N.
#  If that choice is not appropriate, a subset of the data can be 
#  defined and re-analyzed to achieve the same effect.

# last modified 11 January 2020 by Jim Ramsay

# extract dimensions of the data and analysis

p <- length(xfdlist)
N <- dim(xfdlist[[1]]$coef)[2]

#  branch to either scalar or functional dependent variable

if (inherits(y, "numeric"))  {

    #  scalar dependent variable case

    yvec   <- y
    SSECV  <- 0
    errfd  <- c()
    for (m in 1:length(CVobs)) {
      i        <- CVobs[m]
      #  eliminate case i from the weights
      wti <- wt[-i]
      xfdlisti <- vector("list",p)
      for (j in 1:p) {
        xfdj          <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
          betafdParj <- betalist[[j]]
          betafdj    <- betafdParj$fd
          basisj     <- betafdj$basis
          betarangej <- basisj$rangeval
          conbasisj  <- create.constant.basis(betarangej)
          xfdj       <- fd(matrix(xfdj,1,N), conbasisj)
        }
        basisj <- xfdj$basis
        coefj  <- xfdj$coefs
        if (dim(coefj)[1] == 1) coefj <- matrix(coefj[-i],1,N-1)
        else                    coefj <- as.matrix(coefj[,-i])
        xfdlisti[[j]] <- fd(coefj,basisj)
      }
      yveci         <- yvec[-i]
      fRegressListi <- fRegress(yveci, xfdlisti, betalist, wti)
      betaestlisti  <- fRegressListi$betaestlist
      yhati <- 0
      for (j in 1:p) {
        betafdParj <- betaestlisti[[j]]
        betafdj    <- betafdParj$fd
        xfdj       <- xfdlist[[j]]
        bbasisj    <- betafdj$basis
        rangej     <- bbasisj$rangeval
        nfine      <- max(501, bbasisj$nbasis*10+1)
        tfine      <- seq(rangej[1], rangej[2], len=nfine)
        delta      <- tfine[2]-tfine[1]
        betavec    <- eval.fd(tfine, betafdj, 0)
        xveci      <- eval.fd(tfine, xfdj[i], 0)
        yhati      <- yhati + delta*(sum(xveci*betavec) -
                                    0.5*( xveci[1]    *betavec[1] +
                                          xveci[nfine]*betavec[nfine] ))
      }
      errfd[i] <- yvec[i] - yhati;
      SSECV <- SSECV + errfd[i]^2
    }
 } else {

    #  functional dependent variable case

    if (inherits(y, "fdPar")) yfd <- y$fd
    if (inherits(y, "fd"))    yfd <- y
    if (!inherits(y, "fdPar") && !inherits(y, "fd")) 
      stop("Argument y is invalid.")
    SSECV   <- 0
    errcoefs <- c()
    for(i in 1:N){
      # index of case to eliminate
      # eliminate case i from the weights
      wti <- wt[-i]
      # eliminate case i from covariates
      txfdlist <- xfdlist
      for(k in 1:p){
        txfdlist[[k]] <- xfdlist[[k]][-i]
      }
      # eliminate case i from dependent variable
      yfdi <- yfd[-i]
      # carry out the functional regression analysis
      tres <- fRegress(yfdi,txfdlist,betalist,wti)
      #  extract the regression coefficient functions
      betaestlisti <- tres$betaestlist
      #  compute the fit to the data for case i
      yhatfdi <- 0
      for(k in 1:p){
        betafdPark <- betaestlisti[[k]]
        betafdk    <- betafdPark$fd
        xfdk       <- xfdlist[[k]]
        xfdik      <- xfdk[i]
        tempfd     <- xfdik*betafdk
        yhatfdi    <- yhatfdi + tempfd
      }
      #  compute the residual function
      errfdi   <- yfd[i] - yhatfdi
      #  increment the error sum of squares by the integral of the
      #  square of the residual function
      SSECV   <- SSECV + inprod(errfdi,errfdi)
      #  add the coefficients for the residual function
      errcoefs <- cbind(errcoefs,errfdi$coefs)
    }
    #  set up the functional data object for the residual fns
    errfd <- fd(errcoefs,errfdi$basis)
    names(errfd$fdnames)[[3]] <- "Xval Errors"
}
return(list(SSECV=SSECV, errfdcv=errfd))
}


