landmarkreg <- function(unregfd, ximarks, x0marks=xmeanmarks,
                        WfdPar=NULL, ylambda=1e-10) {
  #  Arguments:
  #  UNREGFD ... functional data object for curves to be registered
  #  XIMARKS ... N by NL array of times of interior landmarks for
  #                 each observed curve
  #  XOMARKS ... vector of length NL of times of interior landmarks for
  #                 target curve
  #  WFDPAR  ... a functional parameter object defining a warping function
  #  YLAMBDA ... smoothing parameter to be used in computing the registered
  #                 functions.  For high dimensional bases, local wiggles may be
  #                 found in the registered functions or its derivatives that are
  #                 not seen in the unregistered functions.  In this event, this
  #                 parameter should be increased.
  #  Returns:
  #  FDREG   ... a functional data object for the registered curves
  #  WARPFD  ... a functional data object for the warping functions
  #  WFD     ... a functional data object for the W functions defining the
  #              warping functions
  
  # Warning:  As of March 2022, landmark registration cannot be done using
  # function smooth.basis instead of function smooth.morph.  The 
  # warping function must be strictly monotonic, and we have found that using 
  # smooth.basis too often violates this contraint.  Function 
  # smooth.morph ensures monotonicity.
  
  #  Last modified 29 March 2022  by Jim Ramsay
  
  #  check unregfd
  
  if (!(inherits(unregfd,  "fd"))) stop(
    "Argument unregfd  not a functional data object.")
  
  basisobj <- unregfd$basis
  nbasis   <- basisobj$nbasis
  rangeval <- basisobj$rangeval
  
  #   ---------------------------------------------------------------------
  #                  check ximarks and x0marks
  #   ---------------------------------------------------------------------
  
  #  check ximarks being matrix with ncurve rows and nmarks columns
  
  if (is.numeric(ximarks)) {
    nximarks <- length(ximarks)
    # if ximarks is a vector, coerce it to a single row matrix
    if (is.vector(ximarks))     ximarks <- matrix(ximarks,1,nximarks)
    # if ximarks is a data.frame, coerce it to a matrix
    if (is.data.frame(ximarks)) ximarks <- as.matrix(ximarks)
  } else {
    stop("Argument ximarks is not numeric.")
  }
  
  #  compute row-wise mean of ximarks to serve as x0marks if
  #  needed
  
  xmeanmarks <- apply(ximarks,2,mean)
  
  #  check x0marks and coerce it to be a one-row matrix
  
  if (is.numeric(x0marks)) {
    nx0marks <- length(x0marks)
    if (is.vector(x0marks)) x0marks <- matrix(x0marks,1,nx0marks)
  } else {
    stop("Argument x0marks is not numeric.")
  }
  
  #  check that ximarks and x0marks have same number of columns
  
  if (ncol(ximarks) != length(x0marks)) 
    stop("The number of columns in ximarks is not equal to length of x0marks.")
  
  # check that ximarks are within range of unregfd
  
  if (any(ximarks <= rangeval[1]) || any(ximarks >= rangeval[2]))
    stop("Argument ximarks has values outside of range of unregfd.")
  
  # check that x0marks are within range of unregfd
  
  if (any(x0marks <= rangeval[1]) || any(x0marks >= rangeval[2]))
    stop("Argument x0marks has values outside of range of unregfd.")
  
  #  determine the number of curves to be registered
  
  ncurve   <- dim(ximarks)[1]
  
  #   ---------------------------------------------------------------------
  #                        check WFDPAR
  #   ---------------------------------------------------------------------
  
  #  set up default WfdPar
  
  if (is.null(WfdPar)) {
    Wnbasis   <- length(x0marks) + 2
    Wbasis    <- create.bspline.basis(rangeval, Wnbasis)
    Wfd       <- fd(matrix(0,Wnbasis,ncurve), Wbasis)
    WfdPar    <- fdPar(Wfd, 2, 1e-10)
  }
  
  WfdPar <- fdParcheck(WfdPar, ncurve)
  
  #  set up WFD0 and WBASIS
  
  Wfd0    <- WfdPar$fd
  WLfd    <- WfdPar$Lfd
  Wbasis  <- Wfd0$basis
  Wrange  <- Wbasis$rangeval
  Wlambda <- WfdPar$lambda
  
  #   ---------------------------------------------------------------------
  #                        set up analysis
  #   ---------------------------------------------------------------------
  
  nfine   <- min(c(101,10*nbasis))
  xfine   <- seq(rangeval[1],rangeval[2],length=nfine)
  yfine   <- eval.fd(xfine, unregfd)
  yregmat <- yfine
  hfunmat <- matrix(0,nfine,ncurve)
  hinvmat <- matrix(0,nfine,ncurve)
  Wlambda <- max(Wlambda,1e-10)
  
  xval    <- c(rangeval[1],x0marks,rangeval[2])
  Wnbasis <- Wbasis$nbasis
  Wcoef   <- matrix(0,Wnbasis,ncurve)
  nval    <- length(xval)
  
  #  --------------------------------------------------------------------
  #                  Iterate through curves to register
  #  --------------------------------------------------------------------
  
  if (ncurve > 1) cat("Progress:  Each dot is a curve\n")
  
  for (icurve in 1:ncurve) {
    if (ncurve > 1) cat(".")
    #  set up landmark times for this curve
    yval   <- c(rangeval[1],ximarks[icurve,],rangeval[2])
    #  smooth relation between this curve"s values and target"s values
    #  use monotone smoother
    Wfd  <- smooth.morph(xval, yval, Wrange, WfdPar)$Wfdobj
    hfun <- monfn(xfine, Wfd)
    b    <- (rangeval[2]-rangeval[1])/(hfun[nfine]-hfun[1])
    a    <- rangeval[1] - b*hfun[1]
    hfun <- a + b*hfun
    hfun[c(1,nfine)] <- rangeval
    Wcoefi           <- Wfd$coef
    Wcoef[,icurve]   <- Wcoefi
    hfunmat[,icurve] <- hfun
    
    #  compute h-inverse  in order to register curves
    
    Wcoefi       <- Wfd$coefs
    Wfdinv       <- fd(-Wcoefi,Wbasis)
    WfdParinv    <- fdPar(Wfdinv, WLfd, Wlambda)
    Wfdinv       <- smooth.morph(hfun, xfine, Wrange, WfdParinv)$Wfdobj
    hinv         <- monfn(xfine, Wfdinv)
    b            <- (rangeval[2]-rangeval[1])/(hinv[nfine]-hinv[1])
    a            <- rangeval[1] - b*hinv[1]
    hinv         <- a + b*hinv
    hinv[c(1,nfine)] <- rangeval
    
    #  compute registered curves
    
    yregfd <- smooth.basis(hinv, yfine[,icurve], basisobj)$fd
    yregmat[,icurve] <- eval.fd(xfine, yregfd, 0)
  }
  
  if (ncurve > 1) cat("\n")
  
  #  create functional data objects for the registered curves
  
  regfdPar <- fdPar(basisobj, 2, ylambda)
  regfd    <- smooth.basis(xfine, yregmat, regfdPar)$fd
  regnames <- unregfd$fdnames
  names(regnames)[3] <- paste("Registered",names(regnames)[3])
  regfd$fdnames <- regnames
  
  #  create functional data objects for the warping functions
  
  warpfd                <- smooth.basis(xfine, hfunmat, basisobj)$fd
  warpfdnames           <- unregfd$fdnames
  names(warpfdnames)[3] <- paste("Warped",names(regnames)[1])
  warpfd$fdnames        <- warpfdnames
  
  #  create functional data objects for the inverse warping functions
  
  warpfdinv             <- smooth.basis(xfine, hinvmat, basisobj)$fd
  warpfdnames           <- unregfd$fdnames
  names(warpfdnames)[3] <- paste("Warped",names(regnames)[1])
  warpfdinv$fdnames     <- warpfdnames
  
  #  The core function defining the strictly monotone warping
  
  Wfd <- fd(Wcoef, Wbasis)
  
  return( list(regfd=regfd, warpfd=warpfd, warpfdinv=warpfdinv, Wfd = Wfd) )
}
