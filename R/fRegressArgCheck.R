fRegressArgCheck <- function(yfd, xfdlist, betalist, wt=NULL) {
  #  FREGRESS_ARGCHECK checks the first four arguments for the functions
  #  for function regression, including FREGRESS.
  
  #  Last modified 16 December 2020 by Jim Ramsay
  
  #  --------------------  Check classes of arguments  --------------------
  
  #  check that YFD is of class either 'fd' or 'numeric' and compute sample size N
  
  if (!(is.fdPar(yfd) || is.fd(yfd) || is.numeric(yfd) || is.matrix(yfd))) stop(
    "First argument is not of class 'fdPar', 'fd', 'numeric' or 'matrix'.")
  
  #  As of 2020, if yfd is an fdPar object, it is converted to an fd object.
  #  The added structure of the fdPar class is not used in any of the fRegress codes.
  # The older versions of fda package used yfdPar as the name for the first member.
  
  if (is.fdPar(yfd)) yfd <- yfd$fd
  
  if (inherits(yfd, "fd")) {
    ycoef <- yfd$coefs
    N     <- dim(ycoef)[2]
  } else {
    N <- length(yfd)
  } 
  
  #  check that xfdlist is a list object and compute number of covariates p
  
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
    
  #  extract the range if YFD is functional
  
  if (inherits(yfd, "fd")) {
    rangeval <- yfd$basis$rangeval
  } else {
    rangeval = c(0,1)
  #   allscalar <- TRUE
  #   for (j in 1:p) {
  #     if (inherits(xfdlist[[j]], "fd")) {
  #       rangeval <- xfdlist[[j]]$basis$rangeval            
  #       allscalar <- FALSE
  #       break
  #     }
  #   }
    # if (allscalar) stop(
    #   paste("The dependent variable and all the independent",   
    #         "variables are scalar."))
  }
  
 
  
  #  --------------------  check contents of BETALIST  -------------------
  
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
  
  
  
  
 #  --------------------  check contents of XFDLIST  -------------------
  
  #  If the object is a vector of length N,
  #  it is converted to a functional data object with a
  #  constant basis
  
 
  
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
      xfdlist[[j]] <- fd(matrix(xfdj,1,N), create.constant.basis(betalist[[j]]$fd$basis$rangeval))
    } 
    if (!(inherits(xfdlist[[j]], "fd"     ) || 
          inherits(xfdlist[[j]], "numeric") ||
          inherits(xfdlist[[j]], "matrix" ))) {
      print(paste("XFDLIST[[",j,"]] is not an FD or numeric or matrix object."))
      xerror = TRUE
    }
  }  
  
  
  
  
  
  
  
  
  if (xerror || berror) stop(
    "An error has been found in either XFDLIST or BETALIST.")
  
  #  --------------------  check contents of WEIGHTS  -------------------
  
  if (is.null(wt)) wt = rep(1,N)
  if (length(wt) != N) stop("Number of weights not equal to N.")
  if (any(wt < 0))     stop("Negative weights found.")
  
  #  ---------------------  return the argument list  --------------------
  
  # The older versions of fda package used yfdPar as the name for the first member.
  
  return(list(yfd=yfd, xfdlist=xfdlist, betalist=betalist, wt=wt))
  
}