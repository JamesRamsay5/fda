
#  Generator function of class fdPar

fdPar <- function(fdobj=NULL, Lfdobj=NULL, lambda=0, estimate=TRUE,
                  penmat=NULL){
  
  # Sets up a functional parameter object
  #  Arguments:
  #  FDOBJ    ... A functional data object.
  #               The basis for this object is used to define
  #               the functional parameter, or functional
  #               parameters of FDOBJ has replications.
  #               When an initial value is required for iterative
  #               estimation of a functional parameter, the coefficients
  #               will give the initial values for the iteration.
  #               Note: The direct conversion of a basisfd object
  #               to a fdPar object is no longer allowed.
  #  LFDOBJ   ... A linear differential operator value or a derivative
  #               value for penalizing the roughness of the object.
  #               By default, this is 0.
  #  LAMBDA   ... The penalty parameter controlling the smoothness of
  #               the estimated parameter.  By default this is 0.
  #  ESTIMATE ... If TRUE, the parameter is estimated; if FALSE, the
  #               parameter is held fixed at this value.
  #               By default, this is 1.
  #  PENMAT   ... The penalty matrix.
  #               In repeated calls to SMOOTH_BASIS, if this is
  #               saved, then the penalty does not need evaluating
  #               repeatedly.  Don't use, though, if LFDOBJ or LAMBDA
  #               are changed in the calculation.
  #
  #  An alternative argument list:
  #  The first argument can also be a basis object.  In this case, an
  #  FD object is set up with an empty coefficient matrix.
  #
  #  Return:
  #  FDPAROBJ ... A functional parameter object
  
  #  Last modified 21 January 2024 by Jim Ramsay
  
  #  ----------------------------------------------------------------------
  #         Check that fdobj is a functional data object and
  #         get the number of basis functions
  #  ----------------------------------------------------------------------
  
  if(!inherits(fdobj, 'fd')) {
      print("Converting a basis object to an fdPar object is not allowed.")
      print("Define a functional data object with coefficients instead.")
      stop("")
  } else {
    # the first object is an fd object, and we need nbasis later
    nbasis <- fdobj$basis$nbasis
  }
  
  #  ----------------------------------------------------------------------
  #                            Check parameters
  #  ----------------------------------------------------------------------
  
  #  check Lfdobj, the linear derivative operator
  
  {
    if (is.null(Lfdobj)) {
      #  Lfdobj is NULL, construct a default linear derivative operator according to
      #  according tothe basis system.
      if(fdobj$basis$type=='fourier') {
        #  a fourier basis
        rng <- fdobj$basis$rangeval
        Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)
      } else {
        #  a spline basis
        norder <- {
          # extract the order of the basis system for a spline basis, and
          # otherwise the order is 2 by default
          if (fdobj$basis$type=='bspline') norder.bspline(fdobj$basis)
          else 2
        }
        # the constructed differential operator 
        Lfdobj <- int2Lfd(max(0, norder-2))
      }
    }
    else
      #  A vector of differential orders times their coefficients
      #  is supplied.  convert to the linear differential operator.
      Lfdobj <- int2Lfd(Lfdobj)
  }
  
  if (!inherits(Lfdobj, "Lfd"))
    stop("'Lfdobj' is not a linear differential operator object.")
  
  #  check lambda
  
  if (!is.numeric(lambda)) stop("Class of LAMBDA is not numeric.")
  if (lambda < 0) stop("LAMBDA is negative.")
  
  #  check estimate
  
  if (!is.logical(estimate)) stop("Class of ESTIMATE is not logical.")
  
  #  check penmat if supplied
  
  if (!is.null(penmat)) {
    if (!is.numeric(penmat)) stop("PENMAT is not numeric.")
    #    penmatsize <- size(penmat)
    penmatsize <- dim(penmat)
    if (any(penmatsize != nbasis)) stop("Dimensions of PENMAT are not correct.")
  }
  
  #  ----------------------------------------------------------------------
  #                    set up the fdPar object
  #  ----------------------------------------------------------------------
  
  #  S3 definition
  
  fdParobj <- list(fd=fdobj, Lfd=Lfdobj, lambda=lambda, estimate=estimate,
                   penmat=penmat)
  
  oldClass(fdParobj) <- "fdPar"
  
  fdParobj
  
}

#  ----------------------------------------------------------------------

#  "print" method for "fdPar"

print.fdPar <- function(x, ...)
{
  object <- x
  cat("Functional parameter object:\n\n")
  print("Functional data object:")
  print.fd(object$fd)
  print("Linear differential operator object:")
  print.Lfd(object$Lfd)
  cat(paste("\nSmoothing parameter =",object$lambda,"\n"))
  cat(paste("\nEstimation status =",object$estimate,"\n"))
  if (!is.null(object$penmat)) {
    print("Penalty matrix:")
    print(object$penmat)
  }
}

#  ----------------------------------------------------------------------

#  "summary" method for "fdPar"

summary.fdPar <- function(object, ...)
{
  cat("Functional parameter object:\n\n")
  print("Functional data object:")
  summary.fd(object$fd)
  print("Linear differential operator object:")
  summary.Lfd(object$Lfd)
  cat(paste("\nSmoothing parameter =",object$lambda,"\n"))
  cat(paste("\nEstimation status =",object$estimate,"\n"))
  if (!is.null(object$penmat))
    print(paste("Penalty matrix dimensions:",dim(object$penmat)))
}


