Data2fd <- function(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
                    lambda=3e-8/diff(as.numeric(range(argvals))),
                    fdnames=NULL, covariates=NULL, method="chol") {
  #  Arguments:
  # ARGVALS  A set of N argument values, set by default to equally spaced
  #             on the unit interval (0,1).
  # Y        an array containing values of curves
  #             If the array is a matrix, rows must correspond to argument
  #             values and columns to replications, and it will be assumed
  #             that there is only one variable per observation.
  #             If Y is a three-dimensional array, the first dimension
  #             corresponds to argument values, the second to replications,
  #             and the third to variables within replications.
  #             If Y is a vector, only one replicate and variable are assumed.
  # BASISROBJ A functional basis object.  
  # NDERIV    An order of derivative to be used if a roughness penalty is
  #           set up.
  # LAMBDA.   A positive scalar that multplies the roughness penalty to 
  #           control the amount of smoothness.
  # FDNAMES  A cell of length 3 with names for
  #             1. argument domain, such as "Time"
  #             2. replications or cases
  #             3. the function.
  # COVARIATES  A N by Q matrix Z of covariate values used to augment
  #             the smoothing function, where N is the number of
  #             data values to be smoothed and Q is the number of
  #             covariates.  The process of augmenting a smoothing
  #             function in this way is often called "semi-parametric
  #             regression".  The default is the null object NULL.
  # METHOD      The method for computing coefficients.  The usual method
  #             computes cross-product matrices of the basis value matrix,
  #             adds the roughness penalty, and uses the Choleski
  #             decomposition of this to compute coefficients, analogous
  #             to using the normal equations in least squares fitting.
  #             But this approach, while fast, contributes unnecessary
  #             rounding error, and the qr decomposition of the augmented
  #             basis matrix is prefererable.  But nothing comes for free,
  #             and the computational overhead of the qr approach can be a
  #             serious problem for large problems (n of 1000 or more).
  #             For this reason, the default is "method" = "chol", but if

  #argvals=NULL 
  #y=NULL 
  #basisobj=NULL 
  #nderiv=NULL
  #lambda=3e-8/diff(as.numeric(range(argvals)))
  #fdnames=NULL 
  #covariates=NULL 
  #method="chol"
  
  #. Last modified 24 January 2024 by Jim Ramsay
  #. Oiginally designed by Spencer Graves
  
  #  call function argvalsySwap() to swap arguments of needed
  
  argChk <- argvalsySwap(argvals, y, basisobj)
  
  if(!is.numeric(AV <- argChk$argvals)){
    if(is.null(AV))
      stop('is.null(argChk$argvals); should be numeric')
    cat('argChk$argvals is not numeric.\n')
    cat('class(argChk$argvals) = ', class(AV), '\n')
    print(AV)
  }  
  
  #. arguments now set up for smoothing, call function smooth.basisPar
  
  # print("invoking smooth.basisPar")
  smBasis <- smooth.basisPar(argChk$argvals, argChk$y,
                             fdobj=basisobj, Lfdobj=nderiv, lambda=lambda,
                             fdnames=fdnames,
                             covariates=covariates, method="chol")
  # smooth.basis <- function(argvals, y, basisobj,
  #                          wtvec=NULL,   fdnames=NULL, covariates=NULL,
  #                          method="chol", dfscale=1, returnMatrix=FALSE) {
    # smBasis <- smooth.basisPar(argChk$argvals, argChk$y,
  #                            fdobj=basisobj, Lfdobj=nderiv, lambda=lambda,
  #                            fdnames=fdnames,
  #                            covariates=covariates, method="chol")
  
  # return functional data object extracted from list object smBasis
  
  return(smBasis$fd)
}

#  -------------------------------------------------------------------------

## 2020-01-16:  Spencer Graves makes argvalsySwap 
## an internal function
argvalsySwap = function(argvals=NULL, y=NULL, basisobj=NULL)
{
  # print("1.  if y is NULL")
  ##
  ## 1.  if y is NULL, then (1) if argvals is also missing, stop
  ##                   else (2) redefine argvals as y
  ##
  if(is.null(y)){
    if(is.null(argvals)) 
      #. both argvals and yj are NUll, terminate the analysis
      stop("'y' is missing with no default")
      #   argvals are present, store argvals as y
      cat("'y' is missing, using 'argvals'\n") 
      #  swap argvals and y
      y       <- argvals
      argvals <- NULL 
    }
  # print("2.  carry on")
  ##
  ## 2.  carry on, 
  ##    constructing substitute for missing argument objects if needed
  ##
  dimy <- dim(as.array(y))
  if(is.null(argvals)) { 
    {
    #  argvals is missing, construct it from dimensions of y
    if(is.null(basisobj)) {
      #. basis object is also missing, construct its default object of order 1
      #. see function basisfd() for this.
      basisobj <- create.bspline.basis(basisobj)
    } else {
      if(is.numeric(basisobj)) {
        if(any(basisobj < 0)) {
          stop("Basis object cannot have negative order.")
        } else {
          #  basis object is present but is a numeric value
          if(length(basisobj) > 1) {
            # is a numeric vector of length > 1, 
            # use default basis object with order 1
            basisobj <- create.bspline.basis(basisobj)
          } else {
            # single non-negative value, 
            # make default basis object of order value
            basisobj <- create.bspline.basis(norder=basisobj)
          }
        }
      } else {
        # print("if(inherits(basisobj, 'fd'))")
        # print(inherits(basisobj, 'fd'))
        if(inherits(basisobj, 'fd')) {
          # basis object is actually a function data object,
          # use its basis object
            basisobj <- basisobj$basis
        } else 
          # print("if(inherits(basisobj, 'fdPar'))")
          # print(inherits(basisobj, 'fdPar'))
          if(inherits(basisobj, 'fdPar')) {
            basisobj <- basisobj$fd$basis
          #} else {
            #stop("A basis object cannot be constructed from object provided.")
          #}
        }
      }
    }
    }
    #  argument is a basis object, carry on to set its range
    # print("carry on to set its range")
    a01 <- basisobj$rangeval
    if(is.null(a01)) 
      stop('basisobj does not have a required ',
           'rangeval component.')
    #} else {
      #. rangeval if present, carry on to construct argval if needed
      n <- dimy[1]
      cat(paste("'argvals' is missing;  using seq(", a01[1],
              ", ", a01[2], ", length=", n, ")\n"))  
      # argvals vector is equally spaced between rangeval values
      argvals <- seq(a01[1], a01[2], length=n)
    #}
    #  missing objects now complete, return to smooth.fdPar
      # print("return")
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  # print("3. consistency of argvals and y")
  ##
  ## 3. Arguments argvals, y and basisobj are in all correct as supplied
  ##
  #    Now check consistency of argvals and y
  dima <- dim(as.array(argvals)) 
  { 
    if(length(dimy) == length(dima)) {
      if(any(dimy != dima)) 
        stop("dimensions of 'argvals' and 'y' must be compatible;\n",
             "  dim(argvals) = ", paste(dima, collapse=' x '),
             ";  dim(y) = ", paste(dimy, collapse=' x ') )
        #   argvals and y are consistent now check basisobj
        if(inherits(basisobj, 'fd')) {
          basisobj <- basisobj$basis
        } else {
          if(inherits(basisobj, 'fdPar')) {
            basisobj <- basisobj$fd$basis
          } else {
            #  basisobj is neither an fd or an fdPar object
            #. if an integer, matrix or array define a suitable basis object,
            #. otherwise stop with an error message
            if(inherits(basisobj, 'array')) {
              #. basisobj is a matrix or array
              basisobj <- fd.$basis
            } else { 
              if(inherits(basisobj, 'integer')) {
                basisobj <- create.bspline.basis(argvals, norder=basisobj)
              } else {
                if(is.null(basisobj)) {
                  basisobj <- create.bspline.basis(argvals)
                } else {
                  if(!inherits(basisobj, 'basisfd'))
                    stop("'basisobj' is NOT a functional basis",
                         " object (class 'basisfd');  class = ",
                         class(basisobj)[1])
                }
              }
            }
          }
        }
    }
    #. extract rangeval vector of length 2
    a01  <- basisobj$rangeval
    #. compute range of argvals
    arng <- range(argvals)
    #. check range ofa01 and a02 as consistent argvals
    if ((a01[1]<=arng[1]) && (arng[2]<=a01[2])) {
      return(list(argvals=argvals, y=y, basisobj=basisobj))
    }
    yrng <- range(y)
    if((a01[1]<=yrng[1]) && (yrng[2]<=a01[2])) {
      cat(paste("'argvals' is NOT contained in basisobj$rangeval",
                ", but 'y' is;  swapping 'argvals' and 'y'.\n"))
      return(list(argvals=y, y=argvals, basisobj=basisobj)) 
    }
    stop("Neither 'argvals' nor 'y' are contained in ",
         "basisobj$rangeval")
  }
  # print("4. swap")
  ##
  ## 4.  If(length(dimy) < length(dima)) swap ...
  ##
  if(length(dimy)<length(dima)) {
    cat(paste("Swapping 'y' and 'argvals', because 'y' is ",
              "simpler,\n  and 'argvals' should be;  now ",
              "dim(argvals) = ", paste(dimy, collapse=" x "),
              ";  dim(y) = ", paste(dima, collapse=" x "),"\n" )) 
    y. <- argvals
    argvals <- y
    y <- y.
    #
    d. <- dima
    dima <- dimy
    dimy <- d.
  }   
  #
  if(any(dima != dimy[1:length(dima)]))
    stop("A dimension of 'argvals' does not match 'y':\n",
         "  dim(argvals) = ", paste(dima, collapse=" x "),
         ";  dim(y) = ", paste(dimy, collapse=" x ") )      
  # print("5. swap")
  ##        
  ## 5.  Check compatibility of argvals with basisobj
  ##        
  {
    if(inherits(basisobj, 'fd')) basisobj <- basisobj$basis
    else {
      if(inherits(basisobj, 'fdPar'))
        basisobj <- basisobj$fd$basis
      else {
        if(inherits(basisobj, 'array')){
          fd. <- fd(basisobj)
          basisobj <- fd.$basis
        }
        else { 
          if(inherits(basisobj, 'integer'))
            basisobj <- create.bspline.basis(argvals, norder=basisobj)
          else {
            if(is.null(basisobj))
              basisobj <- create.bspline.basis(argvals)
            else
              if(!inherits(basisobj, 'basisfd'))
                stop("'basisobj' is NOT a functional basis",
                     " object (class 'basisfd');  class = ",
                     class(basisobj)[1])
          }
        }
      }
    }
  }
  a01 <- basisobj$rangeval
  arng <- range(argvals)
  if((a01[1]<=arng[1]) && (arng[2]<=a01[2])) {
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  #
  stop("'argvals' are not contained in basisobj$rangeval")
}

