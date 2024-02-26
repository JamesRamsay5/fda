Data2fd <- function(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
                    lambda=3e-8/diff(as.numeric(range(argvals))),
                    fdnames=NULL, covariates=NULL, method="chol") {
  
  #  Last modified 5 February 2024
  
  # Changes proposed by Spencer Graves 2010.12.08 ...
  # if(is.null(lambda))
  #   lambda <- 1e-9*sd(argChk$y)/diff(range(argChk$argvals))
  #
  # Error in smooth.basis ... argvals is not numeric
  # in R CMD check, cannot replicate line by line.  
  
  #. Tests six situations requiring modification or termination.
  
  argChk <- argvalsySwap(argvals, y, basisobj)
  
  # another rest with messages
  
  if(!is.numeric(AV <- argChk$argvals)){
    # terminal message
    if(is.null(AV)) stop('is.null(argChk$argvals); should be numeric')
    #. otherwise alert message
    cat('argChk$argvals is not numeric.\n')
    cat('class(argChk$argvals) = ', class(AV), '\n')
    print(AV)
  }  
  
  #. Success and smoothing ... S3 object of class fdSmooth is returned
  
  fdSmoothobj <- smooth.basisPar(argChk$argvals, argChk$y,
                                 fdobj=basisobj, Lfdobj=nderiv, lambda=lambda,
                                 fdnames=fdnames,
                                 covariates=covariates, method="chol")
  
}

#  -------------------------------------------------------------------------

## 2020-01-16:  Spencer Graves makes argvalsySwap 
## An internal function that tests for 6 situations that require modification
## with a warning, or a terminal error message

argvalsySwap = function(argvals=NULL, y=NULL, basisobj=NULL) {
  
  if (inherits(basisobj, 'basisfd')) rangeval <- basisobj$rangeval
  
  ##. --------------------------------------------------------------------------
  ## 1.  if(is.null(y)) use argvals
  ##. --------------------------------------------------------------------------
  
  if(is.null(y)){
    if(is.null(argvals)) stop("'y' is missing with no default")
    #   Store argvals as y and alert
    cat("'y' is missing, using 'argvals'\n") 
    y <- argvals
    argvals <- NULL 
  }
  
  ##. --------------------------------------------------------------------------
  ## 2.  if(is.null(argvals)). argvals <- seq(basisobj$rangeval, dim(y)[1])
  ##. --------------------------------------------------------------------------

  dimy <- dim(as.array(y))
  if(is.null(argvals)){
    {
      if(is.null(basisobj)){
        basisobj <- create.bspline.basis(basisobj)
      } else {
        if(is.numeric(basisobj)) {
          if(length(basisobj)>1){
            basisobj <- create.bspline.basis(basisobj)
          } else 
            basisobj <- create.bspline.basis(norder=basisobj)
        }
        else {
          if(inherits(basisobj, 'fd')){
            basisobj <- basisobj$basis
          } else 
            if(inherits(basisobj, 'fdPar'))
              basisobj <- basisobj$fd$basis
        }
      }
    }
    if(is.null(rangeval))
      stop('basisobj does not have a required ',
           'rangeval component.')
    #    
    n <- dimy[1]
    #. alert message
    cat(paste("'argvals' is missing;  using seq(", rangeval[1],
              ", ", rangeval[2], ", length=", n, ")\n"))       
    argvals <- seq(rangeval[1], rangeval[2], length=n)
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  
  ##. --------------------------------------------------------------------------
  ## 3.  if(length(dim(argvals)) == length(dim(y))) ... 
  ##. --------------------------------------------------------------------------

  dima <- dim(as.array(argvals))
  {
    if(length(dimy) == length(dima)){
      if(any(dimy != dima))
        #. terminal message
        stop("dimensions of 'argvals' and 'y' must be compatible;\n",
             "  dim(argvals) = ", paste(dima, collapse=' x '),
             ";  dim(y) = ", paste(dimy, collapse=' x ') )
      #     Check basisobj
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
                    #. terminal message
                    stop("'basisobj' is NOT a functional basis",
                         " object (class 'basisfd');  class = ",
                         class(basisobj)[1])
              }
            }
          }
        }
      }
      arng     <- range(argvals)
      if ((rangeval[1]<=arng[1]) && (arng[2]<=rangeval[2])) {
        return(list(argvals=argvals, y=y, basisobj=basisobj))
      }
      #
      yrng <- range(y)
      if((rangeval[1]<=yrng[1]) && (yrng[2]<=rangeval[2])) {
        #. alert message
        cat(paste("'argvals' is NOT contained in basisobj$rangeval",
                  ", but 'y' is;  swapping 'argvals' and 'y'.\n"))
        return(list(argvals=y, y=argvals, basisobj=basisobj)) 
      }
      #   Terminal message 
      stop("Neither 'argvals' nor 'y' are contained in ",
           "basisobj$rangeval")
    }
  }
  
  ##. --------------------------------------------------------------------------
  ## 4.  If(length(dimy) < length(dima)) swap argvals and y
  ##. --------------------------------------------------------------------------
  
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
  #  error message if argvals and y are inconsistent
  if(any(dima != dimy[1:length(dima)]))
    #  terminal message
    stop("A dimension of 'argvals' does not match 'y':\n",
         "  dim(argvals) = ", paste(dima, collapse=" x "),
         ";  dim(y) = ", paste(dimy, collapse=" x ") ) 
  
  ##. --------------------------------------------------------------------------
  ## 5.  Check compatibility of argvals with basisobj
  ##. --------------------------------------------------------------------------
  
  {
    if(inherits(basisobj, 'fd'))basisobj <- basisobj$basis
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
                #. error message if basisobj incorrect class
                stop("'basisobj' is NOT a functional basis",
                     " object (class 'basisfd');  class = ",
                     class(basisobj)[1])
          }
        }
      }
    }
  }

  ##. --------------------------------------------------------------------------
  ## 6.  Check compatibility of argvals with basisobj$rangeval
  ##. --------------------------------------------------------------------------
  
  # set up a safety zone for argvals out of range by a tiny amount
  delta <- 1e-7*(rangeval[2]-rangeval[1]) # the tiny amount
  arng  <- range(argvals)
  if ((rangeval[1]-arng[1]) <  delta) all(argvals < rangeval[1]) <- rangeval[1]
  if ((rangeval[2]-arng[2]) < -delta) all(argvals > rangeval[2]) <- rangeval[2]
  # test for argvals being out of range by more than delta
  if((rangeval[1] <= arng[1]) && (arng[2] <= rangeval[2])) {
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  } else {
    #  error message 
    stop("There are argvals not contained within basisobj$rangeval")
  }
  
}

