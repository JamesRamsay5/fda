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
  
  argChk <- argvalsySwap(argvals, y, basisobj)
  
  if(!is.numeric(AV <- argChk$argvals)){
    print(AV)
    if(is.null(AV))
      stop('is.null(argChk$argvals); should be numeric')
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
  ## 1.  if(is.null(y)) use argvals for y
  ##. --------------------------------------------------------------------------
  
  if(is.null(y)){
    if(is.null(argvals)) stop("'y' is missing with no default")
    #   Store argvals as y and alert
    cat("'y' is missing, using 'argvals'\n") 
    y <- argvals
    argvals <- NULL 
  }
  
  ##. --------------------------------------------------------------------------
  ## 2.  test for missing argvals, if so construct a sequence
  ##. --------------------------------------------------------------------------

  dimy <- dim(as.array(y))
  if(is.null(argvals)) {
    # the following code block is run if TRUE
    {  # beginning of code block
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
    }  #. end of code block
    # This is executed whether or not the previous was
    #. locate the range from basisobj
    a01 <- basisobj$rangeval
    #  if range is null, error message and stop
    if(is.null(a01))
      stop('basisobj does not have a required rangeval component.')
    n <- dimy[1]
    #  construct the argval sequence
    argvals <- seq(a01[1], a01[2], length=n)
    #  warning message about the swap
    cat(paste("'argvals' is missing;  using seq(", a01[1],
              ", ", a01[2], ", length=", n, ")\n"))
    #. return
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  
  ##. --------------------------------------------------------------------------
  ## 3.  swapping y and argvals 
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
  for (i in 1:length(argvals)) {
    argi <- argvals[i]
    # print(argi)
    if (argi < rangeval[1] && argi >= rangeval[1]-delta) {
      argi <- rangeval[1]
    }
    if (argi > rangeval[2] && argi <= rangeval[2]+delta) {
      argi <- rangeval[2]
    }
  }
  if (any(argvals < rangeval[1])  || any(argvals > rangeval[2])) {
    #  error message 
    stop("There are argvals not contained within interval basisobj$rangeval")
  }
  return(list(argvals=argvals, y=y, basisobj=basisobj))
  
}

