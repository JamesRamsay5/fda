fRegress.character <- function(y, data=NULL, betalist=NULL,
                               wt=NULL, y2cMap=NULL, SigmaE=NULL,
                               method=c('fRegress', 'model'),
                               sep='.', ...) {
  fRegress.formula(y=y, data=data, betalist=betalist,
                   wt=wt, y2cMap=y2cMap, SigmaE=SigmaE,
                   method=method, sep=sep, ...)
}

fRegress.formula <- function(y, data=NULL, betalist=NULL,
                             wt=NULL, y2cMap=NULL, SigmaE=NULL,
                             method=c('fRegress', 'model'),
                             sep='.', ...) {
  
  #  Last modified 1 November 2020 by Jim Ramsay
  
  ##
  ## 1.  get y = left hand side of the formula
  ##
  Formula <- y                    #  character vector
  yName   <- Formula[[2]]         #  name of dependent variable object
  yNm     <- as.character(yName)  #  character strong for name
  #  check name of dependent variable
  if(!inherits(yName, 'name') || (length(yNm) > 1))
    stop('The left hand side of formula must be a simple object; ',
         ' instead, LHS = ', as.character(Formula)[2],
         ', which has class = ', class(yName))
  #
  
  ##
  ## 2.  check the data argument
  ##
  
  dataNames <- names(data)
  #  extract name of dependent variable
  y <- {
    if (yNm %in% dataNames) data[[yNm]]
    else get(yNm)  # Search by name for object yNm
  }
  #  get the range of the  dependent variable  and 
  #  obtain the dimensions of the coefficient matrix if y is of class 'fd'
  #  obtain the dimensions of the y if y is of class 'numeric'
  trng <- NULL
  {
    if(inherits(y, 'fd')){
      ydim <- dim(y$coefs)
      if(is.null(ydim) || (length(ydim)<2)) {
        y$coefs <- as.matrix(y$coefs)
        ydim   <- dim(y$coefs)
      }
      ny   <- ydim[2]
      trng <- y$basis$rangeval
    } else{
      if(inherits(y, 'numeric')){
        ydim <- dim(y)
        if(is.null(ydim))
          ny <- length(y)
        else
          ny <- ydim[1]
      }
      else
        stop('The left hand side of formula must have class ',
             'numeric or fd;  instead is ', class(y))
    }
  }
  
  ##
  ## 3.  check the formula for excessive complexity
  ##
  
  allVars  <- all.vars(Formula)
  xNms     <- allVars[allVars != yNm]
  Terms    <- terms(Formula)
  termLbls <- attr(Terms, 'term.labels')
  oops     <- which(!(termLbls %in% xNms))
  if(length(oops) > 0)
    stop('formula contains terms that fRegress can not handle; ',
         ' the first one = ', termLbls[oops[1]])
  #
  k1 <- length(allVars)
  type <- rep(NA,k1)
  names(type) <- allVars
  nbasis      <- type
  if(inherits(y, 'fd')){
    type[1] <- y$basis$type
    nb      <- y$basis$nbasis
    if(!is.null(nb)) nbasis[1] <- nb
  }
  
  ##
  ## 4.  Inventory the right hand side
  ##
  
  k0       <- length(xNms)
  xfdlist0 <- vector('list', k0)
  names(xfdlist0) <- xNms
  xNames          <- xfdlist0
  nVars           <- rep(NA, k0)
  names(nVars)    <- xNms
  oops <- FALSE
  for(i in 1:k0) {
    xNm <- xNms[i]
    xi <- {
      if(xNm %in% dataNames) data[[xNm]]
      else get(xNm)
    }
    {
      if(class(xi) %in% c('fd', 'fdPar')) {
        xj <- {
          if(inherits(xi, 'fd')) xi
          else xi$fd
        }
        xrng <- xj$basis$rangeval
        {
          if(is.null(trng))
            trng <- xrng
          else
            if(any(xrng != trng)){
              oops <- TRUE
              cat('incompatible rangeval found in ', xNm,
                  '$rangeval = ', paste(xrng, collapse=', '),
                  ' != previous = ', paste(trng, collapse=', '),
                  sep='')
            }
        }
        xdim <- dim(xj$coefs)
        {
          if(is.null(xdim) || (length(xdim)<2)){
            xj$coefs <- as.matrix(xj$coefs)
            xdim <- dim(xj$coefs)
            nxi <- xdim[2]
            nVars[i] <- 1
            xNames[[i]] <- xNm
          }
          else {
            if(length(xdim)<3){
              nxi <- xdim[2]
              nVars[i] <- 1
              xNames[[i]] <- xNm
            }
            else {
              nxi <- xdim[2]
              if(length(xdim)<4){
                nVars[i] <- xdim[3]
                xNmsi <- dimnames(xj$coefs)[[3]]
                {
                  if(is.null(xNmsi))
                    xNames[[i]] <- paste(xNm, 1:xdim[3], sep=sep)
                  else
                    xNames[[i]] <- paste(xNm, xNmsi, sep=sep)
                }
              }
              else {
                oops <- TRUE
                cat(xNm, 'has too many levels:  dim(x$coefs) =',
                    paste(xdim, collapse=', '))
              }
            }
          }
        }
        type[i+1] <- xj$basis$type
        nb <- xj$basis$nbasis
        if(!is.null(nb))nbasis[i+1] <- nb
        xfdlist0[[i]] <- xi
      }
      else {
        if(is.numeric(xi)){
          xdim <- dim(xi)
          {
            if(is.null(xdim) || (length(xdim)<2)){
              nxi <- length(xi)
              nVars[i] <- 1
              xNames[[i]] <- xNm
            }
            else {
              nxi <- xdim[1]
              {
                if(length(xdim)<3){
                  nVars[i] <- xdim[2]
                  xNmsi <- dimnames(xi)[[2]]
                  {
                    if(is.null(xNmsi))
                      xNames[[i]] <- paste(xNm, 1:xdim[2], sep=sep)
                    else
                      xNames[[i]] <- paste(xNm, xNmsi, sep=sep)
                  }
                }
                else{
                  oops <- TRUE
                  cat(xNm, 'has too many levels:  dim(x) =',
                      paste(xdim, collapse=', '))
                }
              }
            }
          }
          xfdlist0[[i]] <- xi
        }
        else {
          if(inherits(xi, 'character'))
            xi <- factor(xi)
          {
            if(inherits(xi, 'factor')) {
              f.i <- formula(paste('~', xNm))
              Xi.df <- data.frame(xi)
              names(Xi.df) <- xNm
              Xi <- (model.matrix(f.i, Xi.df)[, -1, drop=FALSE])
              nxi <- dim(Xi)[1]
              xiNms <- dimnames(Xi)[[2]]
              nVars[i] <- length(xiNms)
              xNmLen <- nchar(xNm)
              xiLvls <- substring(xiNms, xNmLen+1)
              xNames[[i]] <- paste(xNm, xiLvls, sep=sep)
              xfdlist0[[i]] <- Xi
            }
            else{
              oops <- TRUE
              cat('ERROR:  variable', xNm, 'must be of class',
                  'fd, numeric, character or factor;  is', class(xi))
              nxi <- length(xi)
            }
            }
        }
      }
    }
    if(nxi != ny){
      cat('ERROR:  variable', xNm, 'has only',
          nxi, 'observations !=', ny,
          '= the number of observations of y.')
      oops <- TRUE
    }
  }
  if(oops)stop('illegal variable on the right hand side.')
  # If no functions found:
  if(is.null(trng)){
    warning("No functions found;  setting rangeval to 0:1")
    trng <- 0:1
  }
  
  ##
  ## 5.  Create xfdlist
  ##
  
  xL.L0   <- rep(1:k0, nVars)
  xNames. <- c('const', unlist(xNames))
  k <- 1+sum(nVars)
  xfdlist <- vector('list', k)
  names(xfdlist) <- xNames.
  #  create constfd for the intercept
  #  xfdlist[[1]] <- create.constant.basis(trng)
  xfdlist[[1]] <- rep(1, ny)
  i1 <- 1
  for(ix in 1:k0) {
    i0  <- i1+1
    xNm <- xNms[ix]
    xi  <- xfdlist0[[ix]]
    {
      if(inherits(xi, 'fd')) {
        if(nVars[ix]<2) {
          i1            <- i0
          xfdlist[[i0]] <- xi
        }
        else {
          #          i1 <- (i1+nVars[ix])
          for(i in 1:nVars[ix]){
            i1  <- i1+1
            xii <- xi
            xii$coefs <- xi$coefs[,,i, drop=FALSE]
            xfdlist[[i1]] <- xii
          }
        }
      }
      else {
        if(is.numeric(xi)) {
          if(nVars[ix]<2) {
            i1 <- i0
            xfdlist[[i0]] <- xi
          }
          else{
            for(i in 1:nVars[ix]) {
              i1 <- i1+1
              xfdlist[[i1]] <- xi[, i]
            }
          }
        }
      }
    }
  }
  
  ##
  ## 6.  check betalist or set up betalist
  ##
  
  {
    if(inherits(betalist, 'list')) {
      #  betalist is an argument
      if(length(betalist) != k)
        stop('length(betalist) = ', length(betalist),
             ';  must be ', k, ' to match length(xfdlist).')
      betaclass <- sapply(betalist, class)
      oops      <- which(betaclass != 'fdPar')
      if(length(oops)>0)
        stop('If betalist is a list, all components must have class ',
             'fdPar;  component ', oops[1], ' has class ',
             betaclass[oops[1]])
    }
    else {
      # betalist must be set up
      betalist <- vector('list', k)
      names(betalist) <- xNames.
      for(i in 1:k) {
        if(is.numeric(xfdlist[[i]])) {
          #  if xfdlist[[i]] is numeric, basis is set up using that  of dependent variable y
          if(is.numeric(y)) {
            bbasis        <- create.constant.basis(trng)
            bfd           <- fd(basisobj=bbasis)
            betalist[[i]] <- fdPar(bfd)
          } else {
            #  if 'fd' use the basis of dependent variable
            if(inherits(y, 'fd')) {
              bfd           <- with(y, fd(basisobj=basis, fdnames=fdnames))
              betalist[[i]] <- fdPar(bfd)
            }
            else {
              bfd           <- with(y, fd(basisobj=basis, fdnames=fdnames))
              betalist[[i]] <- with(y, fdPar(bfd, Lfd, lambda, estimate, penmat))
            }
          }
        }
        else {
          #  use basis for the independent variable
          xfdi <- {
            if(i>1) xfdlist0[[xL.L0[i-1]]]
            else    xfdlist[[1]]
          }
          if(inherits(xfdi, 'fd')){
            bfd           <- with(xfdi, fd(basisobj=basis, fdnames=fdnames))
            betalist[[i]] <- fdPar(bfd)
          }
          else{
            bfd           <- with(xfdi$fd, fd(basisobj=basis, fdnames=fdnames))
            betalist[[i]] <- with(xfdi, fdPar(bfd, Lfd, lambda,
                                              estimate, penmat))
          }
        }
      }
    }
  }
  
  ##
  ## 7.  extract or set up weight
  ##
  
  {
    if(is.null(wt))
      wt <- rep(1, ny)
    else {
      if(length(wt) != ny)
        stop('length(wt) must match y;  length(wt) = ',
             length(wt), ' != number of y observations = ', ny)
      if(any(wt<0))
        stop('Negative weights found;  not allowed.')
    }
  }
  xiEnd   <- cumsum(nVars)
  xiStart <- c(1, xiEnd[-1])
  fRegressList <- list(y=y, xfdlist=xfdlist, betalist=betalist, wt=wt)
  
  ##
  ## 8.  either output argument list for fRegress() or invoke itcs
  ##
  
  method <- match.arg(method)
  if(method=='model') {
    return(fRegressList)
  } else {
    if(inherits(y, 'fd')) {
      do.call('fRegress.fd',    fRegressList)
    } else {
      do.call('fRegress.double', fRegressList)
    }
  }
}
