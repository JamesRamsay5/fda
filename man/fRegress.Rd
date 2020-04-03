\name{fRegress}
\alias{fRegress}
\alias{fRegress.fd}
\alias{fRegress.fdPar}
\alias{fRegress.numeric}
\alias{fRegress.formula}
\alias{fRegress.character}
\title{
  Functional Regression Analysis
}
\description{
  This function carries out a functional regression analysis, where
  either the dependent variable or one or more independent variables are
  functional.  Non-functional variables may be used on either side
  of the equation.  In a simple problem where there is a single scalar
  independent covariate with values \eqn{z_i, i=1,\ldots,N} and a single
  functional covariate with values \eqn{x_i(t)}, the two versions of the
  model fit by \code{fRegress} are the \emph{scalar} dependent variable
  model

  \deqn{y_i = \beta_1 z_i + \int x_i(t) \beta_2(t) \, dt + e_i}

  and the \emph{concurrent} functional dependent variable model

  \deqn{y_i(t) = \beta_1(t) z_i + \beta_2(t) x_i(t) + e_i(t).}

  In these models, the final term \eqn{e_i} or \eqn{e_i(t)} is a
  residual, lack of fit or error term.

  In the concurrent functional linear model for a functional dependent
  variable, all functional   variables are all evaluated at a common
  time or argument value $t$.  That is, the fit is defined in terms of
  the behavior of all variables at a fixed time, or in terms of "now"
  behavior.

  All regression coefficient functions \eqn{\beta_j(t)} are considered
  to be functional.  In the case of a scalar dependent variable, the
  regression coefficient for a scalar covariate is converted to a
  functional variable with a constant basis.   All regression
  coefficient functions can be forced to be \emph{smooth} through the
  use of roughness penalties, and consequently are specified in the
  argument list as functional parameter objects.
}
\usage{
fRegress(y, xfdlist, betalist, wt=NULL, y2cMap=NULL, SigmaE=NULL,
         method, sep='.', ...)
\method{fRegress}{fd}(y, xfdlist, betalist, wt=NULL,
                  y2cMap=NULL, SigmaE=NULL, method, sep='.', ...)
\method{fRegress}{fdPar}(y, xfdlist, betalist, wt=NULL,
                  y2cMap=NULL, SigmaE=NULL, method, sep='.', ...)
\method{fRegress}{numeric}(y, xfdlist, betalist, wt=NULL,
                  y2cMap=NULL, SigmaE=NULL, method, sep='.', ...)
\method{fRegress}{formula}(y, xfdlist=NULL, betalist=NULL, wt=NULL,
                  y2cMap=NULL, SigmaE=NULL, method=c('fRegress', 'model'), 
                  sep='.', ...)
\method{fRegress}{character}(y, xfdlist=NULL, betalist=NULL, wt=NULL,
                  y2cMap=NULL, SigmaE=NULL, method=c('fRegress', 'model'), 
                  sep='.', ...)
}
\arguments{
  \item{y}{the dependent variable object.  It may be an object of five
  possible classes:
    \itemize{
      \item{character or formula}{
	    a \code{formula} object or a \code{character} object that can be
	    coerced into a \code{formula} providing a symbolic description
	    of the model to be fitted satisfying the following rules:

	    The left hand side, \code{formula} \code{y}, must be either a
	    numeric vector or a univariate object of class \code{fd} or
	    \code{fdPar}.  If the former, it is replaced by 
	    \code{fdPar(y, ...)}.

	    All objects named on the right hand side must be either
	    \code{numeric} or \code{fd} (functional data) or \code{fdPar}.
	    The number of replications of \code{fd} or \code{fdPar}
	    object(s) must match each other and the number of observations
	    of \code{numeric} objects named, as well as the number of
	    replications of the dependent variable object.  The right hand
	    side of this \code{formula} is translated into \code{xfdlist},
	    then passed to another method for fitting (unless \code{method}
	    = 'model'). Multivariate independent variables are allowed in a
	    \code{formula} and are split into univariate independent
	    variables in the resulting \code{xfdlist}.  Similarly,
	    categorical independent variables with \code{k} levels are
	    translated into \code{k-1} contrasts in \code{xfdlist}.  Any
	    smoothing information is passed to the corresponding component
	    of \code{betalist}.
      }
      \item{scalar}{
	    a vector if the dependent variable is scalar.
      }
      \item{fd}{
	    a functional data object if the dependent variable is
	    functional.  
      }
      \item{fdPar}{
	    a functional parameter object if the dependent variable is
	    functional, and if it is desired to smooth the prediction
	    of the dependent variable.
      }
    }
  }
  \item{xfdlist}{
    a list of length equal to the number of independent
    variables (including any intercept). Members of this list are the
    independent variables.  They can be objects of either of these two
    classes:

    \itemize{
      \item{scalar}{
	    a numeric vector if the independent variable is scalar.
      }
      \item{fd}{
	    a (univariate) functional data object.
      }
    }
    In either case, the object must have the same number of replications
    as the dependent variable object.  That is, if it is a scalar, it
    must be of the same length as the dependent variable, and if it is
    functional, it must have the same number of replications as the
    dependent variable.  (Only univariate independent variables are
    currently allowed in \code{xfdlist}.)
  }
  \item{betalist}{
    For the \code{fd}, \code{fdPar}, and \code{numeric} methods,
    \code{betalist} must be a list of length equal to
    \code{length(xfdlist)}.  Members of this list are functional
    parameter objects (class \code{fdPar}) defining the regression
    functions to be estimated.  Even if a corresponding independent
    variable is scalar, its regression coefficient must be functional if
    the dependent variable is functional.  (If the dependent variable is
    a scalar, the coefficients of scalar independent variables,
    including the intercept, must be constants, but the coefficients of
    functional independent variables must be functional.)  Each of these
    functional parameter objects defines a single functional data
    object, that is, with only one replication.

    For the \code{formula} and \code{character} methods, \code{betalist}
    can be either a \code{list}, as for the other methods, or
    \code{NULL}, in which case a list is created.  If \code{betalist} is
    created, it will use the bases from the corresponding component of
    \code{xfdlist} if it is function or from the response variable.
    Smoothing information (arguments \code{Lfdobj}, \code{lambda},
    \code{estimate}, and \code{penmat} of function \code{fdPar}) will
    come from the corresponding component of \code{xfdlist} if it is of
    class \code{fdPar} (or for scalar independent variables from the
    response variable if it is of class \code{fdPar}) or from optional
    \code{\dots} arguments if the reference variable is not of class
    \code{fdPar}.
  }
  \item{wt}{
    weights for weighted least squares
  }
  \item{y2cMap}{
    the matrix mapping from the vector of observed values to the
    coefficients for the dependent variable.  This is output by function
    \code{smooth.basis}.  If this is supplied, confidence limits are
    computed, otherwise not.
  }
  \item{SigmaE}{
    Estimate of the covariances among the residuals.  This can only be
    estimated after a preliminary analysis with \code{fRegress}.
  }
  \item{method}{
    a character string matching either \code{fRegress} for functional
    regression estimation or \code{mode} to create the argument lists
    for functional regression estimation without running it.
  }
  \item{sep}{
    separator for creating names for multiple variables created from single
    variables on the right hand side of the \code{formula} \code{y}.
    This happens with multidimensional \code{fd} objects as well as with
    categorical variables.
  }
  \item{\dots}{ optional arguments }
}
\details{
  Alternative forms of functional regression can be categorized with
  traditional least squares using the following 2 x 2 table:

  \tabular{lcccc}{
    \tab \tab explanatory \tab variable \tab \cr
    response \tab | \tab scalar \tab  | \tab function \cr
    \tab | \tab \tab | \tab \cr
    scalar   \tab | \tab lm     \tab   | \tab fRegress \cr
    \tab | \tab \tab | \tab \cr
    function \tab | \tab fRegress \tab | \tab fRegress or linmod \cr
  }

  \code{linmod} predicts a functional response from a convolution
  integral, estimating a bivariate regression function.


  In the computation of regression function estimates in
  \code{fRegress}, all independent variables are treated as if they are
  functional.  If argument \code{xfdlist} contains one or more vectors,
  these are converted to functional data objects having the constant
  basis with coefficients equal to the elements of the vector.

  Needless to say, if all the variables in the model are scalar, do NOT
  use this function.  Instead, use either \code{lm} or \code{lsfit}.

  These functions provide a partial implementation of Ramsay and
  Silverman (2005, chapters 12-20).
}
\value{
  These functions return either a standard \code{fRegress} fit object or
  or a model specification:

  \item{fRegress fit}{
    a list of class \code{fRegress} with the following components:

    \itemize{
      \item{y}{
	the first argument in the call to \code{fRegress} (coerced to
	\code{class} \code{fdPar})

      }
      \item{xfdlist}{
	the second argument in the call to \code{fRegress}.

      }
      \item{betalist}{
	the third argument in the call to \code{fRegress}.

      }
      \item{betaestlist}{
	a list of length equal to the number of independent variables
	and with members having the same functional parameter structure
	as the corresponding members of \code{betalist}.  These are the
	estimated regression coefficient functions.

      }
      \item{yhatfdobj}{
	a functional parameter object (class \code{fdPar}) if the
	dependent variable is functional or a vector if the dependent
	variable is scalar.  This is the set of predicted by the
	functional regression model for the dependent variable.

      }
      \item{Cmatinv}{
	a matrix containing the inverse of the coefficient matrix for
	the  linear equations that define the solution to the regression
	problem.  This matrix is required for function
	\code{\link{fRegressStderr}} that estimates confidence regions
	for the regression coefficient function estimates.

      }
      \item{wt}{
	the vector of weights input or inferred

      }


   If \code{class(y)} is numeric, the \code{fRegress} object also includes:

      \item{df}{
	equivalent degrees of freedom for the fit.
      }
      \item{OCV}{
      the leave-one-out cross validation score for the model.
      }
      \item{gcv}{
      the generalized cross validation score.
      }


      If \code{class(y)} is either \code{fd} or \code{fdPar}, the
      \code{fRegress} object returned also includes 5 other components:

      \item{y2cMap}{
	an input \code{y2cMap}

      }
      \item{SigmaE}{
	an input \code{SigmaE}

      }
      \item{betastderrlist}{
	an \code{fd} object estimating the standard errors of
	\code{betaestlist}

      }
      \item{bvar}{
	a covariance matrix

      }
      \item{c2bMap}{	a map
      }
    }
  }

  \item{model specification}{
    The \code{fRegress.formula} and \code{fRegress.character}
    functions translate the \code{formula} into the argument list
    required by \code{fRegress}.
    With the default value 'fRegress' for the argument \code{method},
    this list is then used to call the appropriate other
    \code{fRegress} function.

    Alternatively, to see how the \code{formula} is translated, use
    the alternative 'model' value for the argument \code{method}.  In
    that case, the function returns a list with the arguments
    otherwise passed to these other functions plus the following
    additional components:

    \itemize{
      \item{xfdlist0}{
	a list of the objects named on the right hand side of
	\code{formula}.  This will differ from \code{xfdlist} for
	any categorical or multivariate right hand side object.

      }
      \item{type}{
	the \code{type} component of any \code{fd} object on the right
	hand side of \code{formula}.

      }
      \item{nbasis}{
	a vector containing the \code{nbasis} components of variables
	named in \code{formula} having such components

      }
      \item{xVars}{
	an integer vector with all the variable names on the right
	hand side of \code{formula} containing the corresponding
	number of variables in \code{xfdlist}.  This can exceed 1 for
	any multivariate object on the right hand side of class either
	\code{numeric} or \code{fd} as well as any categorical
	variable.

      }
    }
  }
}
\author{
  J. O. Ramsay, Giles Hooker, and Spencer Graves
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009)
    \emph{Functional Data Analysis in R and Matlab}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\seealso{
  \code{\link{fRegressStderr}},
  \code{\link{fRegressCV}},
  \code{\link{linmod}}
}
\examples{
##
## set up the regression of log annual precipitation on temperature functions
##

# compute annual precipitations for 35 weather stations
annualprec <- apply(CanadianWeather$dailyAv[,,"Precipitation.mm"], 2,sum)
# log precipitation
logannualprec <- log10(annualprec)
# We use a small basis because there are only 35 stations
tempbasis  <- create.fourier.basis(c(0, 365), 7)
#  smooth the temperature data with this basis
tempfd <- smooth.basis(day.5,
          CanadianWeather$dailyAv[,,"Temperature.C"], tempbasis)$fd

##
## Manual construction of arguments xfdlist and betalist
##

#  single functional covariate, tempfd

xfdlist <- list(const=rep(1, 35), tempfd=tempfd)

# The intercept must be constant for a scalar response

betabasis1 <- create.constant.basis(c(0, 365))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)

# the temperature regression coefficient

betafd2     <- with(tempfd, fd(basisobj=basis, fdnames=fdnames))
# convert to an fdPar object
betafdPar2  <- fdPar(betafd2)

#  store the constant and temperature functions in a list

betalist <- list(const=betafdPar1, tempfd=betafdPar2)

#  execute the functional regression to get the list of results

precip.temp <- fRegress(logannualprec, xfdlist, betalist)

names(precip.temp)

#  plot the two regression functions

par(ask=TRUE)
plotbeta(precip.temp$betaestlist) 

#  plot log annual precipitation and its fit

par(ask=FALSE)
plot(precip.temp$yhat, precip.temp$yhat, type="l", 
     xlab="log annual precipitation fit",
     ylab="log annual precipitation")
points(precip.temp$yhat, logannualprec, pch="o")

##
## use the simpler setup using the formula interface
##

precip.temp <- fRegress(logannualprec ~ tempfd)

###
### functional response with vector explanatory variables
###

#  set up prediction of temperature from climate region

region <- CanadianWeather$region
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65,
                  axes=list('axesIntervals'))
Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
                 dailyAv[,,'Temperature.C'], daybasis65)$fd)
region.contrasts <- model.matrix(~factor(CanadianWeather$region))
rgnContr3 <- region.contrasts
dim(rgnContr3) <- c(1, 35, 4)
dimnames(rgnContr3) <- list('', CanadianWeather$place, c('const',
   paste('region', c('Atlantic', 'Continental', 'Pacific'), sep='.')) )

const365 <- create.constant.basis(c(0, 365))
region.fd.Atlantic    <- fd(matrix(rgnContr3[,,2], 1), const365)
region.fd.Continental <- fd(matrix(rgnContr3[,,3], 1), const365)
region.fd.Pacific     <- fd(matrix(rgnContr3[,,4], 1), const365)
region.fdlist <- list(const=rep(1, 35),
     region.Atlantic=region.fd.Atlantic,
     region.Continental=region.fd.Continental,
     region.Pacific=region.fd.Pacific)

beta1 <- with(Temp.fd, fd(basisobj=basis, fdnames=fdnames))
beta0 <- fdPar(beta1)
betalist <- list(const=beta0, region.Atlantic=beta0,
             region.Continental=beta0, region.Pacific=beta0)

TempRgn <- fRegress(Temp.fd, region.fdlist, betalist)

#  plot regression functions for baseline, Atlantic, Continental and Pacific

par(ask=TRUE)
plotbeta(TempRgn$betaestlist)

#  plot estimated fits for these regions

par(ask=FALSE)
plot(TempRgn$yhatfdobj)

##
## simpler:  formula interface
##

TempRgn.f <- fRegress(Temp.fd ~ region)

par(ask=FALSE)
plot(TempRgn.f$yhatfdobj)

###
### functional response with
###            (concurrent) functional explanatory variable
###

#  predict knee angle from hip angle;  

gaittime  <- (0:19) + 0.5
gaitrange <- c(0,20)
gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)
harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd, lambda=1e-2)$fd
hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]

# manual set-up

const  <- rep(1, dim(kneefd$coef)[2])
xfdlist  <- list(const=const, hipfd=hipfd)

beta0 <- with(kneefd, fd(basisobj=basis, fdnames=fdnames))
beta1 <- with(hipfd,  fd(basisobj=basis, fdnames=fdnames))

betalist  <- list(const=fdPar(beta0), hipfd=fdPar(beta1))

fRegressout <- fRegress(kneefd, xfdlist, betalist)

# formula interface

knee.hip.f <- fRegress(kneefd ~ hipfd)

#See also the following demos:

#demo('canadian-weather', package='fda')
#demo('gait', package='fda')
#demo('refinery', package='fda')
#demo('weatherANOVA', package='fda')
#demo('weatherlm', package='fda')
}
\keyword{smooth}
