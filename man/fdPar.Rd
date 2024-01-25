\name{fdPar}
\alias{fdPar}
\title{
  Define a Functional Parameter Object
}
\description{
  A functional parameter object is a functional data object such as \eqn{x(t)}
  with an additional object defining a roughness penalty and a smoothing 
  parameter.  The roughness penalty times the smoothing parameter added to
  the error sum of squares for the fit of the functional data object to the
  data results in the fitting function \eqn{x(t)} being smoother than it would
  be if it were defined only its basis functions.
  
  A common case for illustrating smoothing with a roughness penalty 
  combines an error sum
  of squares \deqn{\mbox{SSE}(x) = \sum_{j=1}^n [y_j - x(t_j)]^2} with
  a roughness measure for \eqn{x} defined by its
  second derivative \eqn{D^2x}
  \deqn{\mbox{PEN}(x) = \int [D^2x(t)]^2 d\mbox{t}} \eqn{\mbox{PEN}(x)} 
  multiplied by a positive scalar \eqn{\lambda} defines the composite error 
  sum of squares measure
  \deqn{\mbox{F}(x) = \mbox{SSE}(x) + \lambda \mbox{PEN}(x)}  
  
  Minimizing \eqn{F(x)} results in a satisfactory fit to the data that is also 
  smooth in the sense that its second derivative has smaller variation than if 
  only \eqn{\mbox{SSE}(x)} were minimized.  The larger the smoothing parameter
  \eqn{\lambda} is, the smoother the second derivative of \eqn{x} will be.
  
  See  Ch. 5 in R&S or 5.2 in RH&G  cited below.
  
}
\usage{
fdPar(fdobj=NULL, Lfdobj=NULL, lambda=0, estimate=TRUE,
      penmat=NULL)
%fdPar(fdobj=fd(), Lfdobj=NULL, lambda=0, estimate=TRUE, penmat=NULL)
}
\arguments{
  \item{fdobj}{
    a functional data object, functional basis object, a functional
    parameter object or a matrix.  If it a matrix, it is replaced by
    fd(fdobj).  If class(fdobj) == 'basisfd', it is converted to an
    object of class \code{fd} with a coefficient matrix consisting of a
    single column of zeros.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.

    If \code{NULL}, Lfdobj depends on fdobj[['basis']][['type']]:

    \describe{
      \item{bspline}{
	Lfdobj <- int2Lfd(max(0, norder-2)), where norder =
	norder(fdobj).
      }
      \item{fourier}{
	Lfdobj = a harmonic acceleration operator:

	\code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)}

	where rng = fdobj[['basis']][['rangeval']].
      }
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{lambda}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter.
  }
  \item{estimate}{not currently used.}
  \item{penmat}{
    a roughness penalty matrix.  Including this can eliminate the need
    to compute this matrix over and over again in some types of
    calculations.
  }
}
\details{
  Choosing \eqn{\lambda} will often involve trying a range of values and
  evaluating each fit by eye.  However,there are quantitive evaluations of 
  roughness penalized measures of fit that can be minimized, such as 
  generalized cross-validation (GCV).  See R&S and RH&G for details.
}
\value{
  a functional parameter object (i.e., an object of class \code{fdPar}),
  which is a list with the following components:

  \item{fd}{
    a functional data object (i.e., with class \code{fd})
  }
  \item{Lfd}{
    a linear differential operator object (i.e., with class
    \code{Lfd})
  }
  \item{lambda}{a nonnegative real number}
  \item{estimate}{a logical value that can be dropped because it is not  
  currently used.}
  \item{penmat}{Either NULL or a square, symmetric matrix with penmat[i,j] =
    integral over fd[['basis']][['rangeval']] of basis[i]*basis[j].
    if NULL, the penalty matrix will be set up.  Computing it will
    speed up computation if multiple analyses are envisaged.}
}
\source{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009), 
    Functional data analysis with R and Matlab, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York
}
\seealso{
  \code{\link{cca.fd}},
  \code{\link{density.fd}},
  \code{\link{fRegress}},
  \code{\link{intensity.fd}},
  \code{\link{pca.fd}},
  \code{\link{smooth.fdPar}},
  \code{\link{smooth.basis}},
  \code{\link{smooth.monotone}},
  \code{\link{int2Lfd}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## Simple example
##
#  set up range for density
rangeval <- c(-3,3)
#  set up some standard normal data
x <- rnorm(50)
#  make sure values within the range
x[x < -3] <- -2.99
x[x >  3] <-  2.99
#  set up basis for W(x)
basisobj <- create.bspline.basis(rangeval, 11)
#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0,11,1), basisobj)
WfdParobj <- fdPar(Wfd0)

##
##  smooth the Canadian daily temperature data
##
#    set up the fourier basis
nbasis   <- 365
dayrange <- c(0,365)
daybasis <- create.fourier.basis(dayrange, nbasis)
dayperiod <- 365
harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), dayrange)
#  Make temperature fd object
#  Temperature data are in 365 by 12 matrix tempav
#    See analyses of weather data.
#  Set up sampling points at mid days
daytime  <- (1:365)-0.5
#  Convert the data to a functional data object
daybasis365 <- create.fourier.basis(dayrange, nbasis, dayperiod)
dayfd365    <- fd(matrix(0,365,12), daybasis365)
templambda  <- 1e1
tempfdPar   <- fdPar(dayfd365, Lfdobj=harmaccelLfd, lambda=templambda)
\dontrun{
daytempfd <- with(CanadianWeather, 
                  smooth.basis(daytime, dailyAv[,,"Temperature.C"], 
                  tempfdPar)$fd)
}
par(oldpar)

}

\keyword{smooth}
