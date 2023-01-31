\name{plotreg.fd}
\alias{plotreg.fd}
\title{
  Plot the results of the registration of a set of curves
}
\description{
A function is said to be aligned or registered with a target function if its
salient features, such as peaks, valleys and crossings of fixed thresholds,
occur at about the same argument values as those of the target.
Function \code{plotreg.fd} plots for each curve that is registered (1)
the unregistered curve (blue dashed line), (2) the target curve (red dashed
line) and (3) the registered curve (blue solid line).  It also plots within
the same graphics window the warping function $h(t)$ along with a dashed
diagonal line as a comparison.
}
\usage{
plotreg.fd(reglist)
}
\arguments{
  \item{reglist}{
    a named list that is output by a call to function \code{register.fd}.
    The members of \code{reglist} that are required are:
    \code{regfd}  ... the registered functions,
    \code{Wfd}    ... the functions W(t) defining the warping functions h(t),
    \code{yfd}    ... the unregistered functions, and
    \code{y0fd}   ... the target functions.
  If required objects are missing, REGLIST was probably generated by
  an older verson of REGISTER.FD, and the registration should be redone.
  }
}
\value{
  a series of plots, each containing two side-by-side panels.  Clicking
  on the R Graphics window advances to the next plot.
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 6 & 7.
}
\seealso{
  \code{\link{smooth.monotone}},
  \code{\link{smooth.morph}}
  \code{\link{register.fd}}
}
\examples{
if (!CRAN()) {
#  register and plot the angular acceleration of the gait data
gaittime  <- seq(0.05, 0.95, 0.1)*20
gaitrange <- c(0,20)
#  set up a fourier basis object
gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)
#  set up a functional parameter object penalizing harmonic acceleration
harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitfdPar    <- fdPar(gaitbasis, harmaccelLfd, 1e-2)
#  smooth the data
gaitfd <- smooth.basis(gaittime, gait, gaitfdPar)$fd
#  compute the angular acceleration functional data object
D2gaitfd      <- deriv.fd(gaitfd,2)
names(D2gaitfd$fdnames)[[3]] <- "Angular acceleration"
D2gaitfd$fdnames[[3]] <- c("Hip", "Knee")
#  compute the mean angular acceleration functional data object
D2gaitmeanfd  <- mean.fd(D2gaitfd)
names(D2gaitmeanfd$fdnames)[[3]] <- "Mean angular acceleration"
D2gaitmeanfd$fdnames[[3]] <- c("Hip", "Knee")
#  register the functions for the first 10 boys
#  argument periodic = TRUE causes register.fd to estimate a horizontal shift
#  for each curve, which is a possibility when the data are periodic
nBoys <- 2 # use only 2 boys to save test time.
#  set up the basis for the warping functions
nwbasis   <- 7
wbasis    <- create.bspline.basis(gaitrange,nwbasis,3)
Warpfd    <- fd(matrix(0,nwbasis,nBoys),wbasis)
WarpfdPar <- fdPar(Warpfd)
#  carry out the continuous registration
gaitreglist <- register.fd(D2gaitmeanfd, D2gaitfd[1:nBoys], WarpfdPar,
                           iterlim=4, periodic=TRUE)
# set iterlim=4 to reduce the compute time;
# this argument may not be needed in many applications.
#  plot the results
plotreg.fd(gaitreglist)
#  display horizonal shift values
print(round(gaitreglist$shift,1))
}
}
\keyword{smooth}
