\name{plot.Lfd}
\alias{plot.Lfd}
\title{
  Plot a Linear Differential Operator Object
}
\description{
  Plot the coefficents of the terms of order 0 through \code{m-1} of an
  object of class \code{Lfd} and length \code{m}.
}
\usage{
\method{plot}{Lfd}(x, axes=NULL, ...)
}
\arguments{
  \item{x}{
    a linear differential operator object to be plotted.
  }
  \item{axes}{
    Either a logical or a list or \code{NULL} passed to \code{plot.fd}.

    \describe{
      \item{logical}{
	whether axes should be drawn on the plot
      }
      \item{list}{
	a list used to create custom \code{axes} used to create axes via
	\code{x$axes[[1]]} and  \code{x$axes[-1]}.  The primary example of
	this uses \code{list("axesIntervals", ...)}, e.g., with
	\code{Fourier} bases to create \code{CanadianWeather} plots
      }
    }
  }
  \item{\dots }{
    additional plotting arguments that can be used with function
    \code{plot}
  }
}
%\details{ }
\value{
  invisible(NULL)
}
\section{Side Effects}{
  a plot of the linear differential operator object.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{Lfd}}, \code{\link{plot.fd}}
}
\examples{
#  Set up the harmonic acceleration operator
dayrange  <- c(0,365)
Lbasis  <- create.constant.basis(dayrange,
                  axes=list("axesIntervals"))
Lcoef   <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj  <- fd(Lcoef,Lbasis)
bwtlist <- fd2list(bfdobj)
harmaccelLfd <- Lfd(3, bwtlist)
oldpar <- par(no.readonly=TRUE)
plot(harmaccelLfd)
par(oldpar)
}
\keyword{smooth}
\keyword{hplot}
