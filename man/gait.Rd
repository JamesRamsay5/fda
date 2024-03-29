\name{gait}
\alias{gait}
\title{
  Hip and knee angle while walking
}
\description{
  Hip and knee angle in degrees through a 20 point movement cycle for 39
  boys  
}
\format{
   An array of dim c(20, 39, 2) giving the "Hip Angle" and "Knee Angle"
   for 39 repetitions of a 20 point gait cycle.   
}
\details{
  The components of dimnames(gait) are as follows:
  
  [[1]] standardized gait time = seq(from=0.025, to=0.975, by=0.05) 

  [[2]] subject ID = "boy1", "boy2", ..., "boy39"  

  [[3]] gait variable = "Hip Angle" or "Knee Angle"
}
\examples{
oldpar <- par(no.readonly=TRUE)
plot(gait[,1, 1], gait[, 1, 2], type="b")
par(oldpar)
}
\keyword{datasets}
