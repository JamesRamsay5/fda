#  test problems for fRegress and fRegress.stderr for concurrent model

#  This file can't be put in a test subdirectory, but could be put in demo.

#  Last modified 3 August by Jim Ramsay

# library(fda)
# 
# setwd("R")
# 
# source("fRegress.R")
# source("fRegress.stderr.R")

#  --------------------------------------------------------------------
#  Concurrent functional linear model:
#  Two covariates,  monomial beta, periodic covariate, domain c(-pi,pi)
#  --------------------------------------------------------------------

trng = c(-pi,pi)

#  bases for covariates and regression functions

Xbasis = create.fourier.basis(trng, 3)   #  order 3 fourier basis
Bbasis = create.monomial.basis(trng, 3)  #  order 3 monomial basis

par(mfrow = c(2,1))
par(ask=FALSE)
plot(Xbasis)
plot(Bbasis)

#  set up two covariate objects

N = 5  #  number of replications

theta = seq(-pi, pi, len=N)

#  first covariate:  sine function 2nd coefficient

Xcoef1 = matrix(0,3,N)
Xcoef1[2,] = sin(theta/2)
Xfd1 = fd(Xcoef1, Xbasis)

#  second covariate:  cosine function 2nd coefficient

Xcoef2 = matrix(0,3,N)
Xcoef2[2,] = cos(theta/2)
Xfd2 = fd(Xcoef2, Xbasis)

#  plot covariates: sine and cosines with varying amplitudes

nfine = 51
tfine = seq(-pi,pi,len=nfine)

Xmat1 = eval.fd(tfine, Xfd1)
Xmat2 = eval.fd(tfine, Xfd2)

par(mfrow=c(2,1),ask=FALSE)
matplot(tfine, Xmat1, type="l")
matplot(tfine, Xmat2, type="l")

#  set up two regression functions, first linear, second quadratic

Bcoef1 = matrix(c(1, 1/pi, 0     ),3,1)
Bcoef2 = matrix(c(0, 0,    1/pi^2),3,1)
Bfd1   = fd(Bcoef1, Bbasis)
Bfd2   = fd(Bcoef2, Bbasis)

#  plot the regression functions

par(ask=FALSE)
plot(Bfd1)
plot(Bfd2)

#  set up the functional dependent variable for the concurrent model

#  basis 11 order four bspline functions

ynbasis = 11
Ybasis  = create.bspline.basis(trng, ynbasis)

#  evaluate the right side over 26 observation points

nobs = 26
tobs = seq(-pi,pi,len=nobs)

#  Evaluate the model values over tobs

Xtru1 = eval.fd(tobs, Xfd1)
Btru1 = eval.fd(tobs, Bfd1)
Ytru1 = (Btru1 %*% matrix(1,1,N)) * Xtru1
Xtru2 = eval.fd(tobs, Xfd2)
Btru2 = eval.fd(tobs, Bfd2)
Ytru2 = (Btru2 %*% matrix(1,1,N)) * Xtru2
Ytru = Ytru1 + Ytru2

#  define the true curve as a fd object using the 
#  functional version of the * operator

Ytrufd = Xfd1*Bfd1 + Xfd2*Bfd2

#  add Gaussian noise to model value

sigma = 0.1
Yobs = Ytru + matrix(rnorm(nobs*N),nobs,N)*sigma

#  display estimated noise variance

Yres = Yobs - Ytru
Yvar = sum(Yres^2)/(N*(nobs-ynbasis))

print(sqrt(Yvar))

#  get estimate of error variance 

SigmaE = Yvar*diag(rep(1,nobs))

#  smooth the noisy values to get approximate dependent variable

YsmthList = smooth.basis(tobs, Yobs, Ybasis)
Ysmthfd   = YsmthList$fd
Ysmthfine  = eval.fd(tfine, Ysmthfd)

#  get y2cMap mapping original noisy data into functional coefs

y2cMap = YsmthList$y2cMap

#  compute and plot standard error curve for data smooths

Ybasisobs  = eval.basis(tobs, Ybasis)
YhatVarobs = Ybasisobs %*% y2cMap %*% SigmaE %*% t(y2cMap) %*% t(Ybasisobs)
YhatStderr = matrix(sqrt(diag(YhatVarobs)),nobs,1)

par(ask=FALSE)
plot(tobs, YhatStderr, type="l", lty=1)

#  display the data, the smooth estimates, the true curve values
#  and 95% pointwise confidence bands

Ysmthobs = eval.fd(tobs, Ysmthfd)
par(ask=TRUE)
par(mfrow=c(1,1))
for (i in 1:N) {
  plot(tfine, Ysmthfine[,i], type="l", lty=1, ylim=c(-1,1))
  lines(c(-pi,pi),c(0,0), lty=3)
  lines(tobs, Ysmthobs[,i] + 2*YhatStderr, lty=2)
  lines(tobs, Ysmthobs[,i] - 2*YhatStderr, lty=2)
  points(tobs, Yobs[,i], pch="o", col=4)
}

#  now set up the fRegress analysis

#  covariate list

xfdlist = vector("list",2)
xfdlist[[1]] = Xfd1
xfdlist[[2]] = Xfd2

#  regression function list

betalist = vector("list",2)
betalist[[1]] = Bfd1
betalist[[2]] = Bfd2

#  fRegress fit

fRegressResult <- fRegress(Ysmthfd, xfdlist, betalist)

#  get and display estimated fit to functional dependent variable

yhatfdobj = fRegressResult$yhatfdobj

par(ask=FALSE)
plot(yhatfdobj)

#  use fRegress.predict with original covariates:  should 
#  return unchanged estimated fit
#  pointwise standard error of estimate not possible
#  at this point because fRegress.stderr hasn't been run.

object   = fRegressResult
newdata  = xfdlist 

yhatobj2 = predict.fRegress(fRegressResult, xfdlist)

#  extract  and display regression coefficients

betaestlist = fRegressResult$betaestlist
betafdobj1 = betaestlist[[1]]
betafdobj2 = betaestlist[[2]]

par(ask=FALSE)
par(mfrow=c(2,1))
plot(betafdobj1)
plot(betafdobj2)

#  get standard errors for both regression coefficients and for
#  fit object yhatfdobj.  These  and other objects are added to
#  the fRegress class object produced by fRegress

fRegressResult = fRegress.stderr(fRegressResult, y2cMap, SigmaE)

YhatStderr = fRegressResult$YhatStderr

#  plot five standard error curves

nplot = dim(YhatStderr)[1]
tplot = seq(-pi,pi,len=nplot)

par(ask=FALSE)
matplot(tplot, YhatStderr, type="l")    

#  get the five covariates and regression functions at plotting points

yhatfdplot = eval.fd(tplot, yhatfdobj)

#  plot the fit and its 95% confidence limits, as well as the data

par(ask=TRUE)
par(mfrow=c(1,1))
for (i in 1:N) {
  plot (tplot, yhatfdplot[,i], type="l", col=1, ylim=c(-1,1))
  lines(tplot, yhatfdplot[,i]+2*YhatStderr[,i],col=2,lty=2)
  lines(tplot, yhatfdplot[,i]-2*YhatStderr[,i],col=2,lty=2)
  points(tobs, Yobs[,i], col=4, pch="o")
  lines(c(-pi,pi),c(0,0), lty=3)
  title(paste("i",i))
}

#  now use predict.fRegress to get pointwise standard error of fit

predictResult = predict.fRegress(fRegressResult, xfdlist, se.fit = TRUE)

yhatfdplot = eval.fd(tplot,predictResult$pred)
YhatStderr = predictResult$YhatStderr
