#  test problems for fRegress and fRegress.stderr for scalar on function model

#  Last modified 3 August by Jim Ramsay

# library(fda)
# 
# setwd("R")
# source("fRegress.R")
# source("fRegress.stderr.R")

#  One covariate,  monomial beta, periodic covariate

#  set up rangeval and fine mesh for plotting
trng = c(-pi,pi)
nfine = 101
tfine = seq(-pi,pi,len=nfine)

#  the number of scalar values in the dependent variable

N = 51
tobs = seq(-pi, pi, len=N)

#  the bspline basis for the covariates

nxbasis = 16
Xbasis = create.bspline.basis(trng, nxbasis)

#  the coefficient matrices for the covariates

Xcoef1 = matrix(rnorm(N*nxbasis),nxbasis,N)
Xcoef2 = matrix(rnorm(N*nxbasis),nxbasis,N)

#  the fd objects for the covariates

Xfd1 = fd(Xcoef1,Xbasis)
Xfd2 = fd(Xcoef2,Xbasis)

#  set up and plot the monomial basis for the regression coefficient functions

nbbasis = 3
Bbasis = create.monomial.basis(trng, nbbasis)

#  the random values for the regression coefficient matrices

Bcoef = matrix(rnorm(6),3,2)

#  the regression coefficient functions

Bfd1 = fd(Bcoef[,1], Bbasis)
Bfd2 = fd(Bcoef[,2], Bbasis)

#  the list containing covariate functions

xfdlist <- vector("list", 2)
xfdlist[[1]] <- Xfd1
xfdlist[[2]] <- Xfd2

#  the list containing regression functions

betalist = vector("list",2)
betalist[[1]] <- Bfd1
betalist[[2]] <- Bfd2

#  compute the true or errorless value of the model

XBbasismat = inprod(Xbasis, Bbasis)
Yvectru <- t(Xcoef1) %*% XBbasismat %*% Bcoef[,1] + t(Xcoef2) %*% XBbasismat %*% Bcoef[,2]

#  standard deviation of error or residual values

sigma = 2

#  add random errors to true values

Yvec <- as.vector(Yvectru + rnorm(N) * sigma)

#  Analyze the data

fRegressResult <- fRegress(Yvec, xfdlist, betalist)

#  extract the estimated model values

Yhatvec = fRegressResult$yhatfdobj

#  plot the estimated values against the true values

par(ask=FALSE)
plot(Yvectru, Yhatvec)
lines(Yvectru, Yvectru)

#  extract and plot the estimaed values of regression functions

betaestlist = fRegressResult$betaestlist
Bhatfd1 = betaestlist[[1]]
Bhatfd2 = betaestlist[[2]]

Bvec1 = eval.fd(tfine, Bfd1)
Bvec2 = eval.fd(tfine, Bfd2)

Bhatvec1 = eval.fd(tfine, Bhatfd1$fd)
Bhatvec2 = eval.fd(tfine, Bhatfd2$fd)

# plot estimated (solid) and true (dashed) regression functions

par(ask=FALSE)
par(mfrow=c(2,1))
plot(tfine, Bvec1, type="l", lty=2)
lines(tfine, Bhatvec1, lty=1)
plot(tfine, Bvec2, type="l", lty=2)
lines(tfine, Bhatvec2, lty=1)

#  estimate the variance of the residuals and set up as diagonal matrix

yres = Yvec - Yhatvec
SigmaE = mean(yres^2)
SigmaE = SigmaE * diag(rep(1,N))

#  y2cMap is not needed for scalar dependent variable models

y2cMap = NULL

#  estimate the standard errors of the fitted values

fRegressResult = fRegress.stderr(fRegressResult, y2cMap, SigmaE)
YhatStderr = fRegressResult$YhatStderr

#  plot the estimated values, 95% confidence limits, true values, and data

par(ask=FALSE)
par(mfrow=c(1,1))
plot  (1:N, Yhatvec, type="b",   col=1, lwd=2)
points(1:N, Yhatvec+2*YhatStderr,col=2)
points(1:N, Yhatvec-2*YhatStderr,col=2)
points(1:N, Yvec, col=4, pch="*")
points(1:N, Yvectru, col=3, pch="x")
points(1:N, Yvec, col=4, pch="o")

#  now use predict.fRegress to get pointwise standard error of fit

predictResult = predict.fRegress(fRegressResult, xfdlist, se.fit = TRUE)

#  check that YhatStderr2 is the same as YhatStderr

YhatStderr2 = predictResult$YhatStderr
cbind(YhatStderr,YhatStderr2)



