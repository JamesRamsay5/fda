fRegress.fd <- function(y, xfdlist, betalist, wt=NULL,
                        y2cMap=NULL, SigmaE=NULL, returnMatrix=FALSE, ...){
  yfdPar <- fdPar(y, ...)
  print("inside fRegress.fd")
  fRegress(yfdPar, xfdlist, betalist, wt=wt,
                 y2cMap=y2cMap, SigmaE=SigmaE, returnMatrix=FALSE, ...)
}
