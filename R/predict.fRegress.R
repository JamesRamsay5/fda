predict.fRegress <- function(object, newdata=NULL, se.fit = FALSE,
                             interval = c("none", "confidence", "prediction"),
                             level = 0.95, ...){
  #
  #  Argument object is a list object with the class name "fRegress" and is produced by
  #  a call to function fRegress.
  ##
  ## 1.  fit ???
  ##
  if(is.null(newdata))
    pred <- object$yhatfdobj
  else{
    nx <- length(object$xfdlist)
    Nnew <- length(newdata)
    pred <- rep(0, Nnew)
    for(i in 1:nx){
      xi <- predict(object$xfdlist[[i]], newdata)
      bi <- predict(object$betaestlist[[i]], newdata)
      pred <- pred+bi*xi
    }
  }
  ##
  ## 2.  Need se.fit?
  ##
  int <- match.arg(interval)
  need.se <- (se.fit || (int != "none"))
  if(!need.se)return(pred)
  #
  stop('Need se.fit;  not implemented yet')
}

residuals.fRegress <- function(object, ...){
  object$yfdPar - predict(object, ...)
}
