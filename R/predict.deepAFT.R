### Approximate function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

predict.deepAFT = function(object, newdata, newy = NULL) {
  sfit = object$sfit
  lp   = object$linear.predictors
  risk = object$risk
  if(missing(newdata)) {
    X = object$X
    #residuals = sfit$residuals
    residuals = NULL
  }
  else {
    ### if there is new data
    model = object$model
    X = newdata
    residuals = NULL

    lp  = (model%>%predict(X) + object$mean.ipt)
    risk = exp(-lp)
  }
  result = list(lp = lp, risk = risk, time = object$time, sfit = sfit)

  if(!is.null(newy)) {
    if(missing(newdata)) stop("Error: newdata cannot missing when there is new y.")
    if(length(newy[, 1]) != length(X[, 1]))
      stop("Error: new y shall have the same subjects as the new data.")

    ## sort survival time for new y from smallest to largest
    idx = order(newy[, 1])
    lp  = lp[idx]
    risk = risk[idx]
    newy = newy[idx, ]

    time   = newy[, 1]  #time
    status = newy[, 2]  #status

    #baseline survival function
    aft.time = risk*time
    sf = .appxf(sfit$surv, x=sfit$time, xout=aft.time)
    sf = ifelse(sf>0, sf, min(sf[sf>0]))
    cumhaz = -log(sf)
    residuals = (status - cumhaz)

    ## update lp, riks, time and cumhaz, all ordered by time
    result$lp = lp
    result$risk = risk
    result$time = time
    result$status = status
    result$surv = sf
  }
  result$residuals = residuals

  return(result)
}