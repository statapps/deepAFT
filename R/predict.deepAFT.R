### Approximate function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

predict.deepAFT = function(object, newdata, newy = NULL, ...) {
  sfit = survfit(object)
  lp   = object$predict
  risk = object$risk
  if(missing(newdata)) {
    x = object$x
    #martingale residual
    residuals = residuals(object, type = 'm')
  }
  else {
    ### if there is new data
    model = object$model
    x = newdata
    residuals = NULL

    lp  = (model%>%predict(x) + object$mean.ipt)
    risk = exp(-lp)
  }
  result = list(predictor = lp, risk = risk, time = object$time, sfit = sfit)

  if(!is.null(newy)) {
    if(missing(newdata)) stop("Error: newdata cannot missing when there is new y.")
    if(length(newy[, 1]) != length(x[, 1]))
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
