### Approximate function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

predict.deepAFT = function(object, newdata, newy = NULL, ...) {
  result = summary(object)
  sfit = result$sfit
  if(missing(newdata)) {
    return(result)
  }
  else {
    ### if there is new data
    m = object$model
    x = newdata
    ### if x is a numeric vector, change it to matrix
    if(is.null(dim(x))) x = t(as.matrix(x))

    lp  = (m%>%predict(x) + object$mean.ipt)
    risk = exp(-lp)
    result$predictors = lp
    result$risk = risk
    result$locations = 1/risk
    result$cindex = NULL
    result$residuals = NULL
  }

  if(!is.null(newy)) {
    if(missing(newdata)) stop("Error: newdata cannot missing when there is new y.")
    if(length(newy[, 1]) != length(x[, 1]))
      stop("Error: new y shall have the same subjects as the new data.")

    time   = newy[, 1]  #time
    status = newy[, 2]  #status

    #baseline survival function
    aft.time = risk*time
    sf = .appxf(sfit$surv, x=sfit$time, xout=aft.time)
    sf = ifelse(sf>0, sf, min(sf[sf>0]))
    cumhaz = -log(sf)
    result$residuals = (status - cumhaz)
    if(exists("survConcordance"))
      cindex = survConcordance(newy~risk)
    else 
      cindex = concordance(newy~risk)
    result$cindex = cindex
    result$c.index= cindex$concordance
    class(result) = "summary.deepAFT"
  }
  return(result)
}
