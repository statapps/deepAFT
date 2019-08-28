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
    result$c.index = NULL
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
    
    cindex = try(concordance(newy~lp))
    if(class(cindex)=="try-error") then cindex = NULL
    else {
      result$cindex = cindex
      result$c.index= cindex$concordance
    }
    class(result) = "summary.deepAFT"
  }
  return(result)
}

### Find IPCW for factor x (with 5 or less levels).
.ipcw = function(time, status, x=NULL) {
  n = length(status)
  event = 1 - status
  if(is.null(x)) x = rep(1, n)

  # Fit a Cox model is input 'x' is a matrix;
  if(!is.vector(x)) {
    cfit = coxph(Surv(time, event)~x)
    gt = exp(-predict(cfit, type = 'expected'))
    smin = min(gt[gt>0]) #set surv to a small number if the last obs fails.
    gt = ifelse(gt>0, gt, smin)
    return(gt)
  }
  # Fit a nonparametric model if input x is a factor or a vector 
  if(length(unique(x))>20) xf = cut(x, 5)
  else xf = as.factor(x)
  xn = as.numeric(xf)
  max.xn = max(xn)
  gt = rep(NaN, n)
  if(max.xn > n/20) stop("Too many censoring groups. Please use 5 or less censoring groups.")
  for(i in 1:max.xn) {
    idx = (xn == i)
    gf = survfit(Surv(time, event)~xn, subset = idx)
    ti = time[idx]
    tg = gf$time
    sg = gf$surv
    smin = min(sg[sg>0]) #set surv to a small number if the last obs fails.
    sg = ifelse(sg > 0, sg, smin) 
    si = .appxf(sg, tg, ti)
    gt[idx] = si
  }
  return(gt)
}

