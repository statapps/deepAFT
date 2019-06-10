#### deep learning for AFT

deepAFT = function(x, ...) UseMethod("deepAFT")

deepAFT.formula = function(formula, model, data, control = list(...), 
        method = c("BuckleyJames", "ipcw", "transform")) {
  if (missing(data)) 
    data = environment(formula)

  mf = model.frame(formula=formula, data=data)
  method = match.arg(method)

  x = model.matrix(attr(mf, "terms"), data = mf)
  ### remove intercept term
  x = x[, -1]
  y = model.response(mf)
  if (!inherits(y, "Surv")) 
    stop("Response must be a survival object")

  type = attr(y, "type")
  if (type == "counting") 
    stop("start-stop type Surv objects are not supported")
  if (type == "mright" || type == "mcounting") 
    stop("multi-state survival is not supported")

  class(x) = switch(method, BuckleyJames="default", ipcw="ipcw", transform="transform")

  if (missing(control)) control = survreg.control(...)
    else control =  do.call("survreg.control", control)

  fit = do.call(deepAFT, x, y, model, control, ...)
  return(fit)
}
  
deepAFT.default = function(x, y, model, epochs = 30, batch_size = 32, 
  validation_split = 0.1, verbose = 0, epsilon = 0.01, max.iter = 50) {
  time = y[, 1]
  status = y[, 2]
  n = length(status)
  max.t = max(time)
  id = 1:n
  dat = data.frame(cbind(id = id, time = time, status = status))
  ep = 1
  dati = dat
  ipt = .imputeKM(dat)*ep
  ipt0 = ipt
  mean.ipt = mean(log(ipt))
  
  convergence = FALSE
  for(k in 1:max.iter) {
    ###lgy = log(T), with T = imputed time (ipt)
    lgt = log(ipt) - mean.ipt
    
    history = model%>%fit(x, lgt,
      epochs = epochs.n, batch_size = batch.n,
      validation_split = validation_split, verbose = verbose)

    ep0 = ep
    #linear predictors
    lp = (model%>%predict(x)+mean.ipt)
    ep = exp(lp)

    ### do imputation for censoring time
    et = dat$time/ep
    ### restrict rescaled time to be less than max.t
    et = ifelse(et < max.t, et, max.t)
    dati = data.frame(cbind(id = id, time = et, status = dat$status))
    colnames(dati) = c("id", "time", "status")
    ipt = .imputeKM(dati)*ep
    ### restrict imputed time to be less than max.t
    ipt = ifelse(ipt < max.t, ipt, max.t)
    resid = (log(ipt) - lp)
    cat('MSE = ', mean(resid^2))

    #check convergence
    dif.ep = mean(abs(ep-ep0))
    cat(",  epsilon = ", dif.ep, "\n")
    if(dif.ep < epsilon) {
      convergence = TRUE
      break
    }
  }
  if(!convergence) warning("Maximum iterations reached before converge!") 
  else cat('Algorithm converges after ', k, 'iterations!')
  
  ### create outputs
  object = list(X = x, y = Surv(time, status), model = model, mean.ipt = mean.ipt, 
    linear.predictors = lp, means = apply(x, 2, mean), 
                risk = exp(-lp), sfit = survfit(Surv(time, status)~1), residuals = resid, 
                iterations = k, method = "Buckley-James")
  class(object) = 'deepAFT'
  return(object)
}

deepAFT.control = function(epochs = 30, batch_size = 32,
  validation_split = 0.1, verbose = 0, epsilon = 0.01, iter.max = 50) {
  
  list(epochs = epochs, batch_size = batch_size, validation_split = validation_split, verbose = verbose, epsilon = epsilon, iter.max = iter.max)
}

plot.deepAFT = function(x, type = c('predicted', 'residuals', 'baselineKM'), ...) {
  type = match.arg(type)
  time = x$y[, 1]
  log.time  = log(time)
  if (type == 'predicted') {
    predicted = x$linear.predictors
    plot(log.time, predicted, xlab = 'Log survival time', ylab = 'Predicted log survival time')
    abline(0, 1, lty = 2)
  } else if(type == 'residuals') {
    resid = x$residuals
    plot(log.time, resid, xlab = 'Log survival time', ylab = 'Residuals of linear predictors')
    abline(0, 0, lty = 2)
  } else if(type == 'baselineKM') {
    plot(x$sfit)
    title('Baseline KM curve for X = mean(X)')
  }
}

print.deepAFT = function(x, ...) {
  sfit = survfit(x)
  print(sfit)
}

#### impute KM for AFT
### impute censoring time
## st is a n x 2 matrix, with st[ ,1] as time and st[, 2] as survival function
.imputeFun = function(tc, st) {
  sc = st[st[, 1] > tc, ]
  if(is.vector(sc)) return(sc[1])
  
  ## P(T>tc)
  sm = sum(sc[, 2])
  
  ##conditional probability mass function
  pmf  = sc[, 2]/sm
  ## imputed survival time
  ipt = sum(sc[, 1]*pmf)
  return(ipt)
}

.imputeKM = function(dat) {
  sf = survfit(Surv(time, status)~1, data = dat)
  sv = sf$surv
  #sv[length(sv)] = 0
  st = cbind(sf$time, -diff(c(1, sv)))

  idc = dat[dat$status == 0, 1]
  ipt = dat$time
  d.max = max(dat$time[dat$status > 0])
  c.max = max(dat$time[dat$status < 1])
  for (i in idc) {
    tc = dat$time[i]
    if (tc < d.max)
      ipt[i] = .imputeFun(tc, st)
    else ipt[i] = c.max
  }
  return(ipt)
}

deepAFT.ipcw = function(x, y, model, epochs = 30, batch_size = 32, 
  validation_split = 0.1, verbose = 1) {

  .appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }
  time = y[, 1]
  status = y[, 2]
  n = length(status)
  max.t = max(time)

  # fit km curve for censoring
  G_fit = survfit(Surv(time, 1-status)~1)
  G = status/(.appxf(G_fit$surv, x=G_fit$time, xout = time) + 1e-10)

  lgt = log(time)
  mean.ipt = mean(lgt)
  lgt = lgt-mean.ipt

  history = model%>%fit(x, lgt,
    epochs = epochs.n, batch_size = batch.n, sample_weight = G,
    validation_split = validation_split, verbose = verbose)

  #linear predictors
  lp = (model%>%predict(x)+mean.ipt)
  pred_time = exp(lp)
  ### create outputs
  object = list(x = x, y = y, model = model, mean.ipt = mean.ipt, 
    linear.predictors = lp, means = apply(x, 2, mean),
    risk = exp(-lp), method = "ipcw")
  class(object) = 'deepAFT'
  return(object)
}
