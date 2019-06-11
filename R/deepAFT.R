#### deep learning for AFT
deepAFT = function(x, ...) UseMethod("deepAFT")

deepAFT.formula = function(formula, model, data, control = list(...), 
        method = c("BuckleyJames", "ipcw", "transform"), ...) {
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

  #class(x) = switch(method, BuckleyJames="default", ipcw="ipcw", transform="transform")

  if (missing(control)) control = deepAFTcontrol(...)
    else control =  do.call("deepAFTcontrol", control)
  #fit = do.call("deepAFT", x, y, model, control)

  fit = switch(method, BuckleyJames=deepAFT.default(x, y, model, control), 
                       ipcw = deepAFT.ipcw(x, y, model, control),
                       transform = deepAFT.trans(x, y, model, control))
  return(fit)
}

deepAFT.default = function(x, y, model, control, ...) {
  epochs = control$epochs
  batch.n = control$batch.n
  v_split = control$v_split
  verbose = control$verbose
  max.iter= control$max.iter
  epsilon = control$epsilon
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
    
    history = model%>%fit(x, lgt, epochs = epochs, batch_size = batch.n,
      validation_split = v_split, verbose = verbose)

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
  if(!convergence) warning("Maximum iterations reached before converge!\n") 
  else cat('Algorithm converges after ', k, 'iterations!\n')
  
  ### create outputs
  object = list(X = x, y = Surv(time, status), model = model, mean.ipt = mean.ipt, 
    mu = lp, means = apply(x, 2, mean), 
    risk = exp(-lp), sfit = survfit(Surv(time, status)~1), residuals = resid, 
                iterations = k, method = "Buckley-James")
  class(object) = 'deepAFT'
  return(object)
}

deepAFT.ipcw = function(x, y, model, control, ...){
## epochs = 30, batch_size = 32, validation_split = 0.1, verbose = 1) {
  epochs = control$epochs
  batch.n = control$batch.n
  v_split = control$v_split
  verbose = control$verbose
  
  time = y[, 1]
  status = y[, 2]
  #n = length(status)
  #max.t = max(time)
  
  # fit km curve for censoring, set surv to small number if last obs fails.
  G_fit = survfit(Surv(time, 1-status)~1)
  G_fit$surv = ifelse(G_fit$surv > 0, G_fit$surv, 1e-5)
  G = status/.appxf(G_fit$surv, x=G_fit$time, xout = time)
  
  lgt = log(time)
  mean.ipt = mean(lgt)
  lgt = lgt-mean.ipt
  
  history = model%>%fit(x, lgt,
                        epochs = epochs, batch_size = batch.n, sample_weight = G,
                        validation_split = v_split, verbose = verbose)
  
  #linear predictors
  lp = (model%>%predict(x)+mean.ipt)
  pred_time = exp(lp)
  ### create outputs
  object = list(x = x, y = y, model = model, mean.ipt = mean.ipt, 
                mu = lp, means = apply(x, 2, mean),
                risk = exp(-lp), method = "ipcw")
  class(object) = 'deepAFT'
  return(object)
}

deepAFT.trans = function(x, y, model, control, ...){
  ## epochs = 30, batch_size = 32, validation_split = 0.1, verbose = 1) {
  epochs = control$epochs
  batch.n = control$batch.n
  v_split = control$v_split
  verbose = control$verbose
  
  time = y[, 1]
  status = y[, 2]
  #n = length(status)
  #max.t = max(time)
  
  # fit km curve for censoring, set surv to small number if last obs fails.
  # transformation based on book of Fan J (1996, page 168)
  Gfit = survfit(Surv(time, 1-status)~1)
  St = Gfit$surv
  tm = Gfit$time
  St = ifelse(St> 0, St, 1e-5)
  Gt = 1/.appxf(St, x=tm, xout = time)
  
  # Integrate from 0 to t of 1/G(u)
  dt = diff(c(0, tm))
  iG = cumsum(1/St*dt)
  iGt = .appxf(iG, x=tm, xout = time)

  a = min(((iGt - time)/(time*Gt-iGt))[status==1], na.rm=TRUE)
  phi2 = (1+a)*iGt
  phi1 = phi2 - a*time*Gt
  
  time2 = ifelse(status>0, phi1, phi2)
  
  lgt = log(time2)
  mean.ipt = mean(lgt)
  lgt = lgt-mean.ipt
  
  history = model%>%fit(x, lgt,
                        epochs = epochs, batch_size = batch.n, sample_weight = G,
                        validation_split = v_split, verbose = verbose)
  
  #linear predictors
  lp = (model%>%predict(x)+mean.ipt)
  pred_time = exp(lp)
  ### create outputs
  object = list(x = x, y = y, model = model, mean.ipt = mean.ipt, 
                mu = lp, means = apply(x, 2, mean),
                risk = exp(-lp), method = "transform")
  class(object) = 'deepAFT'
  return(object)
}

deepAFTcontrol = function(epochs = 30, batch.n = 32,
  v_split = 0.1, verbose = 0, epsilon = 0.01, max.iter = 50) {
  
  list(epochs = epochs, batch.n = batch.n, v_split = v_split, verbose = verbose, epsilon = epsilon, max.iter = max.iter)
}

plot.deepAFT = function(x, type = c('predicted', 'residuals', 'baselineKM'), ...) {
  type = match.arg(type)
  time = x$y[, 1]
  log.time  = log(time)
  if (type == 'predicted') {
    predicted = x$mu
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
  cat("Deep AFT model with ", x$method, 'method\n\n')
  risk = as.vector(x$risk)
  y = x$y
  location = (1/as.vector(risk))
  cat("Summary of predict location score exp(mu):\n")
  print(summary(location))
  
  cat("\n\nDistribution of T0 = T/exp(mu):\n")
  
  sfit = survfit(x)
  print(sfit)
  
  cindex = survConcordance(y~risk)
  cat("Concordance index: ", cindex$concordance, "\n")
  #print(cindex)
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
