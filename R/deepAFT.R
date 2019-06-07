#### deep learning for AFT

deepAFT = function(x, ...) UseMethod("deepAFT")

deepAFT.formula = function(formula, model, epochs = 30, batch_size = 32, 
                    validation_split = 0.1, verbose = 0, epsilon = 0.01, max.iter = 50, 
                    data = list(...)) {
  mf = model.frame(formula=formula, data=data)
  
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  
  if (class(y) == "Surv") {
    family = "surv";
    st = sort(y[, 1], decreasing = TRUE, index.return = TRUE)
    idx = st$ix
    y = y[idx, ]
    x = x[idx, ]
  }
  do.call(deepAFT, x, y)
}
  
#deepAFT.default = function(x, y, model, epochs = 30, batch_size = 32, 
#  validation_split = 0.1, verbose = 0, epsilon = 0.01, max.iter = 50) {
deepAFT.default = function(x, y) {
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
  object = list(X = x, y = Surv(time, status), model = model, mean.ipt = mean.ipt, linear.predictors = lp,
                risk = exp(-lp), sfit = survfit(Surv(time, status)~1), residuals = resid, 
                iterations = k)
  class(object) = 'deepAFT'
  return(object)
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

################ testing code, no use for now

### replace function_fit with NN fit
#my_fit = lm(lgt~x)

### print out intermediate results if requested
# if(verbose > 0) {
#   # prediciton using linear model
#   ep2 = exp(predict(my_fit)+mean.ipt)
#   tmp = cbind(status, ep2, ep, ipt, ipt0, dat$time)
#   #outputs
#   colnames(tmp) = c('status', 'LinearModel E(T)', 'NN E(T)', 'Imputed T', '1st imputed T', 'censored T')
#   tmp = tmp[tmp[, 1] == 0, ]
#   print(tmp[1:15, -1])
#   cat("mean mu by linear model = ", mean(ep2), "\n")
#   cat("mean mu by NN           = ", mean(ep), "\n")
#   cat("mean imputed time       = ", mean(ipt), "\n")
# }
#############################
