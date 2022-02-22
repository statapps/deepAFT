### baseline cumulative hazard function and martingale residuals for deepAFT
survfit.deepAFT = function(formula, se.fit=TRUE, conf.int = .95, ...) {
  #baseline survival function S_T0(t), with all covariates value = 0
  #where T0 = T/exp(mu), or log(T0) = log(T) - mu, where risk = exp(-mu)
  y0 = formula$y
  y0[, 1] = y0[, 1]*formula$risk
  sfit = survfit(y0 ~ 1, se.fit=se.fit, conf.int=conf.int)
  return(sfit)
}

### Residuals of deepAFT
residuals.deepAFT = function(object, type=c("martingale", "deviance", "coxSnell"), ...) {
  type = match.arg(type)
  sfit = survfit(object, se.fit = FALSE)
  
  time   = object$y[, 1]*object$risk
  status = object$y[, 2]
  
  m = length(sfit$surv)
  ### in case the last subject fails,  S0(t) = 0
  sfit$surv[m] = sfit$surv[m-1]

  # fit survival function at time Ti
  St = .appxf(sfit$surv, x=sfit$time, xout = time)

  # Cox-Snell residual H(T)
  Ht = -log(St)
  
  rr = status - Ht
  drr = sign(rr)*sqrt(-2*(rr+ifelse(status==0, 0, status*log(status-rr))))
  resid = switch(type, martingale = rr, deviance = drr, coxSnell=Ht)
  return(resid)
}
