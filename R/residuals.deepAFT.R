### baseline cumulative hazard function and martingale residuals for deepAFT
survfit.deepAFT = function(object, se.fit=TRUE, conf.int = .95, ...) {
  #baseline survival function S_T0(t), with all covariates value = 0
  #where T0 = T/exp(mu), or log(T0) = log(T) - mu
  y0 = object$y
  y0[, 1] = y0[, 1]/object$risk
  sfit = survfit(y0 ~ 1, se.fit=se.fit)
  return(sfit)
}

### Residuals of deepAFT
residuals.deepAFT = function(object, type=c("martingale", "deviance"), ...) {
  type = match.arg(type)
  sfit = survfit(object, se.fit = FALSE)
  rr = sfit$residuals
  status = sfit$n.event
  drr = sign(rr)*sqrt(-2*(rr+ifelse(status==0, 0, status*log(status-rr))))
  resid = switch(type, martingale = rr, deviance = drr)
  return(resid)
}
