### Approximate function
set.seed(29)
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

library(survival)
library(tensorflow)
library(keras)
#source("deepAFT.R")

a = 1.5
b = c(1.5, 0.8, 1.7, log(0.3),0.2, 1.2)

cen.time = 50
n = 500
max.iter = 50

alpha_lr = 0.01
decay_rate = 0.005
epochs.n = 100
batch.n = 20

#simulate the data
x1 = rbinom(n, 1, 0.5)
x2 = runif(n)*2 - 1
x3 = rnorm(n, 0, .25)
x4 = runif(n,0,2)
x  = cbind(x1, x2, x3,x4)
xa = cbind(x1, x2, x3^2, x2*x4, x2^x1, x4)
p = length(x[1,])
mu = a + xa%*%b
rm(xa)
etm = exp(mu)
ltm = rnorm(n, mu, 0.5)
tx = exp(ltm)

#hist(tx)
#plot(x, log(tx))
cx = rexp(n, 1/150)
status = ifelse(tx < cx, 1, 0)
time   = ifelse(status, tx, cx)
summary(time)
summary(status)
y = Surv(time, status)

# fit km curve for censoring
G_fit = survfit(Surv(time, 1-status)~1)
G = status/(.appxf(G_fit$surv, x=G_fit$time, xout = time) + 1e-10)

### Keras
model = keras_model_sequential()

### define model layers
model %>% layer_dense(units = 7, activation = 'selu', input_shape = c(p)) %>%
  layer_dense(units = 11, activation = 'relu') %>%
  layer_dropout(rate = 0.05)%>%
  layer_dense(units = 7, activation = 'relu') %>%
  layer_dense(units = 1)

### Compile (Define loss and optimizer)
model %>% compile(loss = 'mse', 
       optimizer = optimizer_adam(lr=alpha_lr, decay=decay_rate))
#summary(model)

epochs.n = 100
bach.n = 12
validation_split = 0.1
verbose = 1

lgt = log(time)
mean.ipt = mean(lgt)
lgt = lgt-mean.ipt
history = model%>%fit(x, lgt,
  epochs = epochs.n, batch_size = batch.n, sample_weight = G,
  validation_split = validation_split, verbose = verbose)

#linear predictors
lp = (model%>%predict(x)+mean.ipt)
pred_time = exp(lp)
plot(pred_time, time)
plot(log(time), log(pred_time))
