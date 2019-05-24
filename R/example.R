#### recursive KM for AFT
rm(list = ls())
library(survival)
library(keras)
source("./deepAFT/R/deepAFT.R")

a = 1.5
b = c(1.5, 0.8, 1.7, log(0.3),0.2, 1.2)

cen.time = 50
n = 500
max.iter = 50

alpha_lr = 0.01
decay_rate = 0.0005
epochs.n = 1000
batch.n = 200

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


### Keras
model = keras_model_sequential()

### define model layers
model %>% layer_dense(units = 7, activation = 'selu', input_shape = c(p)) %>%
  layer_dense(units = 11, activation = 'selu') %>%
  # layer_dropout(rate = 0.1)%>%
  layer_dense(units = 7, activation = 'selu') %>%
  layer_dense(units = 1)

### Compile (Define loss and optimizer)
model %>% compile(loss = 'mse', 
       optimizer = optimizer_adam(lr=alpha_lr, decay=decay_rate))
#summary(model)

#standardize X
#.stndx = function(a) {return((a - mean(a))/sd(a))}
#x = apply(x, 2, .stndx)

y = Surv(time, status)

object = deepAFT(x, y, model, epochs = epochs.n, batch_size = batch.n, 
    validation_split = 0, verbose = 0, max.iter = max.iter, epsilon = 0.05)

cat('censoring = ', (1-mean(status))*100, '%\n')
plot(object, type = 'resid')
plot(object)

### C-index calculation
mu<-object$linear.predictors
id=seq(1,length(time))
c_data<-data.frame(id,time,mu,status)
#c_data<-c_data[which(c_data$status!=0),]
concordantpairs<-0
totalpairs=0
for(i in 1:(length(c_data$time)-1)){
  for(j in (i+1):length(c_data$id)){
    if(c_data$time[i]!=c_data$time[j]){
      totalpairs<- totalpairs+1
      # if(c_data$time[i]==c_data$time[j]){print("YOU HAVE A TIE")} #####How do we account for tied failure times (tied predictions are fine)?
      if(c_data$time[i]<c_data$time[j]){
        if (c_data$mu[i]<c_data$time[j]){
          concordantpairs<<- concordantpairs+1
        }
      }
      else if (c_data$time[i]>c_data$time[j]){
        if(c_data$mu[i]>c_data$mu[j]){
          concordantpairs<<-concordantpairs+1
        }
      }
      else if(c_data$mu[i]==c_data$mu[j]){
        concordance<<-concordance +0.5
      }
    }
  }
}
concordantpairs
totalpairs
concordance=concordantpairs/totalpairs
concordance
plot(log(c_data$time),c_data$mu)
