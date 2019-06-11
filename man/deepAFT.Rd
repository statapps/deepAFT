\name{deepAFT}
\alias{deepAFT}
\alias{deepAFT.default}
\alias{deepAFT.formula}
\alias{deepAFT.ipcw}
\alias{deepAFTcontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deep learning for accelerated failure time data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fit a deep learning survival regression model. 
These are location-scale models for an arbitrary transform of the time variable; 
the most common cases use a log transformation, leading to accelerated failure time models.
}
\usage{
deepAFT(x, ...)

\method{deepAFT}{formula}(formula, model, data, control = list(...), 
              method = c("BuckleyJames", "ipcw", "transform"), ...)

\method{deepAFT}{default}(x, y, model, control, ...)

\method{deepAFT}{ipcw}(x, y, model, control, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula expression as for other regression models. 
    The response is usually a survival object as returned by the
          'Surv' function.  See the documentation fo'Surv', 'lm' and
          'formula' for details.}
  \item{model}{Keras model, see below for details.}
  \item{data}{a data.frame in which to interpret the variables named in the formula.}
  \item{x}{Covariates for AFT model}
  \item{y}{Surv object for AFT model}
  \item{method}{methods to handl censoring data in deep AFT model fit, 'BuckleyJames' for Buckley and James method, 
          'ipcw' for inverse probability censoring weights method.
          'transform' for transformation based on book of Fan and Gijbels (1996, page 168)}
  \item{control}{a list of control values, in the format produced by
          'deepAFTcontrol'. The default value i'deepAFTcontrol()'}
  \item{\dots}{optional arguments}
%%  \item{epochs}{}
%%  \item{batch_size}{}
%%  \item{validation_split}{}
%%  \item{verbose}{}
%%  \item{epsilon}{}
%%  \item{max.iter}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
See "Deep learning with R" for details on how to build a deep learning model. 

The following parameters in 'deepAFTcontrol' will be used to control the model fit process.

'epochs': number of deep learning epochs, default is 30.

'batch.n': batch size, default is 128. 'NaN' may be generated if batch size is too small and there is not event in a batch. 

'v_split': validation split, a numerical value in [0, 1), for portion of samples used for validation, default is 0.

'verbose': vervose = 1 for print out verbose during the model fit, 0 for not print.

'epsilon': espilon for converence check, default is epsilon = 0.001.

'max.iter': number of maximum iteration, default is max.iter = 30.

}
\value{
%%  ~Describe the value returned
An object of class "deepAFT" is returned. The deepAFT object contains the following list components:

%%  If it is a LIST, use
  \item{x}{Covariates for AFT model}
  \item{y}{Surv object for AFT model}
  \item{model}{Keras model}
  \item{mean.ipt}{mean survival or censoring time}
  \item{predictor}{predictor score mu = f(x)} 
  \item{risk}{risk score = exp(predictor)}
  \item{method}{method for deepAFT fitting, either Buckley-James or IPCW}
}
\references{
%% ~put references to the literature/web site here ~
Buckley, J. and James, I. (1979). Linear regression with cencored data. Biometrika, 66, page 429-436.

Norman, P. and Chen, B. E. (2019). DeepAFAT: A nonparametric accelerated failure time model with artifical neural network. Manuscript to be submitted. 

Chollet, F. and Allaire J. J. (2017). Deep learning with R. Manning.
}
\author{
Chen, B. E. and Norman P.
}
\note{For right censored survival time only}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{print.deepAFT}}, \code{\link{survreg}}, \code{\link{ibs}}
}
\examples{
## Example for deep learning model
##
### Keras

#model = keras_model_sequential()
### define model layers

#model %>% layer_dense(units = 7, activation = 'relu', input_shape = c(p)) %>%
#  layer_dense(units = 11, activation = 'relu') %>%
#  layer_dropout(rate = 0.1)%>%
#  layer_dense(units = 7, activation = 'relu') %>%
#  layer_dense(units = 1)

### Compile (Define loss and optimizer)
#model %>% compile(loss = 'mse', 
#       optimizer = optimizer_adam(lr=alpha_lr, decay=decay_rate))
#summary(model)

# fit = deepAFT(Surv(time, status) ~ x) 
# print(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{keras}% use one of  RShowDoc("KEYWORDS")
\keyword{deepAFT}% __ONLY ONE__ keyword per line
