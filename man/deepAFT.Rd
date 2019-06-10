\name{deepAFT}
\alias{deepAFT}
\alias{deepAFT.default}
\alias{deepAFT.formula}
\alias{deepAFT.ipcw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deep learning for accelerated failure time data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fit a deep learning survival regression model. These are
location-scale models for an arbitrary transform of the time
variable; the most common cases use a log transformation, leading
to accelerated failure time models.
}
\usage{
deepAFT(x, y, model, epochs = 30, batch_size = 32, 
  validation_split = 0.1, verbose = 0, epsilon = 0.01, max.iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula expression as for other regression models.  The
          response is usually a survival object as returned by the
          ‘Surv’ function.  See the documentation for ‘Surv’, ‘lm’ and
          ‘formula’ for details.}
  \item{x}{Covariates for AFT model}
  \item{y}{Surv object for AFT model}
  \item{model}{Keras model, see below for details.}
  \item{control}{a list of control values, in the format produced by
          ‘deepAFT.control’. The default value is ‘deepAFT.control()'}
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
}
\value{
%%  ~Describe the value returned
An object of class "deepAFT" is returned.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{survreg}}, \code{\link{ibs.deepAFT}}
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
