\name{deepAFT}
\alias{deepAFT}
\alias{bhm.default}
\alias{bhm.formula}
\alias{bhmFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deep learning for accelerated failure time data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deepAFT(x, y, model, epochs = 30, batch_size = 32, 
  validation_split = 0.1, verbose = 0, epsilon = 0.01, max.iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Covariates for AFT model}
  \item{y}{Surv object for AFT model}
  \item{model}{Keras model, see below for details.}
  \item{epochs}{}
  \item{batch_size}{}
  \item{validation_split}{}
  \item{verbose}{}
  \item{epsilon}{}
  \item{max.iter}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Buckley, J. and James, I. (1979). Linear regression with cencored data. Biometrika, 66, page 429-436.
Norman, P. and Chen, B. E. (2019). DeepAFT: A nonparametric accelerated failure time model with artifical neural network. Manuscript to be submitted. 
}
\author{
%%  ~~who you are~~
Chen, B. E. and Norman P.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
