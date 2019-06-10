\name{predict}
\alias{predict.deepAFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicted Values for a deepAFT Object
}
\description{
Compute predicted values for a deepAFT object
}
\usage{
## S3 method for class 'deepAFT'
\method{predict}{deepAFT}(object, newdata, newy=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the results of a model fit using the deepAFT function}
  \item{newdata}{optional new data at which to do predictions. If absent, predictions are for the dataframe used in the original fit}
  \item{newy}{optional new outcome variable y}
  \item{\dots}{other options used in predict()}
}
\details{
  predict.deepAFT is called to predict object from the deepAFT model \code{\link{deepAFT}}.

  The default method, predict has its own help page. Use methods("predict") to get all the methods for the predict generic.
}

\value{
%%  ~Describe the value returned
predict.deepAFT returns a list of predicted values, prediction error and residuals. 

  \item{lp}{linear predictor of beta(w)*Z, where beta(w) is the fitted regression coefficient and Z is covariance matrix.}
  \item{risk}{risk score, exp(lp). When new y is provided, both lp and risk will be ordered by survival time of the new y.}
  \item{cumhaz}{cumulative hzard function.}
  \item{time}{time for cumulative hazard function. Time from new y will be used is provided}
}
\author{Bingshu E. Chen}
\seealso{
 The default method for predict \code{\link{predict}},
\code{\link{deepAFT}}, \code{\link{survfit.deepAFT}}
}

\examples{
#
#  predict(fit)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ predict }

