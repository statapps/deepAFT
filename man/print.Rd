\name{print}
\alias{print.deepAFT}
\alias{print.summary.deepAFT}
\alias{summary.deepAFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print a summary of fitted deepAFT object
}
\description{
print is used to provide a short summary of outputs from "deepAFT".
}
\usage{
\method{print}{deepAFT}(x, ...)
\method{print}{summary.deepAFT}(x, ...)
\method{summary}{deepAFT}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class returned from deepAFT model fit}
  \item{\dots}{other options used in print()}
  \item{object}{a class of deepAFT object}
}
\details{
  print.deepAFT is called to print object or summary of object from the deep learning AFT models \code{\link{deepAFT}}.  
  summary(fit) provides detail summary of `deepAFT' model fit, including predictors, 
baseline survival function for T0=T/exp(mu), and martingale residules for the fitted model.
  The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}

\value{
  %%  ~Describe the value returned
  An object of class "summary.deepAFT" is returned. The object contains the following list components:

 %%  If it is a LIST, use
   \item{location}{location parameter exp(mu), to predice the mean value of survival time.}
   \item{sfit}{survfit object of the baselie survival function of T0=T/exp(mu).}
   \item{cindex}{Concordance index of the fitted deepAFT model.}
   \item{resid}{martingle residuals of the fitted deepAFT model.}
   \item{method}{the model used to fit the deepAFT model.}
}

\author{Bingshu E. Chen}
\seealso{
 The default method for print \code{\link{print.default}}. Other methods include
\code{\link{survreg}}, 
\code{\link{deepAFT}}, 
\code{\link{summary}}
}
\examples{
#
#  print(fit)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ summary }
