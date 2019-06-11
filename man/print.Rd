\name{print}
\alias{print.deepAFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print a summary of fitted deepAFT object
}
\description{
print is used to provide a short summary of outputs from "deepAFT".
}
\usage{
\method{print}{deepAFT}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class returned from deepAFT model fit}
  \item{\dots}{other options used in print()}
}
\details{
  print.deepAFT is called to print object or summary of object from the deep learning AFT models \code{\link{deepAFT}}.  
  summary(fit) provides detail summary of `deepAFT' model fit, including parameter estimates, 
standard errors, and 95 percent CIs.
  The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}
\author{Bingshu E. Chen}
\seealso{
 The default method for print \code{\link{print.default}}. Other methods include
\code{\link{survreg}}, 
\code{\link{deepAFT}}, 
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
