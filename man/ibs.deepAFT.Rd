\name{ibs}
\alias{ibs}
\alias{ibs.deepAFT}
\alias{ibs.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calcuate integrated Brier Score for deepAFT 
}
\description{
The function ibs is used to calcualte integrated Brier Score for deepAFT.
}
\usage{
ibs(object, ...)
\method{ibs}{default}(object, ...)
\method{ibs}{deepAFT}(object, newdata=NULL, newy = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the results of a deepAFT fit.}
  \item{newdata}{optional argument, if no null, new data and new y will be used for calculation.}
  \item{newy}{optional argument, used together with new data.}
  \item{\dots}{other unused arguments.}
}
\details{
  ibs is called to calculate integrate Brier Score for the deepAFT model \code{\link{deepAFT}}.
}
\author{Bingshu E. Chen}
\seealso{
  \code{\link{deepAFT}}
}
\examples{
### To calculate Brier score for the original fitted data
#
#  ibs(fit)
#
### To calculate Brier score for new data with new outcomes
#
#  ibs(fit, newdata, newy)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Integrated Brier Score }
