\name{ibs}
\alias{ibs}
\alias{ibs.default}
\alias{ibs.deepAFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calcuate integrated Brier Score for deepAFT 
}
\description{
The function ibs is used to calcualte integrated Brier Score for deepAFT.
}
\usage{
\method{ibs}{deepAFT}(object, newdata=NULL, newy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the results of a deepAFT fit}
  \item{newdata}{If no null, new data and new y will be used for calculation}
  \item{newy}{Used together with new data}
}
\details{
  ibs.deepAFT is called to calculate integrate Brier Score for the deepAFT model \code{\link{deepAFT}}.

  The default method, ibs.default has its own help page. Use methods("ibs") to get all the methods for the ibs generic.
}
\author{Bingshu E. Chen}
\seealso{
 The default method for ibs \code{\link{ibs.default}},
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
