\name{residuals}
\alias{residuals.deepAFT}
\title{
Calculate Residuals for a deepAFT Fit. 
}
\description{
Calculates martingle, deviance or Cox-Snell residuals for a previously fitted (deepAFT) model.
}
\usage{
## S3 method for class 'deepAFT'
\method{residuals}{deepAFT}(object, type = c("martingale", "deviance", "CoxSnell"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the results of a (deepAFT) fit.}
  \item{type}{character string indicating the type of residual desired. Possible values are "martingale", "deviance". Only enough of the string to determine a unique match is required.}
  \item{\dots}{other unused arguments.}
}
\details{
  residuals.deepAFT is called to compuate baseline survival function S_T0(t) from the deepAFT model \code{\link{deepAFT}}, where T0 = T/exp(mu), or log(T) = log(T) - mu.

  The default method, residuals has its own help page. Use methods("residuals") to get all the methods for the residuals generic.
}

\value{
For martingale and deviance residuals, the returned object is a vector with one element for each subject.
The row order will match the input data for the original fit.

See \code{\link{residuals}} for more detail about other output values. 
}

\note{For deviance residuals, the status variable may need to be reconstructed.}

\author{Bingshu E. Chen}

\seealso{
 The default method for residuals \code{\link{residuals}},
\code{\link{predict.deepAFT}}, \code{\link{survfit.deepAFT}}, and \code{\link{deepAFT}}.
}
\examples{
#
#  residuals(fit, type = "martingale")
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ residuals }
