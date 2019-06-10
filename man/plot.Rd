\name{plot}
\alias{plot.deepAFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a fitted deep learning AFT model
}
\description{
Several different type of plots can be produced for deep AFT mdels. 
Plot method is used to provide a graphic plot  from "deepAFT" model fit.

Use "methods(plot)" and the documentation for these for other plot methods.
}
\usage{
\method{plot}{deepAFT}(x, type = c("predicted", "residuals", "baselineKM"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class returned from "deepAFT" fit.}
  \item{type}{type of plot in deepAFT object, "predicted" to plot the linear predicted values, "residuals" to plot residuals, "baselineKM" to plot baseline Kaplan-Meier survival curve.}
  \item{\dots}{other options used in plot().}
}
\details{
  plot.deepAFT is called to plot the fitted deep learning AFT model.

  The default method, plot.default has its own help page. Use methods("plot") to get all the methods for the plot generic.
}
\author{Bingshu E. Chen}
\seealso{
 The default method for plot \code{\link{plot.default}}.
\code{\link{deepAFT}}
}
\examples{
#
#  plot(fit)
# 
######## plot for deepAFT object
#
#  plot(fit, type = 'baselineKM')
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
