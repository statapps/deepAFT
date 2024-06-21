\name{deepAFT-package}
\alias{deepAFT-package}
\alias{deepAFT-doc}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deep learning for accelerated failure time data
}
\description{
A deep learning package for survival analysis using accelerated failure time (AFT) models.
}

\details{
"deepAFT" is a R package for deep learning accerelated failure time Models.
Please use the following steps to install the most recent version of 'deepAFT' package:

1. First, you need to install the 'devtools' package. 
   You can skip this step if you have 'devtools' installed in your R. Invoke R and then type

   install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "deepAFT" package from github with R commond

  install_github("statapps/deepAFT")


"deepAFT" uses deep learning methods to fir an accerelated failure time model for survival data with high dimensional  
biomarker variables or other covariates.

A stable version of the "deepAFT" is now part of an R package "dnn" available from the Comprehensive R Archive Network 
(https://CRAN.R-project.org/package=dnn) and can be installed using R command 

  install.packages("dnn")
  library(dnn)
  help("deepAFT")

}
\author{Patrick Norman, Wanlux, Li, Wenyu Jiang and Bingshu E. Chen}

Maintainer: Bingshu E. Chen <bingshu.chen@queensu.ca>}
\references{
Buckley, J. and James, I. (1979). Linear regression with cencored data. Biometrika, 66, page 429-436.

Norman, P. and Chen, B. E. (2018). DeepAFAT: A nonparametric accelerated failure time model with artifical neural network. M. Sc. thesis report. Department of Public Health Sciences, Queen's University, Canada.

Norman, P. Li, W., Jiang, W. and Chen, B. E. (2024). deepAFT: A nonlinear accelerated failure time model with artificial neural network. Statistics in Medicine. Published June 18, 2024.

Chollet, F. and Allaire J. J. (2017). Deep learning with R. Manning.
}

\note{
For now, deepAFT works only for right censored data. Methods for other types of censoring will be developed soon.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{coxph, }
\code{glm, }
\code{survival}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
# fit = deepAFT(Surv(time, status) ~ x)
# print(fit)
#  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{deepAFT}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

