\name{deepAFT-package}
\alias{deepAFT-package}
\alias{deepAFT-doc}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deep learning for accelerated failure time data
}
\description{
A deep learning package for survival analysis using accelerated failure time (AFT) model.
}

\details{
"deepAFT" is a R package for Biomarker Threshold Models.
Please use the following steps to install the most recent version of 'deepAFT' package:

1. First, you need to install the 'devtools' package. 
   You can skip this step if you have 'devtools' installed in your R. Invoke R and then type

   install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "deepAFT" package from github with R commond

  install_github("statapps/deepAFT")


"deepAFT" uses different statistical methods to identify cut-point (thershold parameter) for the 
biomarker in either generalized linear models or Cox proportional hazards model.

A stable version of View the "deepAFT" package is also available from the Comprehensive R Archive Network 
(https://CRAN.R-project.org/package=deepAFT) and can be installed using R command 

  install.packages("deepAFT")

}
\author{Bingshu E. Chen and Patrick Norman

Maintainer: Bingshu E. Chen <bingshu.chen@queensu.ca>}
\references{
 Put references to the literature/web site here ~
}

\note{
%%  ~~further notes~~
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

