\name{WalkingStickFemurs}
\alias{WalkingStickFemurs}
\docType{data}

\title{
Walking Stick Femur Length
}

\description{
Two measures of femur length for each of 25 walking sticks (\emph{Timema cristinae}).
}

\usage{data(WalkingStickFemurs)}
\format{
  A data frame with 50 observations on the following 2 variables.
  \describe{
    \item{\code{specimen}}{a factor with levels \code{1-25}}
    \item{\code{femurlength}}{a numeric vector of femur length}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Nosil, P. and B.J. Crespi. 2006. Experimental evidence that predation promotes divergence in adaptive radiation. \emph{Proceedings of the National Academy of Sciences, USA} 103: 9090-9095.
}

\references{
\url{http://www.sfu.ca/biology/faculty/crespi/pdfs/96-Nosil&Crespi2006PNAS.pdf}
}

\examples{
data(WalkingStickFemurs)
str(WalkingStickFemurs)

aovfit <- aov(femurlength ~ 1 + Error(specimen), data = WalkingStickFemurs)
aovfit
(aov.summary <- summary(aovfit))
MS.groups <- 0.002464
MS.error <- 0.000356
(F <- MS.groups / MS.error)
pf(F, 24, 25, lower.tail = FALSE)

# Among-group variance
(var.among <- (MS.groups - MS.error) / 2)

# Repeatability or Intraclass Correlation
var.among / (var.among + MS.error)

# Can use varcomps() and repeatability()
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps

# The same model can be fit with lme()
require(nlme)
lme.fit <- lme(femurlength ~ 1, random = ~ 1 | specimen,
  data = WalkingStickFemurs)
summary(lme.fit)
VarCorr(lme.fit)
R.lme <- repeatability(lme.fit)
R.lme
}

\keyword{datasets}
