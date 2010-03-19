\name{WalkingStickFemurs}
\alias{WalkingStickFemurs}
\docType{data}

\title{
WalkingStickFemurs
}

\description{
WalkingStickFemurs
}

\usage{data(WalkingStickFemurs)}
\format{
  A data frame with 50 observations on the following 2 variables.
  \describe{
    \item{\code{specimen}}{a numeric vector}
    \item{\code{femurlength}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Nosil, P. and B.J. Crespi. 2006. \emph{Proceedings of the National Academy of Sciences, USA} 103: 9090-9095.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(WalkingStickFemurs)
WalkingStickFemurs

# Convert specimen to factor
WalkingStickFemurs$specimen <- factor(WalkingStickFemurs$specimen)

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
repeatability(vc)

# The same model can be fit with lme()
require(nlme)
lme.fit <- lme(femurlength ~ 1, random = ~ 1 | specimen,
  data = WalkingStickFemurs)
summary(lme.fit)
VarCorr(lme.fit)
}

\keyword{datasets}
