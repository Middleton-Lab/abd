\name{VampireBites}
\alias{VampireBites}
\docType{data}

\title{
Vampire Bat Bites
}

\description{
Numbers of cattle bitten by the cow's estrous cycle.
}

\usage{data(VampireBites)}
\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{estrous}}{a factor with levels \code{no} \code{yes}}
    \item{\code{bitten}}{a factor with levels \code{no} \code{yes}}
    \item{\code{frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Turner, D.C. 1975. \emph{The vampire bat: a field study in behavior and ecology}. Johns Hopkins Press: Baltimore, MD.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(VampireBites)
VampireBites

fisher.test(matrix(VampireBites$Frequency, ncol = 2))

# With G-test
# Source from http://www.psych.ualberta.ca/~phurd/cruft/
source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r")
g.test(matrix(VampireBites$Frequency, ncol = 2))
g.test(matrix(VampireBites$Frequency, ncol = 2))$expected
}

\keyword{datasets}
