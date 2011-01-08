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
    \item{\code{estrous}}{a factor with levels: \code{no} and \code{yes}}
    \item{\code{bitten}}{a factor with levels: \code{no} and \code{yes}}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Turner, D.C. 1975. \emph{The Vampire Bat: a Field Study in Behavior and Ecology}. Johns Hopkins Press: Baltimore, MD.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(VampireBites)
VampireBites

xtabs(count ~ estrous + bitten, data = VampireBites)
fisher.test(xtabs(count ~ estrous + bitten, data = VampireBites))

# With G-test
# Source from http://www.psych.ualberta.ca/~phurd/cruft/
try({
	source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r");
	g.test(xtabs(count ~ estrous + bitten, data = VampireBites));
	g.test(xtabs(count ~ estrous + bitten, data = VampireBites))$expected
	}
)
}

\keyword{datasets}
