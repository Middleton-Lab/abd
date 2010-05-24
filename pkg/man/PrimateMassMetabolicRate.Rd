\name{PrimateMassMetabolicRate}
\alias{PrimateMassMetabolicRate}
\docType{data}

\title{
Primate Metabolic Rates
}

\description{
Body mass (g) and metabolic rate (watts) for 17 species of primates.
}

\usage{data(PrimateMassMetabolicRate)}
\format{
  A data frame with 17 observations on the following 2 variables.
  \describe{
    \item{\code{mass.g}}{a numeric vector}
    \item{\code{bmr.watts}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Heusner, A.A. 1991. Size and power in mammals. \emph{Journal of Experimental Biology} 160: 25-54.
}

\references{
\url{http://jeb.biologists.org/cgi/content/abstract/160/1/25}
}

\examples{
data(PrimateMassMetabolicRate)
str(PrimateMassMetabolicRate)
PrimateMassMetabolicRate
}

\keyword{datasets}
