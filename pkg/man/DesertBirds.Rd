\name{DesertBirds}
\alias{DesertBirds}
\docType{data}
\title{
Desert Bird Census Data
}

\description{
Census data for desert birds.
}

\usage{data(DesertBirds)}

\format{
  A data frame with 43 observations on the following 2 variables.
  \describe{
    \item{\code{species}}{a character vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Sauer, J.R., J.E. Hines, and J. Fallon. 2003. The North American breeding bird survey, results and analysis 1966-2002. Version 2003.1. USGS Patuxent Wildlife Research Center, Laurel, MD.
}

\references{
\url{http://www.mbr-pwrc.usgs.gov/bbs/}
}

\examples{
data(DesertBirds)
str(DesertBirds)
histogram( ~count, DesertBirds,
  xlab = "Abundance"
  )
}

\keyword{datasets}
