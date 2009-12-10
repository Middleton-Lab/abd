\name{DesertBirds}
\alias{DesertBirds}
\docType{data}
\title{
Desert bird census data
}

\description{
Desert bird census data
}

\usage{data(DesertBirds)}

\format{
  A data frame with 43 observations on the following 2 variables.
  \describe{
    \item{\code{Species}}{a character vector}
    \item{\code{Count}}{a numeric vector}
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

hist(DesertBirds$Count, breaks = 12,
	ylab = 'Frequency (Number of Species)',
	xlab = 'Abundance',
	main = '', col = 'red')

}

\keyword{datasets}
