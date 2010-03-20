\name{OstrichTemp}
\alias{OstrichTemp}
\docType{data}

\title{
Ostrich Body and Brain Temperatures
}

\description{
Body and brain temperatures (C) in free-ranging ostriches (\emph{Struthio camelus}) at the the Lichtenburg Game Breeding Centre, Lichtenburg, South Africa.
}

\usage{data(OstrichTemp)}
\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{ostrich}}{a numeric vector identifying ostrich number}
    \item{\code{body.temp}}{a numeric vector of body temperature in degrees C}
    \item{\code{brain.temp}}{a numeric vector of brain temperature in degrees C}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Fuller, A., P.R. Kamerman, S.K. Maloney, G. Mitchell, and D. Mitchell. 2003. Variability in brain and arterial blood temperatures in free-ranging ostriches in their natural habitat. \emph{Journal of Experimental Biology} 206: 1171-1181.
}

\references{
\url{http://jeb.biologists.org/cgi/content/abstract/206/7/1171}
}

\examples{
data(OstrichTemp)
OstrichTemp
}

\keyword{datasets}
