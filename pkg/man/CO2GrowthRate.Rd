\name{CO2GrowthRate}
\alias{CO2GrowthRate}

\docType{data}

\title{
Carbon Dioxide and Growth Rate
}

\description{
Growth rates of the unicellular alga \emph{Chlamydomonas} after 1,000 generations of selection under \code{High} and \code{Normal} levels of carbon dioxide.
}

\usage{data(CO2GrowthRate)}
\format{
  A data frame with 14 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{Normal} and \code{High}}
    \item{\code{growthrate}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Collins S. and G. Bell. 2004. Phenotypic consequences of 1,000 generations of selection at elevated CO2 in a green alga. \emph{Nature} 431: 566-569.
}

\references{
\url{http://www.nature.com/nature/journal/v431/n7008/full/nature02945.html}
}

\examples{
data(CO2GrowthRate)
CO2GrowthRate
}

\keyword{datasets}
