\name{Hurricanes}
\alias{Hurricanes}
\docType{data}

\title{
Hurricane Intensities
}

\description{
Number of hurricanes greater than or equal to category 3 on the Safir-Simpson scale during the 20th century.
}

\usage{data(Hurricanes)}
\format{
  A data frame with 4 observations on the following 2 variables.
  \describe{
    \item{\code{Number.of.hurricanes}}{a numeric vector}
    \item{\code{Frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Blake, E.S., E.N. Rappaport, J.D. Jarrell, and C.W. Landsea. 2005. The deadliest, costliest, and most intense United States tropical cyclones from 1851 to 2006 (and other frequently requested hurricane facts). \emph{NOAA Technical Memorandum NWS TPC-4}.
}

\references{
\url{http://www.aoml.noaa.gov/hrd/Landsea/Blakeetal_noaamemoApr2007.pdf}
}

\examples{
data(Hurricanes)
Hurricanes
}

\keyword{datasets}
