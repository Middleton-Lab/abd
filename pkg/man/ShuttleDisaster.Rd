\name{ShuttleDisaster}
\alias{ShuttleDisaster}
\docType{data}

\title{
ShuttleDisaster
}

\description{
Data on temperature and number of O-ring failures for 23 space shuttle launches.
}

\usage{data(ShuttleDisaster)}
\format{
  A data frame with 23 observations on the following 2 variables.
  \describe{
    \item{\code{temperature}}{a numeric vector}
    \item{\code{failures}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Dalal, S.R., E.B. Fowlkes, and B. Hoadley. 1989. Risk Analysis of the Space Shuttle: Pre-Challenger Prediction of Failure. \emph{Journal of the American Statistical Association} 408: 945-957. 
}

\references{
Tufte, E.R. 1997. \emph{Visual Explanations: Images and Quantities, Evidence and Narrative}. Graphics Press.
}

\examples{
data(ShuttleDisaster)
str(ShuttleDisaster)
ShuttleDisaster
}

\keyword{datasets}
