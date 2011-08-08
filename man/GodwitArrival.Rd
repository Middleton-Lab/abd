\name{GodwitArrival}
\alias{GodwitArrival}
\docType{data}

\title{
Godwit Arrival Dates
}

\description{
Arrival dates for males and females in 10 pairs of Black-tailed Godwits (\emph{Limosa limosa})
}

\usage{data(GodwitArrival)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{\code{female}}{a numeric vector}
    \item{\code{male}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Gunnarsson, T.G., J.A. Gill, T. Sigurbjörnsson, and W.J. Sutherland. 2004. Pair bonds: arrival synchrony in migratory birds. \emph{Nature} 431: 646.
}

\references{
\url{http://en.wikipedia.org/wiki/Black-tailed_godwit}
}

\examples{
data(GodwitArrival)
xyplot(male~female, GodwitArrival, main='Arrival of Godwit pairs')
}

\keyword{datasets}
