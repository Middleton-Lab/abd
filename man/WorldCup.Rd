\docType{data}
\name{WorldCup}
\alias{WorldCup}
\title{World Cup Goals}
\format{A data frame with 7 observations on the following 2 variables.
\describe{ \item{score}{a numeric vector} \item{count}{a
numeric vector} }}
\description{
  Number of goals per team during the 2002 World Cup.
}
\examples{
data(WorldCup)
xyplot(count ~ score, WorldCup, type='h', lwd=4)
barchart(count ~ score, WorldCup, origin=0, horizontal=FALSE)
}
\keyword{datasets}

