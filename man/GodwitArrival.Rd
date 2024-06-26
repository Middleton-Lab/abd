% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{GodwitArrival}
\alias{GodwitArrival}
\title{Godwit Arrival Dates}
\format{
A data frame with 10 observations on the following 2 variables.
\describe{ \item{female}{a numeric vector} \item{male}{a
numeric vector} }
}
\source{
Gunnarsson, T.G., J.A. Gill, T. Sigurbjörnsson, and W.J. Sutherland.
2004. Pair bonds: arrival synchrony in migratory birds. \emph{Nature} 431:
646.
}
\usage{
GodwitArrival
}
\description{
Arrival dates for males and females in 10 pairs of Black-tailed Godwits
(\emph{Limosa limosa})
}
\examples{

xyplot(male~female, GodwitArrival, main='Arrival of Godwit pairs')

}
\references{
\url{http://en.wikipedia.org/wiki/Black-tailed_godwit}
}
\keyword{datasets}
