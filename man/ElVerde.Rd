\docType{data}
\name{ElVerde}
\alias{ElVerde}
\title{Diet Breadth in a Rainforest Community}
\format{A data frame with 38 observations on the following 2 variables.
\describe{ \item{breadth}{a numeric vector}
\item{num.species}{a numeric vector} }}
\source{
  Waide R.B. and W.B. Reagan, eds. 1996. \emph{The Food Web
  of a Tropical Rainforest}. University of Chicago Press,
  Chicago.
}
\description{
  Number of different species (\code{breadth}) in 127
  species (\code{no.species}) in the rainforest community
  at El Verde, Puerto Rico
}
\examples{
data(ElVerde)
ElVerde
xyplot(num.species ~ breadth, ElVerde, type='h',lwd=3)
}
\keyword{datasets}

