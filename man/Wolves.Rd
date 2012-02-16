\docType{data}
\name{Wolves}
\alias{Wolves}
\title{Inbreeding in Wolves}
\format{A data frame with 24 observations on the following 2 variables.
\describe{ \item{inbreeding.coefficient}{a numeric vector}
\item{pups}{a numeric vector} }}
\source{
  Liberg, O.H., H. Andrén, H.-C. Pedersen, H. Sand, D.
  Sejberg, P. Wabakken, M. Åkesson, and S. Bensch. 2005.
  Severe inbreeding depression in a wild wolf (\emph{Canis
  lupus}) population. \emph{Biology Letters} 1: 17-20.
}
\description{
  Inbreeding coefficient and the number of pups produced in
  24 mated pairs of wolves (\emph{Canis lupus}) from
  1983-2002.
}
\examples{
data(Wolves)
Wolves
xyplot( inbreeding.coefficient ~ jitter(pups, amount=0.15), Wolves)
}
\keyword{datasets}

