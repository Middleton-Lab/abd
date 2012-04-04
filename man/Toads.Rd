\docType{data}
\name{Toads}
\alias{Toads}
\title{Right-handed Toads}
\format{A data frame with 19 observations on the following 2 variables.
\describe{ \item{n.toads}{a numeric vector} \item{prob}{a
numeric vector} }}
\description{
  Hypothetical probability of a toad being right-handed
}
\examples{
data(Toads)
Toads
# generate this data manually
cbind(0:18, dbinom(0:18, 18, 0.5))
xyplot(prob~n.toads, Toads, type = 'h', lwd = 4)
barchart(prob~n.toads, Toads, origin=0, horizontal=FALSE)
plotDist('binom', params = list(18,0.5), kind = 'hist')
}
\keyword{datasets}

