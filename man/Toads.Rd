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
xyplot(prob~n.toads, Toads, type='h', lwd=4)
barchart(prob~n.toads, Toads, horizontal=FALSE)
}
\keyword{datasets}

