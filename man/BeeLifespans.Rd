\docType{data}
\name{BeeLifespans}
\alias{BeeLifespans}
\title{Bee Lifespans}
\format{A data frame with 33 observations on the following variable.
\describe{ \item{hours}{a numeric vector} }}
\source{
  \emph{inferred from} Visscher, P.K. and R. Dukas. 1997.
  Survivorship of foraging honey bees. \emph{Insectes
  Sociaux} 44: 1-5.
}
\description{
  Lifespan of 33 foraging honey bees.
}
\examples{
data(BeeLifespans)
histogram(~hours, BeeLifespans, n=10)
densityplot(~hours, BeeLifespans)
}
\keyword{datasets}

