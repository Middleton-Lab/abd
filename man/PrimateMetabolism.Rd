\docType{data}
\name{PrimateMetabolism}
\alias{PrimateMetabolism}
\title{Primate Metabolic Rates}
\format{A data frame with 17 observations on the following 2 variables.
\describe{ \item{mass}{mass (g) } \item{bmr}{metabolic rate
(watts)} }}
\source{
  Heusner, A.A. 1991. Size and power in mammals.
  \emph{Journal of Experimental Biology} 160: 25-54.
}
\description{
  Body mass (g) and metabolic rate (watts) for 17 species
  of primates.
}
\examples{
data(PrimateMetabolism)
str(PrimateMetabolism)
xyplot(bmr ~ mass, PrimateMetabolism)
xyplot(bmr ~ mass, PrimateMetabolism, scales=list(log=TRUE))
}
\references{
  \url{http://jeb.biologists.org/cgi/content/abstract/160/1/25}
}
\keyword{datasets}

