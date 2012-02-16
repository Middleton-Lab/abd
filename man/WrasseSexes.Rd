\docType{data}
\name{WrasseSexes}
\alias{WrasseSexes}
\title{Distribution of Wrasses}
\format{A data frame with 3 observations on the following 3 variables.
\describe{ \item{males}{a numeric vector} \item{females}{a
numeric vector} \item{count}{a numeric vector} }}
\description{
  Number and sex of adult wrasses in a section of the Great
  Barrier Reef.
}
\examples{
data(WrasseSexes)
xtabs(count ~ males + females, WrasseSexes)
}
\keyword{datasets}

