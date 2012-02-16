\docType{data}
\name{Sumo}
\alias{Sumo}
\title{Sumo Wrestling Wins}
\format{A data frame with 16 observations on the following 2 variables.
\describe{ \item{wins}{a numeric vector} \item{count}{a
numeric vector} }}
\source{
  Duggan, M. and S.D. Leavitt. 2002. Winning isn't
  everything: Corruption in sumo wrestling. \emph{The
  American Economic Review} 92: 1594-1605.
}
\description{
  Counts of number of wins for sumo wrestlers.
}
\examples{
data(Sumo)
xyplot(count ~ wins, Sumo, type='h', lwd=4)
}
\keyword{datasets}

