\docType{data}
\name{Rattlesnakes}
\alias{Rattlesnakes}
\title{Temperature Change and Meal Size in Rattlesnakes}
\format{A data frame with 17 observations on the following 2 variables.
\describe{ \item{meal.size}{a numeric vector}
\item{temp.change}{a numeric vector} }}
\source{
  Tattersall, G.J., W.K. Milsom, A.S. Abe, S.P. Brito, and
  D.V. Andrade. 2004. The thermogenesis of digestion in
  rattlesnakes. \emph{Journal of Experimental Biology} 207:
  579-585.
}
\description{
  Temperature change after a meal (% of body mass) in 17
  South American rattlesnakes (\emph{Crotalus durissus}).
}
\examples{
data(Rattlesnakes)
str(Rattlesnakes)
xyplot(meal.size ~ temp.change, Rattlesnakes)
}
\references{
  \url{http://jeb.biologists.org/cgi/content/abstract/207/4/579}
}
\keyword{datasets}

