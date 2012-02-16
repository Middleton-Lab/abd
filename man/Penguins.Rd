\docType{data}
\name{Penguins}
\alias{Penguins}
\title{Penguin Heart Rate}
\format{A data frame with 24 observations on the following 2 variables.
\describe{ \item{group}{a factor with levels \code{BF}, \code{BM},
and \code{MF}} \item{slope}{a numeric vector} }}
\source{
  Green, J. A., P. J. Butler, A. J. Woakes, I. L. Boyd and
  R. L. Holder. 2001. Heart rate and rate of oxygen
  consumption of exercising macaroni penguins.
  \emph{Journal of Experimental Biology} 204: 673–684.
}
\description{
  Slope of regressions of mass-specific metabolic rate on
  heart rate for three groups of Macaroni Penguins.
}
\examples{
data(Penguins)
str(Penguins)
dotplot(slope~group, Penguins)
}
\keyword{datasets}

