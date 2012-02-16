\docType{data}
\name{BodyFatHeatLoss}
\alias{BodyFatHeatLoss}
\title{Heat Loss and Body Fat}
\format{A data frame with 12 observations on the following 2 variables.
\describe{ \item{leanness}{a numeric vector}
\item{lossrate}{a numeric vector} }}
\source{
  Sloan, R.E.G. and W.R. Keatinge. 1973. Cooling rates of
  young people swimming in cold water. \emph{Journal of
  Applied Physiology} 35: 371-375.
}
\description{
  Heat loss during exercise and relative body fat in 12
  boys.
}
\examples{
data(BodyFatHeatLoss)
xyplot(lossrate ~ leanness, BodyFatHeatLoss)
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/pubmed/4732330}
}
\keyword{datasets}

