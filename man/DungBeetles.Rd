% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{DungBeetles}
\alias{DungBeetles}
\title{Heritability of Body Condition in Dung Beetles}
\format{
A data frame with 36 observations on the following 2 variables.
\describe{ \item{id}{a numeric vector}
\item{offspring.condition}{a numeric vector} }
}
\source{
\emph{inferred from} Kotiaho, J.S., L.W. Simmons, and J.L. Tomkins.
2001. Towards a resolution of the lek paradox. \emph{Nature} 410: 684-686.
}
\usage{
DungBeetles
}
\description{
Body condition (\code{offspring.condition}) in 36 dung beetles
(\emph{Onthophagus taurus}) from 12 \code{male}s each mated to 3 different
virgin females.
}
\examples{

str(DungBeetles)
xyplot(offspring.condition ~ factor(id), DungBeetles,
  xlab='Dung Beetle',
  ylab='offspring condition')

}
\references{
\url{http://en.wikipedia.org/wiki/Dung_beetle}

\url{http://www.nature.com/nature/journal/v410/n6829/abs/410684a0.html}
}
\keyword{datasets}
