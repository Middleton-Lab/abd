\docType{data}
\name{Tobacco2}
\alias{Tobacco2}
\title{Flower Length in Tobacco Plants}
\format{A data frame with 617 observations on the following 2 variables.
\describe{ \item{flower.length}{a numeric vector}
\item{generation}{a factor with levels \code{F1} \code{F2}} }}
\source{
  East, E.M. 1916. Studies on size inheritance in
  \emph{Nicotiana}. \emph{Genetics} 1: 164-176.
}
\description{
  Distribution of flow lengths in F1 and F2 populations of
  \emph{Nicotiana}.
}
\examples{
data(Tobacco2)
xtabs(~flower.length + generation, Tobacco2)
bwplot(flower.length ~ generation, Tobacco2)
}
\references{
  \url{http://www.genetics.org/content/vol1/issue2/}

  \url{http://en.wikipedia.org/wiki/Nicotiana}
}
\seealso{
  \code{\link{Tobacco}}
}
\keyword{datasets}

