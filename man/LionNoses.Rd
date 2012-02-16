\docType{data}
\name{LionNoses}
\alias{LionNoses}
\title{Lion Age and Nose Coloration}
\format{A data frame with 32 observations on the following 2 variables.
\describe{ \item{age}{a numeric vector}
\item{proportion.black}{a numeric vector} }}
\source{
  Whitman, K., A.M. Starfield, H.S. Quadling and C. Packer.
  2004. Sustainable trophy hunting of African lions.
  \emph{Nature} 428: 175-178.
}
\description{
  Ages (in years; \code{age}) of 32 male lions and relative
  coloration of their noses (\code{proportion.black}).
}
\examples{
data(LionNoses)
xyplot(age ~ proportion.black, LionNoses)
}
\references{
  \url{http://www.nature.com/nature/journal/v428/n6979/abs/nature02395.html}
}
\keyword{datasets}

