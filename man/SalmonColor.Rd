\docType{data}
\name{SalmonColor}
\alias{SalmonColor}
\title{Pacific Salmon Color}
\format{A data frame with 35 observations on the following 2 variables.
\describe{ \item{species}{a factor with levels \code{kokanee} and
\code{sockeye}} \item{skin.color}{a numeric vector} }}
\source{
  Craig, J.K. and C. Foote. 2001. Countergradient variation
  and secondary sexual color: phenotypic convergence
  promotes genetic divergence in carotenoid use between
  sympatric anadromous and nonanadromous morphs of sockeye
  salmon (\emph{Oncorhynchus nerka}). \emph{Evolution} 55:
  380-391.
}
\description{
  Skin color sockeye and kokanee morphs of the Pacific
  salmon (\emph{Oncorhynchus nerka}) raised in a low
  carotenoid environment.
}
\examples{
data(SalmonColor)
SalmonColor
histogram(~skin.color | species, SalmonColor)
bwplot(skin.color ~ species, SalmonColor)
}
\keyword{datasets}

