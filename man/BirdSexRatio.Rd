\docType{data}
\name{BirdSexRatio}
\alias{BirdSexRatio}
\title{Sex Ratios in Birds}
\format{A data frame with 15 observations of one variable \describe{
\item{corr.coeff}{correlation coefficient of sex ratio in bird
offspring} }}
\source{
  West, S.A. and B.C. Sheldon. 2002. Constraints in the
  evolution of sex ratio adjustment. \emph{Science} 295:
  1695-1688.
}
\description{
  Correlation coefficient of sex ratio in bird offspring.
}
\examples{
data(BirdSexRatio)
histogram(~corr.coeff, BirdSexRatio, n = 10,
  xlab = "Correlation Coefficient")
}
\keyword{datasets}

