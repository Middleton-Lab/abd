\docType{data}
\name{Lynx}
\alias{Lynx}
\title{Population Cycles of Lynx in Canada 1752-1819}
\format{A data frame with 68 observations on the following 2 variables.
\describe{ \item{year}{a numeric vector} \item{pelts}{a
numeric vector} }}
\source{
  Elton, C. and M. Nicholson. 1942. The ten-year cycle in
  numbers of the lynx in Canada. \emph{Journal of Animal
  Ecology} 11: 215-244.
}
\description{
  Number of lynx pelts (\code{pelts}) reported in Canada
  per year from 1752 to 1819.
}
\examples{
data(Lynx)
xyplot(pelts ~ year, Lynx, type=c('p','l'))
}
\keyword{datasets}

