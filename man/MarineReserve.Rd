\docType{data}
\name{MarineReserve}
\alias{MarineReserve}
\title{Marine Reserve Biomass}
\format{A data frame with 32 observations of one variable.  \describe{
\item{biomass.ratio}{} }}
\source{
  Halpern, B.S. 2003. The impact of marine reserves: do
  reserves work and does reserve size matter?
  \emph{Ecological Applications} 13: S117-S137.
}
\description{
  Relative biomass in 32 marine reserves.
}
\examples{
data(MarineReserve)
str(MarineReserve)
histogram(~biomass.ratio, MarineReserve)
}
\references{
  \url{http://www.nceas.ucsb.edu/~halpern/pdf/Halpern_EA_2003.pdf}
}
\keyword{datasets}

