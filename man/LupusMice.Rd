\docType{data}
\name{LupusMice}
\alias{LupusMice}
\title{Autoimmune Reactivity in Lupus-prone Mice}
\format{A data frame with 20 observations on the following 2 variables.
\describe{ \item{treatment}{a factor with levels: \code{enhanced},
\code{sham}, and \code{untreated}} \item{dilution}{a numeric vector
of the dilution level at which reactivity could be detected} }}
\source{
  McGaha, T.L., B. Sorrentino, and J.V. Ravetch. 2005.
  Restoration of tolerance in lupus by targeted inhibitory
  receptor expression. \emph{Science} 307: 590-593.
}
\description{
  Autoimmune reactivity (\code{dilution} at which
  reactivity could be detected) in three \code{treatment}s
  of lupus-prone mice.
}
\examples{
data(LupusMice)
str(LupusMice)
}
\keyword{datasets}

