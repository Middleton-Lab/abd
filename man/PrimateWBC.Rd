\docType{data}
\name{PrimateWBC}
\alias{PrimateWBC}
\title{Primate White Blood Cell Counts and Promiscuity}
\format{A data frame with 9 observations on the following 2 variables.
\describe{ \item{WBC.less}{a numeric vector}
\item{WBC.more}{a numeric vector} }}
\source{
  Nunn, C.L., J.L. Gittleman, and J. Antonovics. 2000.
  Promiscuity and the primate immune system. \emph{Science}
  290: 1168-1170.
}
\description{
  White blood cell (WBC) counts in pairs of closely related
  primate species
}
\examples{
data(PrimateWBC)
xyplot(WBC.more ~ WBC.less, PrimateWBC)
}
\keyword{datasets}

