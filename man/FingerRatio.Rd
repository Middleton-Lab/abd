\docType{data}
\name{FingerRatio}
\alias{FingerRatio}
\title{2D:4D Finger Ratio}
\format{A data frame with 46 observations on the following 2 variables.
\describe{ \item{CAGrepeats}{a numeric vector of the number of CAG
repeats} \item{finger.ratio}{a numeric vector of the ratio of digit
2 to digit 4} }}
\source{
  \emph{inferred from} Manning, J.T., P.E. Bundred, D.J.
  Newton, and B.F. Flanagan. 2003. The second to fourth
  digit ratio and variation in the androgen receptor gene.
  \emph{Evolution and Human Behavior} 24: 399-405.
}
\description{
  The ratio of the lengths of the index finger to the ring
  finger in 46 males and the number of CAG repeats for
  each.
}
\examples{
data(FingerRatio)
str(FingerRatio)
xyplot(finger.ratio ~ CAGrepeats, FingerRatio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio"
)
}
\references{
  \url{http://en.wikipedia.org/wiki/Digit_ratio}
}
\keyword{datasets}

