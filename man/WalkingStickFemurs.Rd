\docType{data}
\name{WalkingStickFemurs}
\alias{WalkingStickFemurs}
\title{Walking Stick Femur Length}
\format{A data frame with 50 observations on the following 2 variables.
\describe{ \item{specimen}{a factor with levels \code{1-25}}
\item{femur.length}{a numeric vector of femur length} }}
\source{
  Nosil, P. and B.J. Crespi. 2006. Experimental evidence
  that predation promotes divergence in adaptive radiation.
  \emph{Proceedings of the National Academy of Sciences
  (USA)} 103: 9090-9095.
}
\description{
  Two measures of femur length for each of 25 walking
  sticks (\emph{Timema cristinae}).
}
\examples{
data(WalkingStickFemurs)
str(WalkingStickFemurs)

aov(femur.length ~ specimen, data = WalkingStickFemurs)
}
\references{
  \url{http://www.sfu.ca/biology/faculty/crespi/pdfs/96-Nosil&Crespi2006PNAS.pdf}
}
\keyword{datasets}

