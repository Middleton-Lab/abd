\docType{data}
\name{LizardSprint}
\alias{LizardSprint}
\title{Sprint Speeds in Canyon Lizards}
\format{A data frame with 68 observations on the following 2 variables.
\describe{ \item{lizard}{a numeric vector} \item{speed}{a
numeric vector} }}
\source{
  \emph{inferred from} Huey, R.B. and A.E. Dunham. 1987.
  The repeatability of locomotor performance in natural
  populations of the lizard \emph{Sceloporus merriami}.
  \emph{Evolution} 42: 1116-1120.
}
\description{
  Sprint speeds (\code{speed}) in 34 canyon lizards
  (\emph{Sceloporous merriami}) measured in successive
  years in Big Bend National Park. Note that \code{lizard}
  is not coded as a factor.
}
\examples{
data(LizardSprint)
histogram(~speed, LizardSprint)
Lizard2 <- aggregate(speed ~ lizard, LizardSprint, mean)
histogram(~speed, Lizard2)
}
\references{
  \url{http://en.wikipedia.org/wiki/Sceloporus_merriami}
}
\keyword{datasets}

