\docType{data}
\name{LizardBite}
\alias{LizardBite}
\title{Bite Force in Collard Lizards}
\format{A data frame with 11 observations on the following 2 variables.
\describe{ \item{bite}{force of bite (N)}
\item{territory}{area of territory} }}
\source{
  Lappin, A. K., and J. F. Husak. 2005. Weapon performance,
  not size, determines mating success and potential
  reproductive output in the collared lizard
  (\emph{Crotaphytus collaris}). \emph{The American
  Naturalist} 166: 426-436.
}
\description{
  Bite force (N) and territory area in 11 male collared
  lizards (\emph{Crotaphytus collaris}).
}
\examples{
data(LizardBite)
str(LizardBite)
xyplot(territory ~ bite, LizardBite)
}
\keyword{datasets}

