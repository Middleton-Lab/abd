\name{LizardSprint}
\alias{LizardSprint}
\docType{data}

\title{
Sprint Speeds in Canyon Lizards
}

\description{
Sprint speeds in 34 canyon lizards (\emph{Sceloporous merriami}) measured in successive years in Big Bend National Park.
}

\usage{data(LizardSprint)}

\format{
  A data frame with 68 observations on the following 2 variables.
  \describe{
    \item{\code{lizard}}{a numeric vector}
    \item{\code{speed}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Huey, R.B. and A.E. Dunham. 1987. The repeatability of locomotor performance in natural populations of the lizard \emph{Sceloporus merriami}. \emph{Evolution} 42: 1116-1120.
}

\references{
\url{http://en.wikipedia.org/wiki/Sceloporus_merriami}
}

\examples{
data(LizardSprint)
histogram(~speed, LizardSprint)
Lizard2 <- aggregate(speed ~ lizard, LizardSprint, mean)
histogram(~speed, Lizard2)
}

\keyword{datasets}
