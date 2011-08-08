\name{FlycatcherPatch}
\alias{FlycatcherPatch}
\docType{data}

\title{
Forehead Patch Size in Collared Flycatachers
}

\description{
Forehead patch size in 30 male Collared Flycatachers measured in two consecutive years.
}

\usage{data(FlycatcherPatch)}
\format{
  A data frame with 30 observations on the following 2 variables.
  \describe{
    \item{\code{patch98}}{a numeric vector}
    \item{\code{patch99}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Griffith, S.C. and B.C. Sheldon. 2001. Phenotypic plasticity in the expression of a sexually selected trait: neglected components of variation. \emph{Animal Behaviour} 61: 987-993.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(FlycatcherPatch)
str(FlycatcherPatch)
xyplot(patch99 ~ patch98, FlycatcherPatch)
}

\keyword{datasets}
