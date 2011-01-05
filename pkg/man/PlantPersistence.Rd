\name{PlantPersistence}
\alias{PlantPersistence}
\docType{data}

\title{
Population Persistence Times
}

\description{
Persistence times (\code{generations}) in the annual plant \emph{Cardamine pensylvanica} in four experimental populations (\code{treatment}).
}

\usage{data(PlantPersistence)}

\format{
  A data frame with 16 observations on the following 2 variables.
  \describe{
    \item{\code{generations}}{a numeric vector}
    \item{\code{treatment}}{a factor with levels: \code{Isolated}, \code{Medium}, \code{Long}, and \code{Continuous}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Molofsky, J. and J.-B. Ferdy. 2005. Extinction dynamics in experimental metapopulations. \emph{Proceedings of the National Academy of Sciences (USA)} 102: 3726-3731.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(PlantPersistence)
xyplot(generations~treatment, PlantPersistence)
}

\keyword{datasets}
