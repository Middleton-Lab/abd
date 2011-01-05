\name{HybridPollenSterility}
\alias{HybridPollenSterility}
\docType{data}

\title{
Sterility in Hybrid Pollens
}

\description{
Genetic distance between pairs of species of the genus \emph{Silene} and proportion of their hybrid offspring that are sterile.
}

\usage{data(HybridPollenSterility)}
\format{
  A data frame with 23 observations on the following 2 variables.
  \describe{
    \item{\code{genetic.distance}}{a numeric vector}
    \item{\code{proportion.sterile}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Moyle, L.C., M.S. Olson, and P. Tiffin. 2004. Patterns of reproductive isolation in three angiosperm genera. \emph{Evolution} 58: 1195-1208.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(HybridPollenSterility)
str(HybridPollenSterility)
xyplot(proportion.sterile ~ genetic.distance, HybridPollenSterility)
}

\keyword{datasets}
