\docType{data}
\name{Pollen}
\alias{Pollen}
\title{Sterility in Hybrid Pollens}
\format{A data frame with 23 observations on the following 2 variables.
\describe{ \item{genetic.distance}{a numeric vector}
\item{proportion.sterile}{a numeric vector} }}
\source{
  Moyle, L.C., M.S. Olson, and P. Tiffin. 2004. Patterns of
  reproductive isolation in three angiosperm genera.
  \emph{Evolution} 58: 1195-1208.
}
\description{
  Genetic distance between pairs of species of the genus
  \emph{Silene} and proportion of their hybrid offspring
  that are sterile.
}
\examples{
data(Pollen)
str(Pollen)
xyplot(proportion.sterile ~ genetic.distance, Pollen)
}
\keyword{datasets}

