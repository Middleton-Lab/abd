\docType{data}
\name{ElectricFish}
\alias{ElectricFish}
\title{Electric Fish}
\format{A data frame with 12 observations on the following 3 variables.
\describe{ \item{tributary}{a character vector}
\item{species.upstream}{a numeric vector of the number of species of
electric fish present upstream of the tributary}
\item{species.downstream}{a numeric vector of the number of species
of electric fish present downstream of the tributary} }}
\source{
  Fernandes, C.C., J. Podos, and J.G. Lundberg. 2004.
  Amazonian ecology: tributaries enhance the diversity of
  electric fishes. \emph{Science} 305: 1960-1962.
}
\description{
  Species abundance of electric fish upstream and
  downstream of the entrance of a tributary in the Amazon
  basin.
}
\examples{
data(ElectricFish)
ElectricFish
require(grid)
xyplot(species.upstream ~ species.downstream, data = ElectricFish,
  panel=function(x, y, ...){
    grid.text(ElectricFish$tributary, x=x, y=y,
      rot = 45,
      gp = gpar(cex=.6),
      default.units = 'native')
    }
  )
}
\references{
  \url{http://www.sciencemag.org/cgi/content/abstract/305/5692/1960}
}
\keyword{datasets}

