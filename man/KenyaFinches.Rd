\docType{data}
\name{KenyaFinches}
\alias{KenyaFinches}
\title{Body Mass and Beak Length in Three Species of Finches in Kenya}
\format{A data frame with 45 observations on the following 3 variables.
\describe{ \item{species}{a factor with levels: \code{CRU.WAXB},
\code{CUTTHROA}, and \code{WB.SPARW}} \item{mass}{mass (g)}
\item{beak.length}{beak length (mm)} }}
\source{
  Schluter, D. 1988. The evolution of finch communities on
  islands and continents: Kenya vs. Galapagos.
  \emph{Ecological Monographs} 58: 229-249.
}
\description{
  Data on body mass and beak length in three species of
  finches: Crimson-rumped waxbill (\code{CRU.WAXB}),
  Cutthroat finch (\code{CUTTHROA}), and White-browed
  sparrow weaver (\code{WB.SPARW}).
}
\examples{
data(KenyaFinches)
table(KenyaFinches$species)
xyplot(beak.length ~ species, KenyaFinches)
bwplot(beak.length ~ species, KenyaFinches)
}
\references{
  \url{http://www.jstor.org/pss/1942538}
}
\keyword{datasets}

