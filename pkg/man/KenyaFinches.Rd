\name{KenyaFinches}
\alias{KenyaFinches}
\docType{data}

\title{
Body Mass and Beak Length in Three Species of Finches in Kenya
}

\description{
Data on body mass and beak length in three species of finches: Crimson-rumped waxbill (\code{CRU.WAXB}), Cutthroat finch (\code{CUTTHROA}), and White-browed sparrow weaver (\code{WB.SPARW}).
}

\usage{data(KenyaFinches)}
\format{
  A data frame with 45 observations on the following 3 variables.
  \describe{
    \item{\code{species}}{a factor with levels: \code{CRU.WAXB}, \code{CUTTHROA}, and \code{WB.SPARW}}
    \item{\code{mass.g}}{a numeric vector}
    \item{\code{beaklength.mm}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Schluter, D. 1988. The evolution of finch communities on islands and continents: Kenya vs. Galapagos. \emph{Ecological Monographs} 58: 229-249.
}

\references{
\url{http://www.jstor.org/pss/1942538}
}

\examples{
data(KenyaFinches)
levels(KenyaFinches$species)

KenyaFinches
}

\keyword{datasets}
