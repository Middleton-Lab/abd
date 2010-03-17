\name{ElectricFish}
\alias{ElectricFish}
\docType{data}

\title{
Electric Fish
}

\description{
Species abundance of electric fish upstream and downstream of the entrance of a tributary in the Amazon basin.
}

\usage{data(ElectricFish)}
\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{Tributary}}{a character vector}
    \item{\code{Species.upstream}}{a numeric vector of the number of species of electric fish present upstream of the tributary}
    \item{\code{Species.downstream}}{a numeric vector of the number of species of electric fish present downstream of the tributary}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Fernandes, C.C., J. Podos, and J.G. Lundberg. 2004. Amazonian ecology: tributaries enhance the diversity of electric fishes. \emph{Science} 305: 1960-1962.
}

\references{
\url{http://www.sciencemag.org/cgi/content/abstract/305/5692/1960}
}

\examples{
data(ElectricFish)
ElectricFish
}

\keyword{datasets}
