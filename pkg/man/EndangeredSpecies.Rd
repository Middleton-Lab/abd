\name{EndangeredSpecies}
\alias{EndangeredSpecies}
\docType{data}

\title{
Endangered and Threatened Species
}

\description{
Frequency of taxon groups on the U.S. Fish and Wildlife Service list of endangered and threatened species (2002).
}
\usage{data(EndangeredSpecies)}

\format{
  A data frame with 11 observations on the following 2 variables.
  \describe{
    \item{\code{taxon}}{a character vector}
    \item{\code{num.species}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
U.S. Fish and Wildlife Service. 2001. Number of U.S. listed species per calendar year.
}

\references{
\url{http://www.fws.gov/endangered/wildlife.html}
}

\examples{
data(EndangeredSpecies)
str(EndangeredSpecies)
EndangeredSpecies
}

\keyword{datasets}
