\name{Earthworms}
\alias{Earthworms}
\docType{data}

\title{
Earthworm Diversity and Soil Nitrogen Levels
}

\description{
Number of earthworm species and total nitrogen content in the soil in 39 hardwood forest plots.}

\usage{data(Earthworms)}
\format{
  A data frame with 39 observations on the following 2 variables.
  \describe{
    \item{\code{worm.species}}{a numeric vector}
    \item{\code{nitrogen}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Gundale, M.J., W.M. Jolly, and T.H. Deluca. 2005. Susceptibility of a northern hardwood forest to exotic earthworm invasion. \emph{Conservation Biology} 19: 1075-1083.
}

\references{
\url{http://www3.interscience.wiley.com/journal/118701215/abstract}
}

\examples{
data(Earthworms)
str(Earthworms)
xyplot(nitrogen ~ worm.species, Earthworms)

}

\keyword{datasets}
