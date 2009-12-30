\name{DietBreadthElVerde}
\alias{DietBreadthElVerde}
\docType{data}

\title{
Diet breadth in a rainforest community
}

\description{
Number of different species (\code{breadth}) in 127 species (\code{no.species}) in the rainforest community at El Verde, Puerto Rico
}

\usage{data(DietBreadthElVerde)}
\format{
  A data frame with 38 observations on the following 2 variables.
  \describe{
    \item{\code{breadth}}{a numeric vector}
    \item{\code{no.species}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Waide R.B. and W.B. Reagan, eds. 1996. \emph{The food web of a tropical rainforest}. University of Chicago Press, Chicago.}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(DietBreadthElVerde)
DietBreadthElVerde
sum(DietBreadthElVerde$no.species)
}

\keyword{datasets}
